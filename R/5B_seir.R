##############################################################################
##############################################################################
######                                                                  ######
######            reopen mapping project -- run SEIR                    ######
######                                                                  ######
##############################################################################
##############################################################################




### wrapper function for running SIR under one policy combination 
runSir <- function(placeCnt, placePop, policy, par, simRef){
  #number of phases
  np<-length(policy)
  
  sim0<-NA
  
  # for each phase
  for (j in 1:np){
    contact <- as.vector(policy[[j]])
    # load contact matrix from synthetic population
    
    #### superimpose composition
    ## load composition
    CDataPop <-loadContactData(placePop,contact)
    
    ## load contact
    CDataCnt <-loadContactData(placeCnt,contact)

    types<-CDataCnt$types[,(colnames(CDataCnt$types) %in% c("n"))==FALSE]
    
    ## merge, keep all types in the contact matrix
    types <- merge(types, CDataPop$types, by=c("naics","age","sick","wfh","shift","work_poi","active_emp"),all.x = TRUE) 
   
    ## revert to original order, keep only merged types
    types <- types[order(types$ego.x), ]
    valid_types <- is.na(types$n)==FALSE
    types <- types[valid_types,]
    
    #redefine ego indicators
    types$ego <- 1:nrow(types)
    types<<-types
    
    globalTypeVectors(types)
    Cmat<<-CDataCnt$Cmat[valid_types,valid_types]
    stopifnot(mean(dim(Cmat))==nrow(types))
    
    # load default parameters and set transmission rate
    vpar<-vparameters0
    vpar["beta"]<-setBeta(contact, par, j) 
    
    ### initial condition
    inits<-initialCondition(TTT[j],sim0)
    vt <- seq(0,diff(TTT)[j],1) 
    
    # run SIR
    sim_j = as.data.frame(lsoda(inits, vt, SEIIRRD_model, vpar))
    
    # add tag to keep track of which type is working
    sim_j <- tagActiveEmp(sim_j)
    
    # aggregate results
    if (j>1){
      sim0<-rbind(sim0[1:TTT[j],], sim_j)
    }else{
      sim0<-sim_j
    }
    
    # check population static
    check<-rowSums(sim0[,c("S5","E5","Ia5","Ins5","Ihc5","Rq5","Rnq5","Rqd5","D5")])
    stopifnot(abs(min(check)-max(check))<0.1)
    
    # output compartments to plot
    if (j==np){
      state2plot<-c("Deaths","Cases","Employment Loss")
      what2plot<-list(D=c("D"),C=c("Ihc","Rq","Rqd","D"),Emp=c("Ihc","Rq","D"))
      var <-matrix(0, length(what2plot), max(TTT)+1)
      for (i in 1:length(what2plot)){
        var[i,] <- extractSeveralState(what2plot[[i]],sim0)
      }
    }
  }
  return(var)
}




# compare SIR for two MSAs -----------------------------------------------------
plotSIR2MSA <- function(fn, var, small, legendMSA) {
  
  # state to plot
  name2plot<-c("Deaths","Cases","EmpLoss")
  
  # different y axis range, depending on MSAs and run version
  colvec <-c(6,7,1)
  colvec2<-c(4,5,2)
  if (Generic==2){
    yub<-c(0.9,70,40)
  }else{
    if (small==1){
      yub<-c(0.001,0.2,0.1)
    }else{
      yub<-c(0.8,60,30)    
    }    
  }

  for (i in 1:length(colvec)){
    
    if (fn!="")  fn1 <- file.path(outPath, "figure", paste(name2plot[i], fn, sep=""))
    if (fn!="")  png(fn1)
    color  <- mycol[colvec[i]]
    color2 <- mycol[colvec2[i]]
    idx <- (i-1)*3
    ### plot two regions
    plot(0:max(TTT),var[i,]
         ,type="l",ylab="Percent of population",xlab="",ylim=c(0,yub[i]),lwd=2,col=color)
    lines(0:max(TTT),var[3+i,],
          type="l",lty=2,lwd=2,col=mycol[colvec[i]])
    lines(0:max(TTT),var[6+i,],
          type="l",lty=5,lwd=2,col=mycol[colvec2[i]])
    lines(0:max(TTT),var[9+i,],
          type="l",lty=1,lwd=2,col=mycol[colvec2[i]])
    
    legend("topleft",
           legend=legendMSA,
           col=c(color,color,color2,color2),
           lty=c(1,2,5,1), 
           horiz=F,lwd=1.3,bty="n",cex=1)
    
    if (fn!=""){
      dev.off()
      print(paste("  saved plot:",fn1))
    }
  }
  

}


#####################################
# run SIR
#####################################
start_time <- Sys.time()


### get max number of policy scenarios across MSAs
P <- checkLoad(policyParm) 
np <- length(grep("Scenario",colnames(P),perl=T))


# we need to run all 4 MSAs together, in the following order
msaListCntct<<-c("5600","1600","5600","1600","6920","3760","6920","3760")
msaListCmpst<<-c("5600","5600","1600","1600","6920","6920","3760","3760")


#NP throughout
nopolicy <- "_W4-S3-N3-B2-R2-P2-F2-E2-M1"
policy<-rep(nopolicy, 4)



#### foreach MSA
m_num<-1
var_agg<-c()
for (m in msaListCntct){
  # load calibrated parameters for this location
  par<-calibratedPar(m, Generic)
  
  ### msa policy and date
  msaPD <- loadPolicyDates(m)
  
  # number of phases for this MSA
  np <- length(msaPD$refPolicy)
  
  # timing of each policy
  TTT <<-c(msaPD$TVec[1:np]-msaPD$TVec[1], msaPD$TVec[np+1])
  
  # initial condition
  initNumIperType<<-par$I0

  # sir across phases, output time series of compartments that we want to plot
  var<-runSir(m, msaListCmpst[m_num], policy, par, NA)
    
  # consolidate outputs
  var_agg<-rbind(var_agg, var)

  
  # compare results for two MSAs
  if (m_num %% 4 ==0){
    if (m_num<=4){
      legendMSA <- c("NYC Contact - NYC Population", "Chicago Contact - NYC Population",
                     "NYC Contact - Chicago Population", "Chicago Contact - Chicago Population")
      small<-0
    }else{
      legendMSA <- c("Sacramento Contact - Sacramento Population", "Kansas City Contact - Sacramento Population",
                     "Sacramento Contact - Kansas City Population", "Kansas City Contact - Kansas City Population")
      small<-1
    }
    
    fnEnd <- paste('_', m, '_', msaListCntct[m_num-1], verTag, ".png", sep="")
    plotSIR2MSA(fnEnd, var_agg, small,legendMSA)
  }

  # reset aggregation
  if (m_num %% 4 ==0){
    var_agg<-c()
  }  
  m_num<-m_num+1
  
}
# rm(Cmat, par)








