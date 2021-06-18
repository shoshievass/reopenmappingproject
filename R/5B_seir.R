##############################################################################
##############################################################################
######                                                                  ######
######            reopen mapping project -- run SEIR                    ######
######        compare contact and composition across MSAs               ######
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
    vpar["beta"]<-setBeta(contact, par) 
    
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
      # what2plot<-list(D=c("D"),C=c("Ihc","Rq","Rqd","D"),Emp=c("Ihc","Rq","D"), Ia=c("hospital"))
      what2plot<-list(D=c("D"))
      var <-matrix(0, length(what2plot), max(TTT)+1)
      for (i in 1:length(what2plot)){
        var[i,] <- extractSeveralState(what2plot[[i]],sim0)
        if (i==4){
          ### number of people in hospital, assuming each inpatient stay for 10 days
          hosp <- cumsum(var[i,])
          var[i,] <- c(hosp[1:HOSPDAYS], 
                       hosp[(HOSPDAYS+1):length(hosp)]-hosp[1:(length(hosp)-HOSPDAYS)])
        }
      }
    }
  }
  return(var)
}


# compare SIR for two MSAs -----------------------------------------------------
plotSIR2MSA <- function(fn, var, small, legendMSA) {
  
  # state to plot
  # name2plot<-c("Deaths","Cases","EmpLoss","Hospital")
  # legendpos<-c("bottomright", "bottomright", "topright", "topright")
  
  name2plot<-c("Deaths")
  legendpos<-c("bottomright")
  
  
  # different y axis range, depending on MSAs and run version
  colvec <-c(6,6,6,6)
  colvec2<-c(2,2,2,2)
  
  
  if (Generic==2){
    yub<-c(1,80,60,1.2)
  }else{
    if (small==1){
      yub<-c(0.06,20,10,0.1)
    }else{
      yub<-c(0.8,70,50,0.8)    
    }    
  }
  l<-length(name2plot)
  
  for (i in 1:l){
    
    if (fn!="")  fn1 <- file.path(outPath, "figure", paste(name2plot[i], fn, sep=""))
    if (fn!="")  pdf(fn1, family="Palatino")
    par(mar=c(2,5,1,1))
    
    color  <- mycol[colvec[i]]
    color2 <- mycol[colvec2[i]]

    ### plot two regions
    plot(0:max(TTT),var[i,]
         ,type="l",ylab="Percent of population",xlab="",ylim=c(0,yub[i]),lwd=2,col=color,cex.lab=psize, cex.axis=0.7*psize)
    lines(0:max(TTT),var[l+i,],
          type="l",lty=2,lwd=2,col=mycol[colvec[i]])
    lines(0:max(TTT),var[l*2+i,],
          type="l",lty=5,lwd=2,col=mycol[colvec2[i]])
    lines(0:max(TTT),var[l*3+i,],
          type="l",lty=1,lwd=2,col=mycol[colvec2[i]])
    
    legend(legendpos[i],
           legend=legendMSA,
           col=c(color,color,color2,color2),
           lty=c(1,2,5,1), 
           horiz=F,lwd=2, bty="n",cex=0.6*psize)
    
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
if (Generic<2){
  msaListCntct<<-c("5600","1600","5600","1600",
                   "6920","3760","6920","3760")
  msaListCmpst<<-c("5600","5600","1600","1600",
                   "6920","6920","3760","3760")
}else{
  # msaListCntct<<-c("5600","1600","5600","1600",
  #                  "6920","3760","6920","3760",
  #                  "5600","3760","5600","3760",
  #                  "1600","3760","1600","3760",
  #                  "5600","6920","5600","6920")
  # msaListCmpst<<-c("5600","5600","1600","1600",
  #                  "6920","6920","3760","3760",
  #                  "5600","5600","3760","3760",
  #                  "1600","1600","3760","3760",
  #                  "5600","5600","6920","6920")
  msaListCntct<<-c("5600","1600","5600","1600",
                   "6920","3760","6920","3760")
  msaListCmpst<<-c("5600","5600","1600","1600",
                   "6920","6920","3760","3760")
}

#NP throughout but change in beta
policy<-c("_W4-S3-N3-B2-R2-P2-F2-E2-M1",
          "_W4-S3-N3-B2-R2-P2-F2-E2-M2",
          "_W4-S3-N3-B2-R2-P2-F2-E2-M3",
          "_W4-S3-N3-B2-R2-P2-F2-E2-M4")

## WFH
# policy<-c("_W4-S3-N3-B2-R2-P2-F2-E2-M1", 
#           "_W2-S1-N1-B1-R2-P1-F1-E2-M2", 
#           "_W2-S1-N1-B1-R2-P1-F1-E2-M3", 
#           "_W2-S1-N1-B1-R2-P1-F1-E2-M4")

## AS
# policy<-c("_W4-S3-N3-B2-R2-P2-F2-E2-M1", 
#           "_W3-S2-N1-B2-R2-P2-F2-E2-M2", 
#           "_W3-S2-N1-B2-R2-P2-F2-E2-M3", 
#           "_W3-S2-N1-B2-R2-P2-F2-E2-M4")

## 60+
# policy<-c("_W4-S3-N3-B2-R2-P2-F2-E2-M1", 
#           "_W4-S3-N1-B2-R2-P2-F2-E1-M2", 
#           "_W4-S3-N1-B2-R2-P2-F2-E1-M3", 
#           "_W4-S3-N1-B2-R2-P2-F2-E1-M4")

verTagLocal <- paste(verTag, "NP", sep="_")


#### foreach MSA
m_num<-1
var_agg<-c()
for (m in msaListCntct){
  # load calibrated parameters for this location
  par   <-calibratedPar(m, Generic)

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
      legendMSA <- c("NYC Contact -\nNYC Population", "\nChicago Contact -\nNYC Population\n",
                     "NYC Contact -\nChicago Population", "\nChicago Contact -\nChicago Population\n")
      small<-0
    }else if (m_num<=8){
      legendMSA <- c("Sacramento Contact -\nSacramento Population", "\nKansas City Contact -\nSacramento Population\n",
                     "Sacramento Contact -\nKansas City Population", "\nKansas City Contact -\nKansas City Population\n")
      small<-1
    }else if (m_num<=12){
      legendMSA <- c("NYC Contact - NYC Population", "Kansas City Contact - NYC Population",
                     "NYC Contact - Kansas City Population", "Kansas City Contact - Kansas City Population")
      small<-0
    }else if (m_num<=16){
      legendMSA <- c("Chicago Contact - Chicago Population", "Kansas City Contact - Chicago Population",
                     "Chicago Contact - Kansas City Population", "Kansas City Contact - Kansas City Population")
      small<-0
    }else if (m_num<=20){
      legendMSA <- c("NYC Contact - NYC Population", "Sacramento Contact - NYC Population",
                     "NYC Contact - Sacramento Population", "Sacramento Contact - Sacramento Population")
      small<-0
    }else{
      stop("not coded up figure format")
    }
    
    fnEnd <- paste('_', m, '_', msaListCntct[m_num-1], verTagLocal, ".pdf", sep="")
    plotSIR2MSA(fnEnd, var_agg, small,legendMSA)
  }

  # reset aggregation
  if (m_num %% 4 ==0){
    var_agg<-c()
  }  
  m_num<-m_num+1
  
}
rm(Cmat, par)








