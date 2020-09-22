##############################################################################
##############################################################################
######                                                                  ######
######            reopen mapping project -- run SEIR                    ######
######                                                                  ######
##############################################################################
##############################################################################



### wrapper function for running SIR under one policy combination 
runSir <- function(place, policy, par, simRef){
  #number of phases
  np<-length(policy)
  
  sim0<-NA
  
  # for each phase
  for (j in 1:np){
    contact <- as.vector(policy[[j]])
    # load contact matrix from synthetic population
    Cmat<<-loadData(place, contact)
    
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
plotSIR2MSA <- function(fn, var, var_ref, small, legendMSA) {
  
  # state to plot
  name2plot<-c("Deaths","Cases","EmpLoss")
  
  # different y axis range, depending on MSAs and run version
  colvec<-c(6,7,1)
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
    
    ### plot two regions
    plot(0:max(TTT),var[i,]
         ,type="l",ylab="Percent of population",xlab="",ylim=c(0,yub[i]),lwd=2,col=mycol[colvec[i]])
    lines(0:max(TTT),var_ref[i,],
          type="l",lty=2,lwd=2,col=mycol[colvec[i]])
    
    legend("topleft",
           legend=legendMSA,
           col=c(mycol[colvec[i]],mycol[colvec[i]]),
           lty=1:2, 
           horiz=F,lwd=1.3,bty="n",cex=1.3)
    
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
msaList<<-c("5600","1600","6920","3760")

#NP throughout
nopolicy <- "_W4-S3-N3-B2-R2-P2-F2-E2-M1"
policy<-rep(nopolicy, 4)



#### foreach MSA
m_num<-1
for (m in msaList){
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
  var<-runSir(m, policy, par, NA)
    
  # reference MSA
  if (m_num %% 2 ==1){
    var_ref<-var
  }
  
  # compare results for two MSAs
  if (m_num %% 2 ==0){
    if (m_num<3){
      legendMSA <- c("Chicago", "NYC")
    }else{
      legendMSA <- c("Kansas City", "Sacramento")
    }
    
    fnEnd <- paste('_', m, '_', msaList[m_num-1], verTag, ".png", sep="")
    plotSIR2MSA(fnEnd, var, var_ref, floor(m_num/2)-1,legendMSA)
  }
  m_num<-m_num+1
}
rm(Cmat, par)








