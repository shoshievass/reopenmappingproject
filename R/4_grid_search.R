##############################################################################
##############################################################################
######                                                                  ######
######           reopen mapping project -- grid search                  ######
######                                                                  ######
##############################################################################
##############################################################################


## a formal estimation procedure is currently under development 

### parse calibrated parameter
gridPar <- function(parm, R0scale){

  #transmission rate
  beta1<-parm[[1]]/R0scale
  beta2<-parm[[2]]/R0scale
  
  ### initial condition
  I0<-exp(parm[[3]])
  
  return(list(beta1=beta1,beta2=beta2,I0=I0))
}
  
### set up grid points
gridPoints <- function(lb, ub, step, j, g0){

  if (j==1){
    #initial round with coarse grid
    list1<-seq(lb[1],ub[1],step[1,1])
    list2<-seq(lb[2],ub[2],step[1,2])
    list3<-seq(lb[3],ub[3],step[1,3])
  }else{
    #finer grid
    list1<-seq(max(lb[1],g0[1]-step[j-1,1]),min(ub[1],g0[1]+step[j-1,1]),step[j-1,1])
    list2<-seq(max(lb[2],g0[2]-step[j-1,2]),min(ub[2],g0[2]+step[j-1,2]),step[j-1,2]) 
    list3<-seq(max(lb[3],g0[3]-step[j-1,3]),min(ub[3],g0[3]+step[j-1,3]),step[j-1,3])
  }
  return(expand.grid(list1,list2,list3))
}


### set up contact matrix and beta for each phase in the calibration
setCmatBeta <- function(gsPar, t, CmatList, betaList){
  
  #which version of contact matrix and beta to use
  cmatVer <- gsPar[[paste("cali",t,"_C",sep="")]]
  betaVer <- gsPar[[paste("cali",t,"_beta",sep="")]]
  
  
  ### edit global parameters
  vpar <- vparameters0
  vpar["beta"]<-betaList[betaVer]
  return(list(Cmat=CmatList[[cmatVer]], vpar=vpar))
}

##########################################################################
######## main function for grid search 
##########################################################################
gridSearch <- function(place, covid){
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(paste("!! Starting grid search for MSA", place))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

  ### grid search input parmaeter
  gsPar <- checkLoad(gsParm)
  gsPar <- as.vector(gsPar[gsPar$msa==place,])
  gsCol <- names(gsPar)

  #number of time periods
  T1<-gsPar$T1
  covid<-covid[covid$t>T1,]
  nt<-dim(covid)[1]
  
  ### timing
  ### start/end time of each phase: phase 1, phase 2 with beta1, phase 2 with beta2
  TTTcali <-c(0,gsPar$caliT2-T1,gsPar$caliT3-T1,nt-1)

  ### time range of the sample used for estimation
  tRange<-seq(gsPar$T_range_start, gsPar$T_range_end)-T1
  
  #show several lines in calibration plot
  tVertL<-c(gsPar$T2,gsPar$T3,min(tRange)+T1,max(tRange)+T1)
  
  ### fit log death or death
  logD <- gsPar$logDeath
  logTag <- ifelse(logD, "Log ", "")
  
  ### number of periods for predicted death
  #   = 0 if not using predicted death to evaluate out of sample fit
  TpredD <- gsPar$Toos_predDeath
  
  ### no policy, essential only, and cautious reopening contact matrix
  CmatList <- list(loadData(place, refPhase1), 
                   loadData(place, refPhase2), 
                   loadData(place, refPhase3))
  eig<-largestEigenvalue(CmatList[[1]])
  print(paste(" largest eigenvalue of contact matrices=", format(largestEigenvalue(CmatList[[1]]),digits=4),
              format(largestEigenvalue(CmatList[[2]]),digits=4),
              format(largestEigenvalue(CmatList[[3]]),digits=4)))
  
  
  ### lower/upper bound for beta1, beta2 and initial condition
  lb   <-unlist(gsPar[grep("lb",   gsCol, perl=T)])
  ub   <-unlist(gsPar[grep("ub",   gsCol, perl=T)])
  step1<-unlist(gsPar[grep("step", gsCol, perl=T)])
  
  ### j rounds of grid search with incrementally small steps
  step<-rbind(step1, step1*0.1, step1*0.01, step1*0.001)
  
  g0<-NA
  for (j in 1:4){
    ### grid
    gList<-gridPoints(lb, ub, step, j, g0)
    
    ### placeholder for simulated death/case
    ng<-dim(gList)[1]
    fitDeath<-matrix(0,ng,nt)
    fitCase <-matrix(0,ng,nt)
    
    start_time <- Sys.time()
    for (i in 1:ng){
      ### parameters
      parm<-gridPar(gList[i,], infectDuration*eig)
      betaList<-c(parm$beta1, parm$beta2)
      
      ### initial condition
      initNumIperType<<-parm$I0
      
      ### for three periods phase 0-old beta, phase 1-old beta, phase 1-new beta
      sim0<-NA
      for (t in 1:3){
        inits<-initialCondition(TTTcali[t],sim0)
        vt <- seq(0,diff(TTTcali)[t],1) 
        
        # set contact matrix and parameter in each period
        CB <- setCmatBeta(gsPar, t, CmatList, betaList)
        Cmat<<-CB$Cmat

        # RUN SIR
        sim_j = as.data.frame(lsoda(inits, vt, SEIIRRD_model, CB$vpar))
        
        if (t>1){
          sim0<-rbind(sim0[1:TTTcali[t],], sim_j)
        }else{
          sim0<-sim_j
        }      
      }
      
      ### simulated death and cases (per 100k) 
      fitDeath[i,]<-extractState("D",sim0)*1e3
      fitCase[i,] <-extractSeveralState(c("Ihc","Rq","Rqd","D"),sim0)*1e3
      if ((i %% 10)==0){
        print(paste(i,"/",ng,
                    " beta1=",format(parm$beta1,digits=4),
                    " beta2=",format(parm$beta2,digits=4),
                    " I0=",   format(parm$I0   ,digits=4), sep=""))
      }
    }
    
    end_time <- Sys.time()
    print(end_time - start_time)
    
    
    dead<-covid$deathper100k
    case<-covid$caseper100k
    
    ## log transform of deaths/cases
    if (logD==1){
      dead<-log(dead)
      case<-log(case)
      fitDeath<-log(fitDeath)
      fitCase <-log(fitCase)
    }
    if (ng>1){
      ## fit death, min squared loss
      err<-fitDeath - matrix(1,ng,1)%*%dead
      sse<-rowSums(err[,tRange]^2)
      
      #best fit
      gstar<-which.min(sse)
    }else{
      gstar<-1
    }
    
    #optimized parameter
    g0<-as.double(gList[gstar,])
    parm<-gridPar(g0, infectDuration*eig)
      
    ### plot comparison
    par(mfrow=c(1,2))
    plotCali(covid$t, dead, fitDeath[gstar,], 0, tVertL, logTag, TpredD)
    plotCali(covid$t, case,  fitCase[gstar,], 1, tVertL, logTag, 0)
    
    ## check within grid search boundary
    if(min((g0<ub) * (g0>lb))!=1){
      print(rbind(lb,g0,ub))
    }
    stopifnot(min((g0<ub) * (g0>lb))==1)
  }
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("grid search both rounds done!!")
  print(paste(" beta1=", format(parm$beta1,digits=4),
              " beta2=", format(parm$beta2,digits=4),
              " I0=",    format(parm$I0,digits=4),sep=""))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

  
  ### export csv
  parmOut<-matrix(c(parm$beta1,parm$beta2,parm$I0),1,3)
  colnames(parmOut)<-c("beta1","beta2","I0")
  fn <- paste(parmPath, paste(caliParm, place, datv, ".csv", sep=""), sep="/")
  write.table(parmOut, file=fn, sep=",",col.names=TRUE,row.names=FALSE)
  print(paste("export calibrated parameters :",fn))

  ### show calibration result
  par(mfrow=c(1,1))
  fn <- paste(outPath, "figure",
              paste('calibrate_beta_I0_pred_d_msa', place, ".pdf", sep=""), sep="/")
  pdf(fn)
  plotCali(covid$t, dead, fitDeath[gstar,], 0, tVertL, logTag, TpredD)
  dev.off()
  print(paste("  saved plot:",fn))
}


### load covid death and cases data
COVID <-checkLoad(deathData)

#### foreach MSA run grid search to calibrate parameters
for (m in msaList){
  
  ### death data for this MSA
  covid<-COVID[COVID$msa==m,c("t","deathper100k","caseper100k")]
  
  ### estimate parameter
  gridSearch(m, covid)
}

# rm(COVID, covid, Cmat, gridSearch)
