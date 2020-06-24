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


##########################################################################
######## main function for grid search 
##########################################################################
gridSearch <- function(place, covid){
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(paste("!! Starting grid search for MSA", place))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  
  #number of time periods
  nt<-dim(covid)[1]
  
  ### global parameter values
  vpar <-vparameters0
  vpar1<-vpar
  vpar2<-vpar
  
  ### contact matrix with no policy
  Cmat1<-loadData(place, refPhase1)
  eig<-largestEigenvalue(Cmat1)
  
  ### essential only contact matrix
  Cmat2<-loadData(place, refPhase2)
  
  
  ### grid search input parmaeter
  gsPar <- checkLoad(gsParm)
  gsPar <- as.vector(gsPar[gsPar$msa==place,])
  gsCol <- names(gsPar)
  
  ### lower/upper bound for beta1, beta2 and initial condition
  lb   <-unlist(gsPar[grep("lb",   gsCol, perl=T)])
  ub   <-unlist(gsPar[grep("ub",   gsCol, perl=T)])
  step1<-unlist(gsPar[grep("step", gsCol, perl=T)])
  
  ### start/end time of each phase: phase 1, phase 2 with beta1, phase 2 with beta2
  TTTcali <-c(0,TTT[2],gsPar$T_beta2,nt-1)
  
  ### time range of the sample used for estimation
  tRange<-seq(gsPar$T_range_start, rangeToT)

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
      vpar1["beta"] <-parm$beta1
      vpar2["beta"] <-parm$beta2
      
      ### initial condition
      initNumIperType<<-parm$I0
      
      ### for three periods phase 0-old beta, phase 1-old beta, phase 1-new beta
      sim0<-NA
      for (t in 1:3){
        inits<-initialCondition(TTTcali[t],sim0)
        vt <- seq(0,diff(TTTcali)[t],1) 
        
        # set contact matrix and parameter in each period
        if (t==1){
          Cmat <<-Cmat1
          vpar_j<-vpar1
        }else if (t==2){
          Cmat <<-Cmat2
          vpar_j<-vpar1
        }else{
          Cmat <<-Cmat2
          vpar_j<-vpar2        
        }
        
        # RUN SIR
        sim_j = as.data.frame(lsoda(inits, vt, SEIIRRD_model, vpar_j))
        
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
    
    dead<-log(covid$deathper100k)
    case<-log(covid$caseper100k)
    if (ng>1){
      ## fit log death, min squared loss
      err<-log(fitDeath) - matrix(1,ng,1)%*%dead
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
    par(mfrow=c(1,1))
    plot(covid$t, log(fitCase[gstar,]), type="l", ylim=c(-4,10))
    lines(covid$t, case, type="l", lwd="2", col="blue")
    lines(covid$t, dead, type="l", lwd="2", col="red")
    lines(covid$t, log(fitDeath[gstar,]),type="l")
 
    abline(v=TTTcali[2], col="gray")
    
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
  
  
  #export csv  
  parmOut<-matrix(c(parm$beta1,parm$beta2,parm$I0),1,3)
  colnames(parmOut)<-c("beta1","beta2","I0") 
  fn <- paste(parmPath, paste(caliParm, place, datv, ".csv", sep=""), sep="/")
  write.table(parmOut, file=fn, sep=",",col.names=TRUE,row.names=FALSE)
  print(paste("export calibrated parameters :",fn))
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

rm(COVID, covid, Cmat, gridSearch)
