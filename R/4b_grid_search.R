##############################################################################
##############################################################################
######                                                                  ######
######           reopen mapping project -- grid search                  ######
######                                                                  ######
##############################################################################
##############################################################################


## a formal estimation procedure is currently under development 

### parse calibrated parameter
gridPar <- function(parm, R0scale, phase){
  if (phase==1){
    ### initial condition
    I0<-exp(parm[[1]])
    if(I0<=0) stop("initial infected fraction needs to be >0!")
    
    #transmission rate
    beta0<-parm[[2]]/R0scale
    beta1<-beta0
  }else if (phase<np){
    I0<-initNumIperType
    beta0<-parm[[1]]/R0scale
    beta1<-beta0
  }else{
    I0<-initNumIperType
    beta0<-parm[[1]]/R0scale
    beta1<-parm[[2]]/R0scale    
  }

  if(min(beta0,beta1)<0) stop("transmission rate needs to be >=0!")
  return(list(beta0=beta0,beta1=beta1,I0=I0))
}

### set up grid points
gridPoints <- function(lb, ub, step, j, g0){
  
  npar<-length(lb)
  stopifnot(npar<=2)
  if (j==1){
    if (npar>1){
      list1<-seq(lb[1],ub[1],step[1,1])
      list2<-seq(lb[2],ub[2],step[1,2])
      grid<-expand.grid(list1,list2)
    }else{
      list1<-seq(lb[1],ub[1],step[1])
      grid<-expand.grid(list1)
    }
  }else{
    #finer grid
    if (npar>1){
      list1<-seq(max(lb[1],g0[1]-step[j-1,1]),min(ub[1],g0[1]+step[j-1,1]),step[j-1,1])
      list2<-seq(max(lb[2],g0[2]-step[j-1,2]),min(ub[2],g0[2]+step[j-1,2]),step[j-1,2]) 
      grid<-expand.grid(list1,list2)
    }else{
      list1<-seq(max(lb[1],g0[1]-step[j-1]),min(ub[1],g0[1]+step[j-1]),step[j-1])
      grid<-expand.grid(list1)
    }
  }
  return(grid)
}




### set up grid points
runSEIR <- function(g, R0scale, phase, t, CmatList, TTTcali, sim0){
  
  ### parameters
  parm    <-gridPar(g, R0scale, phase)
  
  ### initial condition
  initNumIperType<<-parm$I0
  
  inits<-initialCondition(TTTcali[t],sim0)
  vt <- seq(0,diff(TTTcali)[t],1) 
  
  # set contact matrix and parameter in each period
  Cmat<<-CmatList[[t]]
  
  vpar <- vparameters0
  vpar["beta0"] <- parm$beta0
  vpar["beta1"] <- parm$beta1
  
  # RUN SIR
  sim_j = as.data.frame(lsoda(inits, vt, SEIIRRD_model, vpar))
  return(sim_j)
}



### demean a time series, mostly for cases
# adjust for -Inf if we take log of 0
demeanSeries <- function(X0){

  if (is.vector(X0)){
    #remove infinity
    X<-X0[!is.infinite(X0)]
    X1<-X0[is.infinite(X0)]
    
    #demean
    X <- X - mean(X)
    X <- c(X1,X)
  }else{
    #remove infinity
    xmean <- colMeans(X0)
    X<-X0[,!is.infinite(xmean)]
    X1<-X0[,is.infinite(xmean)]
    
    #demean
    n <- dim(X)[2]
    X <-X - rowMeans(X)%*%matrix(1,1,n)
    X <- cbind(X1,X)
  }
  return(X)
}



# plot fit in caliration -----------------------------------------------------
plotCali <- function(fn, xdata, xfit, DC, tVertL) {
  
  if(length(xdata)!=length(xfit))  stop("data and predicted series do not have the same length!")
  
  nt<-length(xdata)
  t<-0:(nt-1)
  
  ### plot cases or deaths
  if (DC==0){
    ylbl <- "Death Per 100 000 of Population"
  }else{
    ylbl <- "Log Case Per 100 000 of Population"
  }
  
  ## do we export figure
  if (fn!="") pdf(fn)
  
  ### plot fit
  xdata_finite <- xdata[is.finite(xdata)]
  xfit_finite  <- xfit[is.finite(xfit)]
  plot(t, xfit, 
       ylim=c(min(0, min(xdata_finite),min(xfit_finite)), max(max(xdata_finite), max(xfit_finite))), 
       xlab="", ylab=ylbl, 
       lwd=1.5, xaxt = "n", type="l",col="red", lty=5, cex.lab=psize)
  # x-axis show dates
  t2show<-seq(0,max(t),10)
  axis(side  = 1, at = t2show, label=format(t2show+TNAUGHT, "%m/%d"))
  
  lines(t[1:nt], xdata[1:nt], type="l", lwd=1.5, col="red")
  
  ### indicate key dates
  nline<-length(tVertL)
  lcolors<-c(rep("gray",nline-2),rep("black",2))
  abline(v=tVertL,col=lcolors)
  
  if (fn!="") dev.off()
}


##########################################################################
######## main function for grid search 
##########################################################################
gridSearch <- function(m, covid){
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(paste("!! Starting grid search for MSA", m))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  
  ### grid search input parmaeter
  gsPar <- checkLoad(gsParm)
  gsPar <- as.vector(gsPar[gsPar$msa==m,])
  if (dim(gsPar)[1]==0) stop(paste("missing MSA ", m, " in grid search parameter file", sep=""))
  gsCol <- names(gsPar)
  
  ### msa policy and date
  msaPD <- loadPolicyDates(m)
  
  ### initial time
  T1<-msaPD$TVec[1]
  
  ### data
  covid<-covid[covid$t>T1,]
  nt<-dim(covid)[1] 
  if(nt<gsPar$T_range_end) stop("not enough data for calibration sample period")
  
  ## number of phases
  np<<-length(msaPD$refPolicy)
  
  ### contact matrices
  CmatList<-list()
  eigvList<-c()
  for (i in 1:np){
    if (msaPD$TVec[i]<nt){
      CmatList[[i]]<-loadData(m, msaPD$refPolicy[i])
      eigvList[i]  <-largestEigenvalue(CmatList[[i]])
      print(paste("policy in phase", i, ":",  msaPD$refPolicy[i]))
    }
  }
  R0scale<-infectDuration*eigvList[1]
  np<-length(CmatList)
  print(paste(" largest eigenvalue of contact matrices=", paste(format(eigvList,digits=4), collapse=", ")))  
  
  ### timing
  TTTcali <-msaPD$TVec-T1
  TTTcali[np+2]<-nt-1
  
  ### append one more round of contact matrix
  CmatList[[np+1]]<-CmatList[[np]]
  
  ### time range of the sample used for estimation
  tRange<-seq(gsPar$T_range_start, gsPar$T_range_end)-T1
  
  #show several lines in calibration plot
  tVertL<-c(TTTcali[2:np],min(tRange)+T1,max(tRange)+T1)
  
  ### death and cases in the data
  dead<-covid$deathper100k;  try(if(any(dead<0)) stop("Error in death data."))
  case<-covid$caseper100k;   try(if(any(case<0)) stop("Error in case data."))
  
  ### lower/upper bound for beta1, beta2 and initial condition
  lb   <-unlist(gsPar[grep("lb",   gsCol, perl=T)])
  ub   <-unlist(gsPar[grep("ub",   gsCol, perl=T)])
  step1<-unlist(gsPar[grep("step", gsCol, perl=T)])
  
  ### j rounds of grid search with incrementally small steps
  step<-rbind(step1, step1*0.1, step1*0.01, step1*0.001)
  
  
  ### we have different transmission beta parameter for each phase. 
  # in phase 1 there is also initial condition
  # we search for parameters of one phase at a time to reduce the dimension
  # when searching for parameters for phase t, we use both phase t and phase t+1 to evaluate fit
  # so we make sure the parameter gives us good fit out of sample
  par_star<-c()
  
  sim0<-NA
  for (t in 1:np){
    ### range of time in this phase
    ttt <- c(TTTcali[t]:TTTcali[t+2])+1
    fitRange <- max(min(tRange),min(ttt)):min(max(ttt),max(tRange))
    print(paste("Phase", t,"/",np,": estimation time range:", 
                min(fitRange), "to", max(fitRange)))

    ### select parameters in this phase
    if (t==1){
      par_select<-c(t:(t+1))
    }else if (t==np){
      par_select<-c((t+1):(t+2))
    }else{
      par_select<-c(t+1)
    }

    ## each round of grid search 
    start_time <- Sys.time()
    g0<-NA
    for (j in 1:dim(step)[1]){
      ### grid
      gList<-gridPoints(lb[par_select], ub[par_select], step[,par_select], j, g0)
      
      ### placeholder for simulated death/case
      ng<-dim(gList)[1]
      fitDeath<-matrix(0,ng,max(ttt))
      fitCase <-matrix(0,ng,max(ttt))
      
      ### for each parameter 
      for (i in 1:ng){
        sim_0_tt<-sim0
        for (tt in t:(t+1)){
          sim_j <- runSEIR(gList[i,], R0scale, t, tt, CmatList, TTTcali, sim_0_tt)
          if (tt>1){
            sim_j<-rbind(sim_0_tt[1:TTTcali[tt],], sim_j)
          }
          sim_0_tt<-sim_j
        }
        
        ### simulated death and cases (per 100k) 
        fitDeath[i,]<-extractState("D",sim_j)*1e3
        fitCase[i,] <-extractSeveralState(c("Ihc","Rq","Rqd","D"),sim_j)*1e3
        if ((i %% 10)==0){
          print(paste(i,"/",ng, sep=""))
        }
      }
      
      ## fit death, min squared loss
      err<-fitDeath - matrix(1,ng,1)%*%dead[1:max(ttt)]
      sse<-rowSums(err[,fitRange]^2) 
      
      ## also use cases in the last phase to min out of sample error
      if (t==np){
        ## error fitting demeaned log cases
        err_case<-demeanSeries(log(fitCase)) - matrix(1,ng,1)%*%demeanSeries(log(case))
        
        ## for each weight
        cw_list<-seq(0,200,25)
        oos_err<-matrix(0,length(cw_list),2)
        for (cw in 1:length(cw_list)){
          sse_both<-sse + cw_list[cw]*rowSums(err_case[,fitRange]^2)
          # selected parameter to minimize out death and cases
          oos_err[cw,2]<-which.min(sse_both)
          #out of sample fit
          oos_err[cw,1]<-sum(err[oos_err[cw,2],max(fitRange):nt]^2)
        }
        gstar<-oos_err[which.min(oos_err[,1]),2]
      }else{
        gstar<-which.min(sse)
      }
      
      #best fit
      g0   <-as.double(gList[gstar,])
      deadfit <- fitDeath[gstar,]
      casefit <-  fitCase[gstar,]
      
      ### plot comparison
      par(mfrow=c(1,2))
      plotCali("", dead[1:max(ttt)], deadfit, 0, tVertL)
      plotCali("", log(case[1:max(ttt)]), log(casefit), 1, tVertL)
      
      ## check within grid search boundary
      if(min((g0<ub[par_select]) * (g0>lb[par_select]))!=1) {
        print(rbind(lb[par_select],g0,ub[par_select]))
        stop(paste("grid search hit boundary, revise range of grid search in", gsParm))
      }
    }
    
    ### get output for the optimal parameter
    sim_j <- runSEIR(gList[gstar,], R0scale, t, t, CmatList, TTTcali, sim0)

    ### record optimized phase 1 to t and print calibrated parameter
    if (t>1){
      sim0<-rbind(sim0[1:TTTcali[t],], sim_j)
      if (t<np){
        print(paste(" beta_",t,"=", format(gList[gstar,],digits=4),sep=""))
      }else{
        print(paste(" beta_",t,"(young/healthy)=", format(gList[gstar,1],digits=4),
                    " beta_",t,"(old/sick)=", format(gList[gstar,2],digits=4),sep=""))
      }
    }else{
      sim0<-sim_j
      print(paste(" I0=", format(gList[gstar,1],digits=4),
                  " beta_",t,"=", format(gList[gstar,2],digits=4),sep=""))   
      parm<-gridPar(gList[gstar,], R0scale, t)
      
    }
    
    ## collect calibrated parameter (after transformation)
    parm<-gridPar(gList[gstar,], R0scale, t)
    if (t==1){
      par_star<-c(par_star, parm$I0, parm$beta0)
    }else if (t<np){
      par_star<-c(par_star, parm$beta0)
    }else{
      par_star<-c(par_star, parm$beta0, parm$beta1)
    }
    
    end_time <- Sys.time()
    print(end_time - start_time)
  }
  
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print("grid search done!!")
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  
  
  
  ### export calibration results as csv
  parmOut<-matrix(par_star,1,dim(step)[2])
  colnames(parmOut)<-c("I0",paste("beta",1:(dim(step)[2]-1),sep=""))
  checkWrite(file.path(calibratedParPath, paste(caliParm, m, ".csv", sep="")),
             parmOut, "calibrated parameters")
  
  
  ### plot calibration result
  fn <- file.path(outPath, "figure", paste("calibrate_beta_I0_death_msa", m, ".pdf", sep=""))
  plotCali(fn, dead, deadfit, 0, tVertL)
  print(paste("  saved plot:",fn))
  
  fn <- file.path(outPath, "figure", paste("calibrate_beta_I0_case_msa", m, ".pdf", sep=""))
  plotCali(fn, log(case), log(casefit), 1, tVertL)
  print(paste("  saved plot:",fn))
  
  par(mfrow=c(1,2))
  plotCali("", dead, deadfit, 0, tVertL)
  plotCali("", log(case), log(casefit), 1, tVertL)
  
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