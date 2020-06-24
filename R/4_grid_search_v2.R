##############################################################################
##############################################################################
######                                                                  ######
######           reopen mapping project -- grid search                  ######
######                                                                  ######
##############################################################################
##############################################################################


## a formal estimation procedure is currently under development 

#parse calibrated parameter
gridPar <- function(parm, R0scale){
  
  #transmission rate
  beta1<-parm[[1]]/R0scale
  beta2<-parm[[3]]/R0scale
  
  ### initial condition
  I0<-exp(parm[[2]])
  return(list(I0=I0,beta1=beta1,beta2=beta2))
}

# plot comparison of cases and deaths -----------------------------------------------------
plotCali <- function(t, fitCase, fitDeath, case, dead, t2) {
  par(mfrow=c(1,1))
  
  plot(t, fitCase, 
       xlab="", ylab="log death and cases per 100k", 
       type="l", col="blue", lty=5, ylim=c(-4,10),cex.lab=psize, cex.axis=psize)
  lines(t, case, type="l", lwd="2", col="blue")
  lines(t, dead, type="l", lwd="2", col="red")
  lines(t, fitDeath,type="l",col="red", lty=5)
  
  abline(v=t2[1],  col="gray")
  abline(v=t2[2], col="black")  
  
}

### grid search 
grid_search <- function(place, covid, rangeToT){
  
  
  place<-m
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
  
  
  ### grid search
  # back out beta from the R0 and tau
  if (place=="5600"){
    lb<-c(5, -7, 0)
    ub<-c(30, -1, 25)
    step1<-c(5, 1,  5)
    name  <-"NYC"
    TTTcali <-c(0,17,nt-1)
    tRange<-seq(5,rangeToT)
  }
  if (place=="1600"){
    lb<-c(20, -9,  4)
    ub<-c(50, -4, 16)
    step1<-c(5, 1,  4)
    name  <-"Chicago"
    TTTcali <-c(0,15,nt-1)
    tRangeD<-seq(13,rangeToT) #to be generalized to be first death + 2
    tRangeC<-seq(offset+2,rangeToT)
  }
  if (place=="6920"){
    lb<-c( 0, -5,  0)
    ub<-c(10,  0,  5)
    step1<-c(2, 1,  1)
    name  <-"Sacramento"
    TTTcali <-c(0,14,nt-1)
    tRangeD<-seq(offset+2,rangeToT) #to be generalized to be first death + 2
    tRangeC<-seq(offset+2,rangeToT)
  }
  
  # j rounds of grid search
  step<-rbind(step1, step1*0.1, step1*0.01, step1*0.001)
  for (j in 1:4){
    if (j==1){
      #coarser grid
      rList<-seq(lb[1],ub[1],step[1,1])
      iList<-seq(lb[2],ub[2],step[1,2])
      tList<-seq(lb[3],ub[3],step[1,3])
    }else{
      #finer grid
      rList<-seq(max(lb[1],g0[1]-step[j-1,1]),min(ub[1],g0[1]+step[j-1,1]),step[j-1,1])
      iList<-seq(max(lb[2],g0[2]-step[j-1,2]),min(ub[2],g0[2]+step[j-1,2]),step[j-1,2]) 
      tList<-seq(max(lb[3],g0[3]-step[j-1,3]),min(ub[3],g0[3]+step[j-1,3]),step[j-1,3])
    }
    gList<-expand.grid(rList,iList,tList)
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
      for (t in 1:2){
        inits<-initialCondition(TTTcali[t],sim0)
        vt <- seq(0,diff(TTTcali)[t],1) 
        
        if (t==1){
          Cmat <<-Cmat1
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
      
      #death (per 100k) in data
      fitDeath[i,]<-log(extractState("D",sim0)*1e3)
      fitCase[i,] <-log(extractSeveralState(c("Ihc","Rq","Rqd","D"),sim0)*1e3)
      if ((i %% 10)==0){
        print(paste(i,"/",ng,
                    " beta1=",format(parm$beta1,digits=4),
                    " beta2=",format(parm$beta2,digits=4),
                    " I0=",   format(parm$I0   ,digits=4), sep=""))
      }
    }
    
    end_time <- Sys.time()
    print(end_time - start_time)
    
    #death and cases in the data
    dead<-log(covid$deathper100k)
    case<-log(covid$caseper100k)
    if (ng>1){
      ## fit log death, min squared loss
      errD<-fitDeath - matrix(1,ng,1)%*%dead
      errC<-fitCase  - matrix(1,ng,1)%*%case
      sse<-rowSums(errC[,tRangeC]^2) + 5 * rowSums(errD[,tRangeD]^2)
      
      #best fit
      gstar<-which.min(sse)
    }else{
      gstar<-1
    }
    
    #optimized parameter
    g0<-as.double(gList[gstar,])
    parm<-gridPar(g0, infectDuration*eig)
    
    ### plot comparison of cases and deaths
    
    plotCali(covid$t,fitCase[gstar,],fitDeath[gstar,], case, dead, c(TTTcali[2], rangeToT)) 
    

    
    
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
              " I0=",    format(parm$I0,digits=4),
              sep=""))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  
  
  #export csv  
  # parmOut<-matrix(c(parm$beta1,parm$beta2,parm$I0),1,3)
  # colnames(parmOut)<-c("beta1","beta2","I0") 
  # fn <- paste(parmPath, paste(caliParm, place, datv, "_v2.csv", sep=""), sep="/")
  # write.table(parmOut, file=fn, sep=",",col.names=TRUE,row.names=FALSE)
  # print(paste("export calibrated parameters :",fn))
}


#### time range of sample used for estimation
rangeToT<-60

### load covid death and cases data
COVID <-checkLoad(deathData)

#### foreach MSA run grid search to calibrate parameters
for (m in msaList){
  
  ### death data for this MSA
  covid<-COVID[COVID$msa==m,c("t","deathper100k","caseper100k")]
  
  offset<-0
  nt<-dim(covid0)[1]
   
  if (offset>0){
    covid$caseper100k[(offset+1):nt] <- 10 * covid0$caseper100k[1:(nt-offset)]
    covid$caseper100k[1:offset]<-0
  }
  
  ### estimate parameter
  grid_search(m, covid, rangeToT)
}



# rm(rangeToT, COVID, covid, Cmat, grid_search)
