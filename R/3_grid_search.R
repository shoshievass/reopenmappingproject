##############################################################################
##############################################################################
######                                                                  ######
######           reopen mapping project -- grid search                  ######
######                                                                  ######
##############################################################################
##############################################################################


#####################################
# Grid Search
#####################################

rangeToT<-60
place<-"5600"

### load covid death and cases data
COVID <-read.csv("covid_case_death_jh.csv",header=TRUE)
covid<-COVID[COVID$msa==place,c("t","deathper100k","caseper100k")]
nt<-dim(covid)[1]

### regular parameters

vpar<-vparameters0
vpar1<-vpar
vpar2<-vpar

contact<-"_regular"
Cmat1<-loadData(place, contact)
eig<-largestEigenvalue(Cmat1)

### social distance parameters
contact<-"_socialdistance"
Cmat2<-loadData(place, contact)


### grid search
# back out beta from the R0 and tau
if (place=="5600"){
  lb<-c(20, -1,  5)
  ub<-c(50,  2,  20)
  step1<-c(5, 1,  5)
  step<-rbind(step1, step1*0.1, step1*0.01, step1*0.001)
  #" R0=25.55 beta=0.0154 I0=3.034 lnI0=1.11 R02=19.45 beta2=0.01172" change beta at 30
  #" R0=29.95 beta=0.01859 I0=2.435 lnI0=0.89 R02=15.55 beta2=0.00965" psi=0.7, change beta at 25
  # lb<-c(29.95, 0.89,  15.55)
  # ub<-c(29.95, 0.89,  15.55)
  name  <-"NYC"
  TTT <<-c(0,15,25,nt-1)
  tRange<-seq(10,rangeToT)
}
if (place=="1600"){
  lb<-c(20, -5,  0)
  ub<-c(60,  2,  40)
  step1<-c(10, 1,  10)
  step<-rbind(step1, step1*0.1, step1*0.01, step1*0.001)
  #" R0=38.9 beta=0.00319 I0=0.1212 lnI0=-2.11 R02=21.1 beta2=0.00173" replica, change beta at 20
  #" R0=40.9 beta=0.003354 I0=0.1212 lnI0=-2.11 R02=18.9 beta2=0.00155" replica, change beta at 25
  # 
  # lb<-c(40.9, -2.11,  18.9)
  # ub<-c(40.9, -2.11,  18.9)
  name  <-"Chicago"
  TTT <<-c(0,15,25,nt-1)
  tRange<-seq(13,rangeToT)
}
if (place=="6920"){
  lb<-c(0, -1,  0)
  ub<-c(10,  2,  5)
  step1<-c(2, 1,  1)
  step<-rbind(step1, step1*0.1, step1*0.01, step1*0.001)
  #" R0=6.02 beta=0.001285 I0=1.116 lnI0=0.11 R02=1.78 beta2=0.0003801" replica 2 beta at vector psi, break at t=20
  #" R0=5.98 beta=0.001275 I0=1.116 lnI0=0.11 R02=0.89 beta2=0.0001897" replica 2 beta at psi=0.7, break at t=25
  #" R0=6.22 beta=0.001326 I0=1 lnI0=0 R02=1.89 beta2=0.0004029" replica 2 beta at psi=0.7, break at t=23, set initial condition at 1
  
  # lb<-c(6.22, 0,  1.89)
  # ub<-c(6.22, 0,  1.89)
  name  <-"Sacramento"
  
  TTT <<-c(0,15,23,nt-1)
  tRange<-seq(5,rangeToT)
}
runOK<-0


# j rounds of grid search
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
  if (place=="6920"){
    #for sacramento, we fix initial value
    gList<-expand.grid(rList,tList)
    gList<-cbind(gList[,1], 0, gList[,2])
  }else{
    gList<-expand.grid(rList,iList,tList)
  }
  ng<-dim(gList)[1]
  
  fitDeath<-matrix(0,ng,nt)
  fitCase <-matrix(0,ng,nt)
  
  start_time <- Sys.time()
  for (i in 1:ng){
    ### parameters
    beta <-gList[i,1]/(infectDuration*eig) 
    beta2<-gList[i,3]/(infectDuration*eig)
    vpar1["beta"] <-beta
    vpar2["beta"] <-beta2

    
    ### initial condition
    I0<-exp(gList[i,2])
    initNumIperType<<-I0
    
    ### for three periods phase 0-old beta, phase 1-old beta, phase 1-new beta
    sim0<-NA
    for (t in 1:3){
      inits<-initialCondition2(nc,N,TTT[t],sim0)
      vt <- seq(0,diff(TTT)[t],1) 
      
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
        sim0<-rbind(sim0[1:TTT[t],], sim_j)
      }else{
        sim0<-sim_j
      }      
    }
    
    #store death and cases (per 100k)
    fitDeath[i,]<-extractState("D",sim0)*1e3
    fitCase[i,] <-extractSeveralState(c("Ihc","Rq","Rqd","D"),sim0)*1e3
    if ((i %% 10)==0){
      print(paste(i,"/",ng,
                " beta=",format(beta,digits=4),
                " beta2=",format(beta2,digits=4),
                " I0=",   format(I0,digits=4),
                " cum death=",format(fitDeath[i,nt],digits=0), sep=""))
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
  
  g0<-as.double(gList[gstar,])
  beta <-g0[1]/(infectDuration*eig) 
  beta2<-g0[3]/(infectDuration*eig) 
  names(beta)="beta"
  
  
  ### plot comparison
  par(mfrow=c(1,1))
  
  plot(covid$t, dead, type="l", lwd="2", col="red")
  lines(covid$t, log(fitDeath[gstar,]),type="l")
  abline(v=TTT[2], col="gray")
  
  ## check within grid search boundary
  stopifnot(min((g0<ub) * (g0>lb))==1)
  
  runOK<-ifelse(j==4,1,0)
}


print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("grid search both rounds done!!")
print(paste(
            " R0=",   format(g0[1],digits=4),
            " beta=", format(beta, digits=4),
            " I0=",   format(exp(g0[2]),digits=4),
            " lnI0=", format(g0[2],digits=4),
            " R02=",  format(g0[3],digits=4),
            " beta2=",format(beta2,digits=4),
            sep=""))
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

# 
# if (max(ub-lb)==0 | runOK==1){
#   #fit death and case
#   fn <- paste('calibrate_R0_I0_beta_', place, ver, ".pdf", sep="")
#   pdf(paste(outPath, "figure", fn, sep="/"))
#   plot(covid$t, log(fitCase[gstar,]),
#        xlab="", ylab="log death and cases per 100k", type="l",ylim=c(-4,10),cex.lab=psize, cex.axis=psize)
#   lines(covid$t, dead, type="l", col="red")
#   lines(covid$t, case, type="l", col="blue")
#   lines(covid$t, log(fitDeath[gstar,]),type="l")
#   abline(v=TTT[2], col="gray")
#   abline(v=max(tRange), col="orange")
#   abline(v=min(tRange), col="orange")
#   text(max(covid$t)*0.6, min(dead[tRange])+0.2*(max(dead[tRange])-min(dead[tRange])),
#        paste(name,
#              ": beta0(%)=",format(round(beta*10000)/100,digits=3),
#              ", beta1(%)=",format(round(beta2*10000)/100,digits=3),
#              ", I0=",format(exp(g0[2]),digits=3),
#              sep=""))
#   dev.off()
# 
#   #simpler version
#   fn <- paste('calibrate_beta_I0_', place, ver, ".pdf", sep="")
#   pdf(paste(outPath, "figure", fn, sep="/"))  
#   plot(covid$t, log(fitDeath[gstar,]),
#        xlab="", ylab="Log Death Per 100 000 of Population", lwd=1.5, 
#        xaxt = "n", type="l",col="red", lty=5, ylim=c(-4,6), cex.lab=psize)
#   axis(side  = 1, at = seq(0,90,15))
#   lines(covid$t, dead, type="l", lwd=1.5, col="red")
#   abline(v=TTT[2], col="gray", lwd=1.5)
#   abline(v=max(tRange))
#   abline(v=min(tRange))
#   dev.off()
# }







