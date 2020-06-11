##############################################################################
##############################################################################
######                                                                  ######
######                 COVID Project SIR                                ######
######                                                                  ######
##############################################################################
##############################################################################



rm(list = ls())
dir<-Sys.getenv("HOME")

if (dir=="/Users/hanyang") {
  proj <- paste(dir, "Dropbox", "COVID Project", sep="/")
  dataPath <- paste(proj, "stata", "data", sep="/")
  codePath <- paste(proj, "R", sep="/")
  outPath  <- paste(proj, "R", "output", sep="/")
} 
setwd(dataPath)

## set up global varibales and functions
source(paste(codePath,"002_programs.R"  ,sep='/'))
source(paste(codePath,"001_config.R"  ,sep='/'))


#####################################
# Grid Search
#####################################

FROMT<<-15
OPENT<<-0

rangeToT<-60

place<-"6920"



### load covid death and cases data
COVID <-read.csv("covid_case_death_jh.csv",header=TRUE)
covid<-COVID[COVID$msa==place,c("t","deathper100k","caseper100k")]
nt<-dim(covid)[1]
TT<-nt-FROMT


### regular parameters
contact<-"_regular"
Cmat1<-loadData(place, contact)
vpar <-vparameters0

### compute largest eigen value from contact matrix
eig<-largestEigenvalue(Cmat1)

### social distance parameters
contact<-"_socialdistance"
Cmat2<-loadData(place, contact)


### grid search
# back out beta from the R0 and tau
if (place=="5600"){
  lb<-c(15, 2)
  ub<-c(30,10)
  step1<-c(0.5,0.5)
  step2<-step1*0.1
  #" R0=28.4 I0=2.95 beta=0.01695"
  # selected
  lb<-c(28.4, 2.95)
  ub<-c(28.4, 2.95)
  name  <-"NYC"
  tRange<-seq(10,rangeToT)
}
if (place=="1600"){
  #starting condition FRED
  lb<-c(8, 0.1)
  ub<-c(20,  2)
  #starting condition Replica
  lb<-c(30, 0.02)
  ub<-c(45,  1.6)
  step1<-c(0.5,0.1)
  step2<-step1*0.1
  #" R0=17.7 I0=0.52 beta=0.01085" FRED
  #" R0=35   I0=0.2 beta=0.002817" Replica
  # selected
  lb<-c(35, 0.2)
  ub<-c(35, 0.2)
  name  <-"Chicago"
  tRange<-seq(15,rangeToT)
}
if (place=="6920"){
  lb<-c(2, 1)
  ub<-c(10,14)
  step1<-c(0.5,0.5)
  step2<-step1*0.1
  #" R0=4.05 I0=2.55 beta=0.003302"
  #" R0=3.75 I0=2.2 beta=0.0007885" replica
  lb<-c(3.75, 2.2)
  ub<-c(3.75, 2.2)
  name  <-"Sacramento"
  tRange<-seq(5,rangeToT)
}


# two rounds of grid search
for (j in 1:2){
  if (j==1){
    #coarser grid
    rList<-seq(lb[1],ub[1],step1[1])
    iList<-seq(lb[2],ub[2],step1[2])
  }else{
    #finer grid
    rList<-seq(g0[1]-step1[1],g0[1]+step1[1],step2[1])
    iList<-seq(g0[2]-step1[2],g0[2]+step1[2],step2[2])   
  }
  gList<-expand.grid(rList,iList)
  ng<-dim(gList)[1]
  
  fitDeath<-matrix(0,ng,nt)
  fitCase <-matrix(0,ng,nt)
  
  start_time <- Sys.time()
  for (i in 1:ng){
    ### regular first
    
    ### parameters
    beta<-gList[i,1]/(infectDuration*eig) 
    vpar["beta"]<-beta
    
    ### initial condition
    initNumIperType<<-gList[i,2]
    inits<-initialCondition2(nc,N,0,NA)
    vt <- seq(0,FROMT,1) 
    
    #run SIR for first FROMT days
    Cmat <<-Cmat1
    sim0 = as.data.frame(lsoda(inits, vt, SEIIRRD_model, vpar))
  
    # shutdown
    inits<-initialCondition2(nc,N,FROMT,sim0)
    vt <- seq(0,TT-1,1)
  
    #run SIR during shutdown
    Cmat <<-Cmat2
    sim1 = as.data.frame(lsoda(inits, vt, SEIIRRD_model, vpar))
  
    #combine SIR results
    sim2<-rbind(sim0[1:FROMT,],sim1)
  
    #store death (per 100k)
    fitDeath[i,]<-extractState("D" ,sim2)*1e3
    fitCase[i,] <-extractSeveralState(c("Ihc","Rq","Rqd"),sim2)*1e3
    print(paste(i,"/",ng,
                ": R0=",format(gList[i,1],digits=4),
                " beta=",format(beta,digits=4),
                " cum death=",format(fitDeath[i,nt],digits=0), sep=""))
  }
  
  end_time <- Sys.time()
  print(end_time - start_time)
  

  
  if (ng>1){
    ## fit based on death
    err<-log(fitDeath) - matrix(1,ng,1)%*%log(covid$deathper100k)
    sse<-rowSums(err[,tRange]^2)
    
    #best fit
    gstar<-which.min(sse)
  }else{
    gstar<-1
  }
  
  g0<-as.double(gList[gstar,1:2])
  beta<-g0[1]/(infectDuration*eig) 
  names(beta)="beta"
  
  
  ### plot comparison
  par(mfrow=c(1,1))
  dead<-log(covid$deathper100k)
  case<-log(covid$caseper100k)
  plot(covid$t, dead, type="l", lwd="2", col="red")
  lines(covid$t, log(fitDeath[gstar,]),type="l")
  abline(v=15, col="gray")
  
  ## check within grid search boundary
  stopifnot(  ((g0[1]<ub[1])*(g0[2]<ub[2])*(g0[1]>lb[1])*(g0[2]>lb[2]))==1)
}


print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("grid search both rounds done!!")
print(paste(" R0=",format(g0[1],digits=4),
            " I0=",format(g0[2],digits=4),
            " beta=",format(beta,digits=4),
            sep=""))
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


# ### produce plots when we settled on parameters
# if (max(ub-lb)==0){
#   # fn <- paste('calibrate_R0_I0_', place, ver, ".pdf", sep="")
#   # pdf(paste(outPath, "figure", fn, sep="/"))
#   # plot(covid$t[tRange], dead[tRange],
#   #      xlab="", ylab="log death per 100k", type="l", lwd="2", col="red", cex.lab=psize, cex.axis=psize)
#   # lines(covid$t[tRange], log(fitDeath[gstar,tRange]),type="l")
#   # abline(v=15, col="gray")
#   # text(max(tRange)*0.7, min(dead[tRange])+0.2*(max(dead[tRange])-min(dead[tRange])),
#   #      paste(name,
#   #            ": beta(%)=",format(round(beta*10000)/100,digits=4),
#   #            ", I0=",format(g0[2],digits=4),sep=""), cex=psize)
#   # dev.off()
#   
  #fit death and case
  fn <- paste('calibrate_R0_I0_', place, ver, ".pdf", sep="")
  pdf(paste(outPath, "figure", fn, sep="/"))
  plot(covid$t, log(fitCase[gstar,]),
       xlab="", ylab="log death and cases per 100k", type="l",ylim=c(-4,10),cex.lab=psize, cex.axis=psize)
  lines(covid$t, dead, type="l", col="red")
  lines(covid$t, case, type="l", col="blue")
  lines(covid$t, log(fitDeath[gstar,]),type="l")
  abline(v=FROMT, col="gray")
  abline(v=max(tRange), col="orange")
  abline(v=min(tRange), col="orange")
  text(max(covid$t)*0.7, min(dead[tRange])+0.2*(max(dead[tRange])-min(dead[tRange])),
       paste(name,
             ": beta(%)=",format(round(beta*10000)/100,digits=3),
             ", I0=",format(exp(g0[2]),digits=3),
             sep=""), cex=)
  dev.off()
# }