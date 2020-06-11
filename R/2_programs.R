##############################################################################
##############################################################################
######                                                                  ######
######        reopen mapping project -- functions                      ######
######                                                                  ######
##############################################################################
##############################################################################



#####################################
#### main SIR functions
#####################################


# SEIIRRD: with exposed, two state for I and R -----------------------------------------------------
SEIIRRD_model=function(t, x, vparameters){
  ncompartment = 9
  nbin = length(x)/ncompartment
  #susceptible, exposed, infected asymptomatic, infected never syptomatic, infected syptomatic, 
  #recovered quarantined, recovered finished quarantined, recovered never quarantined, dead
  S    = as.matrix(x[1:nbin])
  E    = as.matrix(x[(1*nbin+1):(2*nbin)])
  Ia   = as.matrix(x[(2*nbin+1):(3*nbin)])
  Ins  = as.matrix(x[(3*nbin+1):(4*nbin)])
  Ihc  = as.matrix(x[(4*nbin+1):(5*nbin)])
  Rq   = as.matrix(x[(5*nbin+1):(6*nbin)])
  Rqd  = as.matrix(x[(6*nbin+1):(7*nbin)])
  Rnq  = as.matrix(x[(7*nbin+1):(8*nbin)])
  D    = as.matrix(x[(8*nbin+1):(9*nbin)])
  
    S[S  <0] = 0
    E[E  <0] = 0
   Ia[Ia <0] = 0
  Ins[Ins<0] = 0
  Ihc[Ihc<0] = 0
   Rq[Rq <0] = 0
  
  with(as.list(vparameters),{
    N = S+E+Ia+Ins+Ihc+Rq+Rqd+Rnq+D
    
    #mortality rate for symptomatic/asymptomatic
    mortV    = DELTAhc[ageSickVec]
    recvV    = GAMMA[ageSickVec]
    tauV     = TAU[ageSickVec]
    epsilonV = EPSILON[ageSickVec]
    
    #infected per active employed in healthcare
    healthNeed <- as.vector( sum(Ihc) / sum((S+E+Ia+Ins+Rqd+Rnq) * healthVec))
      
    #transition
    dS   = -as.matrix(S*beta*scalBETA)*(as.matrix(Cmat)%*%as.matrix((Ia+Ins)/N)) - as.matrix(S*betaH)*healthNeed*healthVec
    dE   = -epsilonV*as.matrix(E) -dS 
    dIa  = +epsilonV*as.matrix(E) -tauV*as.matrix(Ia)                              
    dIns =                        +tauV*as.matrix(Ia)*(1-psi) -gamma*as.matrix(Ins)  
    dIhc =                        +tauV*as.matrix(Ia)*psi     -recvV*as.matrix(Ihc) -mortV*as.matrix(Ihc)
    dRq  = +recvV*as.matrix(Ihc)   -eta*as.matrix(Rq)
    dRqd =                         +eta*as.matrix(Rq)
    dRnq = +gamma*as.matrix(Ins) 
    dD   = +mortV*as.matrix(Ihc)
    
    # output in the same order as the model compartments are at the beginning of the function
    out=c(dS,dE,dIa,dIns,dIhc,dRq,dRqd,dRnq,dD)
    list(out)
  })
}


# initial condition based on other model otuput -----------------------------------------------------
initialCondition <- function(nc,N,fromT,sim0) {
  #set up initial conditions 
  if (fromT==0){
    # set a constant infection fraction in each type
    Ia_0   <- pmin(N,rep(initNumIperType,nc))
    S_0    <- N-Ia_0
    E_0    <- rep(0,nc)
    Ins_0  <- rep(0,nc)
    Ihc_0  <- rep(0,nc)
    Rq_0   <- rep(0,nc)
    Rqd_0  <- rep(0,nc)
    Rnq_0  <- rep(0,nc)
    D_0    <- rep(0,nc)
  }else{
    # use results from previously run scenario
    
    #start from initT
    iniT<-fromT
    
    #column names
    coln<-colnames(sim0)
    
    #from t=iniT means (iniT+1)'th row of sim1 because sim1 starts from 0
    S_0   = as.double(sim0[iniT+1,grep("S[0-9]"  ,coln, perl=T)])
    E_0   = as.double(sim0[iniT+1,grep("E[0-9]"  ,coln, perl=T)])
    Ia_0  = as.double(sim0[iniT+1,grep("Ia[0-9]" ,coln, perl=T)])
    Ins_0 = as.double(sim0[iniT+1,grep("Ins[0-9]",coln, perl=T)])
    Ihc_0 = as.double(sim0[iniT+1,grep("Ihc[0-9]",coln, perl=T)])
    Rq_0  = as.double(sim0[iniT+1,grep("Rq[0-9]" ,coln, perl=T)])
    Rqd_0 = as.double(sim0[iniT+1,grep("Rqd[0-9]",coln, perl=T)])
    Rnq_0 = as.double(sim0[iniT+1,grep("Rnq[0-9]",coln, perl=T)])
    D_0   = as.double(sim0[iniT+1,grep("D[0-9]"  ,coln, perl=T)])
    stopifnot(length(S_0==nc) && length(Ia_0)==nc && length(D_0)==nc)
  }
  inits=c(S=S_0,E=E_0,Ia=Ia_0,Ins=Ins_0,Ihc=Ihc_0,Rq=Rq_0,Rqd=Rqd_0,Rnq=Rnq_0,D=D_0)

  #return initial condition
  return(inits)
}



# calibrated parameters  -----------------------------------------------------
# store previously calibrated parameters
calibratedPar <- function(place) {
  if (place=="5600"){
    I0<-2.95
    beta<-0.01695
    
    #FRED time-varying beta
    I0<-2.435
    beta<-0.01859
    beta2<-0.00965
  }
  if (place=="1600"){
    #fred
    I0<-0.52
    beta<-0.01085

    #replica    
    I0<-0.2
    beta<-0.002817   
    

    #replica, time varying beta
    I0<-0.1212
    beta<-0.003354
    beta2<-0.00155
    
  }
  if (place=="6920"){
    #fred
    I0<-2.55
    beta<-0.003302
    
    #replica
    I0<-2.2
    beta<-0.0007885
    
    #replica, time varying beta
    I0<-1
    beta<-0.001326
    beta2<-0.0004029
  }
  return(list(I0=I0,beta=beta,beta2=beta2))
}




# load contact matrix data -----------------------------------------------------
loadData <- function(place,contact) {
  
  fn <- paste("C2_msa", place, contact, datv,".csv",sep="")
  
  # check input file exist
  if (max(list.files()==fn)==0){
    stop(paste("missing contact matrix: ", fn, sep=""))
  }
  
  # load files
  CData <-read.csv(fn, header=TRUE)
  
  ## define some global variables of data dimensions and types
  
  
  # number of classes, total population, contact matrix
  nc   <<- dim(CData)[1]
  pop  <<- sum(CData$n)
  scl  <<- 100/pop
  N    <<- CData$n
  types<<-cbind(CData[,c("ego","naics","age","shift","sick")])

  #age X sick (age*10 + sick)
  ageSickVec <<-match(CData$age*10 + CData$sick, typeAgeSick)

  #healthcare worker
  healthVec <<-CData$naics==62
  
  #contact matrix
  Cmat <- as.matrix(CData[,grepl("rate", colnames(CData))])
  
  return(Cmat=Cmat)
}





#####################################
#### other helper functions 
#####################################


# indicate whether the naics is open  -----------------------------------------------------
tagOpenNaics <- function(sim1,contact,place) {
  # load industry policy configuration
  naicsPolicy <-read.csv(paste(parmPath,"naics2essentialpolicy.csv",sep="/"),header=TRUE)
  policy <- gsub("_","", contact)
  
  # add open status to sim1
  for (i in naicsPolicy$naics){
    sim1[[paste("naics",i,sep="")]]<-naicsPolicy[[policy]][naicsPolicy$naics==i]
  }
  return(sim1)
}  



# indicate whether the age is allowed to go to workplace  ------------------------------------------------
tagOpenAge <- function(sim1,contact) {
  policy <- gsub("_","", contact)
  
  # add age open status to sim1
  for (i in unique(types$age)){
    # in the isolate 60+ policy, people 60+ cannot work
    sim1[[paste("age",i,sep="")]]<-ifelse(policy=="shutold" && i>=4,0,1)
  }
  return(sim1)
}  


### compute largest eigen value from contact matrixs -----------------------------------------------------
largestEigenvalue <- function(Cmat) {
  eig <- eigen(Cmat)
  eig <- max(Re(eig$values)) 
  return(eig)
}


# make sure x has 2-dimension and then row sum  -----------------------------------------------------
rowSumMat <- function(x) {
  return(rowSums(as.matrix(x)))
}
# choose desired state by letter (S,I,R,..) and return all classes -------------------
pickState <- function(stateLetter,coln) {
  s<-paste(stateLetter, "[0-9]",sep="")
  return(grep(s,coln, perl=T))
}

# extract the time series of one state according to regular expression (aggregate across all classes) ---------------
extractState <- function(stateLetter,simRun) {
  coln <- colnames(simRun)
  return(scl*rowSumMat(simRun[,pickState(stateLetter,coln)]))
}

# extract the time series of several state according to regular expression (aggregate across all classes) ---------------
extractSeveralState <- function(stateLetterList,simRun) {
  coln <- colnames(simRun)
  s<-0
  for (i in 1:length(stateLetterList)){
    s<-s+scl*rowSumMat(simRun[,pickState(stateLetterList[i],coln)])
  }
  return(s)
}


# plot all state in a list with specified color---------------
plotSeveralLines <- function(what2plot,colvec,simRun,ltype) {
  t<-dim(simRun)[1]-1
  for (i in 1:(length(colvec))){
    lines(0:t,extractSeveralState(what2plot[[i]],simRun),
          type="l",lty=ltype,lwd=2,col=mycol[colvec[i]])
  }
}






#####################################
#### plot/output functions 
#####################################

# plot SIR -----------------------------------------------------
plotSIR <- function(sim1,sim2) {
  
  ### aggregate results across classes
  coln <- colnames(sim1)
  TT<-dim(sim1)[1]-1
  
  # state to plot
  state2plot<-c("Susceptible","Exposed","Infected","Recovered")
  colvec<-c(6,7,1,4)
  what2plot<-list(S="S",E="E",I=c("Ia","Ins","Ihc"),R=c("Rq","Rqd","Rnq"))
  
  plot(0:TT,extractState(what2plot[[1]],sim1)
       ,type="l",ylab="Percent of population",xlab="",ylim=c(0,100),lwd=2,col=mycol[colvec[1]])
  plotSeveralLines(what2plot[-1],colvec[-1],sim1,1)
  
  # secondary outputs
  if (min(is.na(sim2))==0){
    plotSeveralLines(what2plot,colvec,sim2,2)
  }
  # indicate start of policy intervention and reopen
  timingVerticalLine("gray")
  
  legend("topright",
         legend=state2plot,
         col=mycol[colvec],lwd=2,bty="n",cex=1.2)
}




# plot SIR health rated variables -----------------------------------------------------
plotSIRHealth <- function(sim1,sim2) {
  
  ### aggregate results across classes
  coln <- colnames(sim1)
  TT<-dim(sim1)[1]-1
  
  # state to plot
  state2plot<-c("Infected Symptomatic", "Infected Asymptomatic","Employment Loss","Quarantine after Treatment")
  what2plot<-list(Ihc="Ihc",I=c("Ia","Ins"),Emp=c("Ihc","Rq","D"), Rq="Rq")
  colvec<-c(1,6,2,7)
  
  # primary outputs
  plot(0:TT,extractSeveralState(what2plot[[1]],sim1)
       ,type="l",ylab="Percent of population",xlab="",ylim=c(0,40),lwd=2,col=mycol[colvec[1]])
  plotSeveralLines(what2plot[-1],colvec[-1],sim1,1)
  
  
  # secondary outputs
  if (min(is.na(sim2))==0){
    plotSeveralLines(what2plot,colvec,sim2,2)
  }
  # indicate start of policy intervention and reopen
  timingVerticalLine("gray")
  
  # death on secondary axis
  par(new = T)
  plot(0:TT,extractState("D",sim1),type="l",lwd=2,col=mycol[5],
       ylim=c(0,0.8), axes=F, xlab=NA, ylab=NA)
  # secondary outputs
  if (min(is.na(sim2))==0){
    lines(0:TT,extractState("D",sim2),type="l",lty=2,lwd=2,col=mycol[mycolLength+5])
  }
  axis(side  = 4, cex.lab=1.2, cex.axis=1.2)
  mtext(side = 4, line = 2, 'Death percentage', cex=1.2)
  
  legend("topleft",
         legend=c(state2plot,"Death (RHS)"),
         col=mycol[c(colvec,5)],lwd=2,bty="n", cex=1.2)
}


# indicate start of policy intervention and reopen  -----------------------------------------------------
timingVerticalLine <- function(coltxt) {
  abline(v=TTT[2], col=coltxt)
  abline(v=TTT[3], col=coltxt)
}


# plot fraction infected by industry  -----------------------------------------------------
plotIbyNaics <- function(sim1) {
  
  #colnames
  coln<-colnames(sim1)
  TT<-dim(sim1)[1]-1
  
  # foreach industry
  for(i in 1:length(naics2plot)){
    ni<-which(types$naics==naics2plot[i])
    #initial population
    popi<-sum(N[ni])
    #fraction infected
    I <- (100/popi) * (
        rowSumMat(sim1[,coln %in% paste("Ia" ,ni,sep="")]) +
        rowSumMat(sim1[,coln %in% paste("Ins",ni,sep="")]) +
        rowSumMat(sim1[,coln %in% paste("Ihc",ni,sep="")])
    )
    if (i==1){
      plot(0:TT,I
           ,type="l",ylab="Percent infected",xlab="",ylim=c(0,80),lwd=2,col=mycol[i])
    }else{
      lines(0:TT,I,type='l',lwd=2,col=mycol[i])
    }
  }
  legend("topright", legend=naicsName,col=mycol[c(1:7)],horiz=F,lwd=2,bty="n", cex=1.2)
}



# export SIR outputs by each type  -----------------------------------------------------
exportSIR <- function(sim1,place,contact,pcombo) {
  
  #colnames 
  coln<-colnames(sim1)
  
  #time periods
  TT<-dim(sim1)[1]
  
  #types and placeholder for SIR output
  out<-types
  
  #types
  ego  <-types$ego
  naics<-types$naics
  
  # foreach outputs
  compartList    <-c(COMPART, "naics")
  nj<-length(compartList)
  for(j in 1:nj){
    
    outj<-matrix(0,nc,TT)   
    
    # foreach type
    for(i in 1:nc){
      #select output in this NAICS
      if (j==nj){
        x<-sim1[[paste("naics",naics[i],sep="")]]
        # some age bins are not working
        if (naics[i] %in% ESSENTIAL==FALSE){
          x<-pmin(x, sim1[[paste("age",types$age[i],sep="")]])
        }
      }else{
        x<-sim1[[paste(compartList[j],ego[i],sep="")]]
      }
      #record outputs
      outj[i,] <- x
    }
    #export csv
    fn <- paste(outPath, "csv", 
                paste("sir_", compartList[j], "_" ,
                      place, pcombo, ver, ".csv", sep=""), sep="/")
    write.table(cbind(out,outj), file=fn,sep=",",col.names=FALSE,row.names=FALSE)
    print(paste("  export output:",fn))
  }
}



# plot SIR across senarios -----------------------------------------------------
packagePlot <- function(sim1,place,contact, sim2) {
  
  ### produce plots
  fn <- paste('SIR_dcm_', place, contact, ver, ".png", sep="")
  pdf(paste(outPath, "figure", fn, sep="/"))
  plotSIR(sim1,sim2)
  dev.off()
  print(paste("  ...saved",fn))
  
  fn <- paste('SIR2_dcm_', place, contact, ver, ".png", sep="")
  pdf(paste(outPath, "figure", fn, sep="/"))
  plotSIRHealth(sim1,sim2)
  dev.off()
  print(paste("  ...saved",fn))
  
  fn <- paste('Infected_naics_dcm_', place, contact, ver, ".png", sep="")
  pdf(paste(outPath, "figure", fn, sep="/"))
  plotIbyNaics(sim1)
  dev.off()
  print(paste("  ...saved",fn))
}



# color schemes -----------------------------------------------------
t_col <- function(color, percent = 50, name = NULL) {
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  ## Save the color
  invisible(t.col)
}