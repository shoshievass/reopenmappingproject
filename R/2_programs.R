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
    
    # output in the same order as the model compartments 
    out=c(dS,dE,dIa,dIns,dIhc,dRq,dRqd,dRnq,dD)
    list(out)
  })
}


# initial condition based on other model otuput -----------------------------------------------------
initialCondition <- function(fromT,sim0) {
  #number of types and population
  N <- types$n
  nc<- length(N)
  
  #set up initial conditions 
  if (fromT==0){
    # set a constant infection fraction in each type
    Ia_0   <- N*pmin(initNumIperType/1e3,1)
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


#####################################
#### other helper functions 
#####################################

# calibrated parameters  -----------------------------------------------------
# load previously calibrated parameters
calibratedPar <- function(place) {
  #load calibrated parameter file
  fn <- file.path(calibratedParPath, 
                  paste(caliParm, place, ".csv", sep=""))
  pData <-checkLoad(fn)
  
  ## extract parameters
  I0<-as.vector(min(pData$I0))
  beta1<-as.vector(min(pData$beta1))
  beta2<-as.vector(min(pData$beta2))

  return(list(I0=I0,beta1=beta1,beta2=beta2))
}


# load contact matrix data -----------------------------------------------------
loadData <- function(place,contact) {
  
  #remove mask assumption from contact code
  contact <- substr(contact, 1, gregexpr("-M",contact)[[1]][1]-1)
  
  # load files
  fn <- file.path(contactMatrixPath, 
              paste(ctMatData, place, contact, ".csv", sep=""))
  CData <-checkLoad(fn)

  ## define some global variables of types of agents
  # types
  types<<- cbind(CData[,c("ego","naics","age","sick","wfh","shift","active_emp","n")])
  
  #age X sick (age*10 + sick)
  ageSickVec <<-match(CData$age*10 + CData$sick, typeAgeSick)

  #healthcare worker
  healthVec <<-CData$naics==heathNAICS
  
  #contact matrix
  Cmat <- as.matrix(CData[,grepl("rate", colnames(CData))])
  
  return(Cmat=Cmat)
}


# check if file exist and load  -----------------------------------------------------
checkLoad <- function(fn) {
  if (!file.exists(fn)){
    stop(paste("missing input file : ", fn, sep=""))
  }
  return(read.csv(fn, header=TRUE))
}


# write csv, check if file already exist  -----------------------------------------------------
checkWrite <- function(fn, out, msg){
  if (file.exists(fn)) print(paste("overwrite existing file"))
  write.table(out, file=fn, sep=",",col.names=TRUE,row.names=FALSE)
  print(paste("save",msg,":",fn))
}


# get code directory  -----------------------------------------------------
getCodePath <- function(filename){
  out <- file.path(codePath,filename)
  return(out)
}


# check contact matrix source and load -----------------------------------------------------
loadCmat <- function(cFn){
  
  fnReplica<-file.path(dataPath, paste(cFn, "_replica.csv",sep=""))
  fnFred   <-file.path(dataPath, paste(cFn, "_fred.csv"   ,sep=""))
  
  #first try replica then fred
  if(file.exists(fnReplica)){
    C<-read.csv(fnReplica, header=TRUE) 
    print("load replica contact matrix")
  }else if(file.exists(fnFred)){
    C<-read.csv(fnFred, header=TRUE)
    print("load fred contact matrix")
  }else{
    stop(paste("missing", cFn))
  } 
  return(C)
}




# load policy scenario and dates -----------------------------------------------------
loadPolicyDates <- function(m){
  
  # load policy and dates for this msa
  P <- checkLoad(policyParm) 
  P <- P[P$MSA==m,]
 
  # variable names
  colNm <- colnames(P)
  colNm <- colNm[!(colNm %in% "Tend")]
  
  # policy scenarios and dates
  policies <- grep("Scenario",colNm,perl=T)
  dates    <- grep("T",colNm,perl=T)
  if(length(policies)!=length(dates)) stop("number of policies and number dates not match!")
  
  # for each phase
  np<-length(policies)
  policyVec<-c()
  TVec<-c()
  for (i in 1:np){
    poli_i <- as.character(P[[colNm[policies[i]]]])
    date_i <- P[[colNm[dates[i]]]]
    
    #policy name tag
    policyVec<-c(policyVec,paste("_",poli_i,policyCaliEM[min(i,length(policyCaliEM))],sep=""))
    
    #start of this policy
    TVec<-c(TVec,date_i)
  }
  
  # end dates
  TVec<-c(TVec, P$Tend)
  if(min(diff(TVec))<=0) stop("policy start dates not increasing")
  
  return(list(refPolicy=policyVec, TVec=TVec))
}



# indicate fraction of people in this type that are actively working  -----------------------------------------------------
tagActiveEmp <- function(sim1) {
  # add open status to sim1
  for (i in 1:length(types$ego)){
    sim1[[paste("active",i,sep="")]]<-types$active_emp[i]
  }
  return(sim1)
}  



# generate string policy tag names from indicator for definition of each contact level -------------------------------------
policyTagString <- function(policyVec) {
  # number of policies
  n <- dim(policyVec)[2]
  policyVec <- matrix(as.matrix(policyVec), ncol=n)
  
  # policy tag
  poli <- ""
  for (i in 1:n){
    tag <- ifelse(i==1, "_", "-")
    poli<- paste(poli, tag, policyLetterCode[i], policyVec[,i], sep="")
  }
  return(poli)
}


# extract policy indciator  -------------------------------------
parsePolicyTag <- function(policyTag, poli) {
  #position of letter code "poli"
  pos <- gregexpr(poli,policyTag)[[1]][1] + 1
  return(as.numeric(substr(policyTag, pos, pos)))
}


# policy combo matrix  -------------------------------------
policyMatrix <- function(refPolicy,np){
  # we gradually incorporate reference policies, so we have the benchmark at each phase
  X <-matrix(0,np,np)
  for (i in 1:np){
    for (j in 1:np){
      X[np+1-i,j]<-min(i,j)
    }
  }
  
  # number of policy to string
  return(matrix(refPolicy[X],np,np))
}



# enumerate all policies to plot, including references  -------------------------------------
genPolicy <- function(reopenPolicy,refPolicy) {
  
  ## number of phases
  np <-length(refPolicy)
  
  ## reference policies
  refp <- policyMatrix(refPolicy,np)
  
  ## different reopen policies
  policyList <<- policyTagString(reopenPolicy)
  nr <- length(policyList)
  
  ## combinations of reference and reopen policy in the last phase
  policyFull <- t(unname(rbind(matrix(rep(refPolicy[1:(np-1)],nr),(np-1),nr),policyList)))

  ## all policy combo to run
 if (genRef4Area==1){
    # generate more reference policies for area plots
    policyCombo<<-rbind(refp,policyFull)
  }else{
    policyCombo<<-rbind(refp[1,],policyFull)
    
  }
 
  return(policyCombo)
}


# set transmission rate beta in SIR model  -------------------------------------
setBeta <- function(contact, par, j) {
  #transmission rate
 
  
  if (fixBETA==1){
    #fix at initial (high) level
    beta<-par$beta1
  }else if (fixBETA==2){
    #fix at reduced level
    beta<-par$beta2        
  }else{
    #any policy assumption on transmission rate
    beta_i <- parsePolicyTag(contact, "M")
    
    if (is.na(beta_i)){
      # consider reduced transmission rate in phase 2 and 3
      beta<-ifelse(j>1,par$beta2,par$beta1)
    }else{
      if (beta_i==1){
        #reduced beta (status quo in phase 2 and 3)
        beta<-par$beta2
      }else if (beta_i==2){
        #more beta than status quo
        beta<-par$beta2*(2/3) + par$beta1*(1/3)
      }else if (beta_i==3){
        #less beta than initial
        beta<-par$beta2*(1/3) + par$beta1*(2/3)
      }else if (beta_i==4){
        #initial beta (status quo in phase 1)
        beta<-par$beta1
      }else(
        stop("only support M in 1-4")
      )
    }
  }
  return(beta)
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
  return((100/sum(types$n))*rowSumMat(simRun[,pickState(stateLetter,coln)]))
}

# extract the time series of several state according to regular expression (aggregate across all classes) ---------------
extractSeveralState <- function(stateLetterList,simRun) {
  coln <- colnames(simRun)
  s<-0
  for (i in 1:length(stateLetterList)){
    s<-s+(100/sum(types$n))*rowSumMat(simRun[,pickState(stateLetterList[i],coln)])
  }
  return(s)
}


# extract final outcome of one state by age ---------------
extractFinalByAge <- function(stateLetter,simRun) {
  coln <- colnames(simRun)
  
  # age groups
  # we skip age group 0: 0-4, since we don't have that in replica
  ageVec<-unique(types$age[types$age>0])
  ageOut<-rep(0,length(ageVec))
  tt <- dim(simRun)[1]
  
  #for each age cohort
  for (i in 1:length(ageVec)){
    ni<-which(types$age==ageVec[i])
    #sum axross all types of this age
    ageOut[i] <-  rowSumMat(simRun[tt,coln %in% paste(stateLetter, ni,sep="")]) 
  }
  return(ageOut)
}



# plot all state in a list with specified color---------------
plotSeveralLines <- function(what2plot,colvec,simRun,ltype) {
  t<-dim(simRun)[1]-1
  for (i in 1:(length(colvec))){
    lines(0:t,extractSeveralState(what2plot[[i]],simRun),
          type="l",lty=ltype,lwd=2,col=mycol[colvec[i]])
  }
}


# compute and print key health and employment outcome of interest---------
calOutcome<-function(simRun){
 
  #1. cumualtive deaths and cases
  scal <- sum(types$n)/100
  deaths <-max(extractState("D",simRun)*scal)
  print(paste("Deaths:", round(deaths)))
  
  cases <-max(extractSeveralState(c("Ihc","Rq","Rqd","D"),simRun)*scal)
  print(paste("Cases:", round(cases)))
  
  #death/cases by age
  print("Deaths in each age group")
  deathByAge <- extractFinalByAge("D",simRun)
  print(deathByAge)
  
  print("Cases in each age group")
  caseByAge <- extractFinalByAge("Ihc", simRun) + 
               extractFinalByAge("Rq" , simRun) +
               extractFinalByAge("Rqd", simRun) + deathByAge
  print(caseByAge)
  
  #2. count employment less
  coln <- colnames(simRun)
  
  ## includes compartment that are working
  s<-simRun[,pickState("S",coln)]
  for (i in c("E","Ia","Ins","Rqd","Rnq")){
    s<-s+simRun[,pickState(i,coln)]
  }
  ## consider wether these individuals are actively working
  active<-s * simRun[,pickState("active",coln)]
  
  ## compute total employment loss in days
  t<-max(TTT)
  empLoss <- sum(types$n[types$naics>0])*t - sum(active[2:(t+1),])
  print(paste("Employment loss (1000days):", 
              round(empLoss/1e3)))  
  
  ## cases
  case<-max(extractSeveralState(c("Ihc","Rq","Rqd","D"),simRun))*1e3
  print(paste("Cumulative cases (per100k):", 
              round(case))) 
  
  ## stack outputs
  out <- matrix(c(deaths, cases, empLoss, sum(types$n), deathByAge, caseByAge),
                1,4+length(caseByAge)*2)
  
  return(out)
}



#####################################
#### plot/output functions 
#####################################




# plot SIR -----------------------------------------------------
plotSIR <- function(fn, sim1,sim2) {
  
  if (fn!="")  png(fn)

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
         col=mycol[colvec],horiz=F,lwd=1.3,bty="n",cex=0.8)
  
  if (fn!=""){
    dev.off()
    print(paste("  saved plot:",fn))
  }
}


# plot SIR health rated variables -----------------------------------------------------
plotSIRHealth <- function(fn, sim1,sim2) {
  
  if (fn!="")  png(fn)
  
  ### aggregate results across classes
  coln <- colnames(sim1)
  TT<-dim(sim1)[1]-1
  
  # state to plot
  state2plot<-c("Employment Loss","Infected Symptomatic", "Infected Asymptomatic","Quarantine after Treatment")
  what2plot<-list(Emp=c("Ihc","Rq","D"), Ihc="Ihc",I=c("Ia","Ins"),Rq="Rq")
  colvec<-c(2,1,6,7)
  
  # primary outputs
  plot(0:TT,extractSeveralState(what2plot[[1]],sim1)
       ,type="l",ylab="Percent of population",xlab="",lwd=2,col=mycol[colvec[1]])
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
       axes=F, xlab=NA, ylab=NA)
  # secondary outputs
  if (min(is.na(sim2))==0){
    lines(0:TT,extractState("D",sim2),type="l",lty=2,lwd=2,col=mycol[mycolLength+5])
  }
  axis(side  = 4, cex.lab=1.2, cex.axis=1.2)
  mtext(side = 4, line = 2, 'Death percentage', cex=1.2)
  
  legend("topleft",
         legend=c(state2plot,"Death (RHS)"),
         col=mycol[c(colvec,5)],horiz=F,lwd=1.3,bty="n",cex=0.8)
  
  if (fn!=""){
    dev.off()
    print(paste("  saved plot:",fn))
  }
}


# indicate start of policy intervention and reopen  -----------------------------------------------------
timingVerticalLine <- function(coltxt) {
  for (i in 2:(length(TTT)-1)){
    abline(v=TTT[i], col=coltxt)
  }
}


# plot fraction infected by industry  -----------------------------------------------------
plotIbyNaics <- function(fn, sim1) {
  
  if (fn!="")  png(fn)
  
  #colnames
  coln<-colnames(sim1)
  TT<-dim(sim1)[1]-1
  
  # foreach industry
  for(i in 1:length(naics2plot)){
    ni<-which(types$naics==naics2plot[i])
    #initial population
    popi<-sum(types$n[ni])
    #fraction infected
    I <- (100/popi) * (
        rowSumMat(sim1[,coln %in% paste("Ia" ,ni,sep="")]) +
        rowSumMat(sim1[,coln %in% paste("Ins",ni,sep="")]) +
        rowSumMat(sim1[,coln %in% paste("Ihc",ni,sep="")])
    )
    if (i==1){
      plot(0:TT,I
           ,type="l",ylab="Percent infected",xlab="",lwd=2,col=mycol[i])
    }else{
      lines(0:TT,I,type='l',lwd=2,col=mycol[i])
    }
  }
  legend("topright", legend=naicsName,col=mycol[c(1:7)],horiz=F,lwd=1.3,bty="n",cex=0.8)
  
  if (fn!=""){
    dev.off()
    print(paste("  saved plot:",fn))
  }
}



# export SIR outputs by each type  -----------------------------------------------------
exportSIR <- function(sim1,place,contact,pcombo) {
  
  #colnames 
  coln<-colnames(sim1)
  
  #time periods
  TT<-dim(sim1)[1]
  
  #types and placeholder for SIR output
  out<-types[,!(colnames(types) %in% c("active_emp"))]
  
  #types
  ego<-types$ego
  
  #number of types
  nc<-length(ego)
  
  # foreach outputs
  compartList    <-c(COMPART, "active")
  nj<-length(compartList)
  for(j in 1:nj){
    
    outj<-matrix(0,nc,TT)   
    
    # foreach type, record SIR output at each T
    for(i in 1:nc){
      x<-sim1[[paste(compartList[j],ego[i],sep="")]]
      #record outputs
      outj[i,] <- x
    }
    colnames(outj)<-paste("T",c(0:(TT-1)),sep="")
    
    # export csv
    fn <- file.path(outPath, "csv", 
                paste("sir_", compartList[j], "_" ,
                      place, pcombo, 
                      "_T", paste(TTT[-1], collapse="-"), 
                      verTag, ".csv", sep=""))
    checkWrite(fn, cbind(out,outj), "SIR outputs") 
    
  }
}



# plot SIR across senarios -----------------------------------------------------
packagePlot <- function(sim1,place,contact, sim2) {
  
  ### produce plots, just print or save
  if (outputSIR==1){
    fnEnd <- paste(place, contact, verTag, ".png", sep="")
    fn1 <- file.path(outPath, "figure", paste('SIR_dcm_', fnEnd, sep=""))
    fn2 <- file.path(outPath, "figure", paste('SIR2_dcm_', fnEnd, sep=""))
    fn3 <- file.path(outPath, "figure", paste('Infected_naics_dcm_', fnEnd, sep=""))
  }else{
    fn1 <-""
    fn2 <-""
    fn3 <-""
  }
  par(mfrow=c(2,2))
  plotSIR(fn1,sim1,sim2)
  plotSIRHealth(fn2,sim1,sim2)
  plotIbyNaics(fn3,sim1)
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

# generate color scheme  -----------------------------------------------------
gen_col <- function() {
  ## set of colors
  mycol <<- c(t_col("red",     perc = 0),
              t_col("orange",  perc = 0),
              t_col("yellow2", perc = 0),
              t_col("green",   perc = 0),
              t_col("dimgray", perc = 0),
              t_col("blue",    perc = 0),
              t_col("purple",  perc = 0))
  mycolLength<<-7
  psize<<-1.25
}