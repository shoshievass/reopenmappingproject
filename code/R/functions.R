#####################################
#### main SIR functions
#####################################


### LEVEL 1: MAIN FUNCTION WRAPPER
### wrapper function for running SIR under one policy combination 
runSir <- function(place_hist, place_size, c_matrix, policy, par){
  #number of phases
  np<-length(policy)
  
  sim0<-NA
  
  # If T=T_beta
  if (np==2){
  betas <- c(par$beta_0, par$beta_1)
  }
  
  # If T\neq T_beta then three phases
  if (np==3){
    betas <- c(par$beta_0, par$beta_mid ,par$beta_1)
  }
  
  # for each phase
  for (j in 1:np){
    contact <- as.vector(policy[[j]])
    
    # load contact matrix from synthetic population
    Cmat<<-loadData(place_hist, place_size, c_matrix, contact)
    
    
    # load default parameters and set transmission rate
    vpar<-vparameters0
    vpar["beta"] <- betas[j]
    print(vpar["beta"])
    
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
    # stopifnot(min(min(sim0[,grep("S[0-9]"  ,coln, perl=T)]))<-1)
    
    # organize outputs across three phasess
    if (j==np){
      pcombo<-paste(policy, collapse="" )
      
      
      # generate csv outputs
      if (outputSIR==1) exportSIR(sim0,place_hist, place_size,c_matrix,contact,pcombo,par)
    }
  }
  return(sim0)
}

### LEVEL 2: FUNCTIONS CALLED BY runSir directly 
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
    # healthNeed <- as.vector( sum(Ihc) / sum((S+E+Ia+Ins+Rqd+Rnq) * healthVec))
    
    # # transmission separately by types
    # betaV = beta0*(oldSickInd==0) + beta1*(oldSickInd==1)
    
    #transition
    # dS   = -as.matrix(S)*betaV*(as.matrix(Cmat)%*%as.matrix((Ia+Ins)/N)) - as.matrix(S*betaH)*healthNeed*healthVec
    # print(dim(Ia))
    # print(dim(Ins))
    # print(dim(Cmat))
    
    dS   = -as.matrix(S)*beta*(as.matrix(Cmat)%*%as.matrix((Ia+Ins)/N)) 
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



# export SIR outputs by each type  -----------------------------------------------------
exportSIR <- function(sim1, place_hist, place_size, c_matrix, contact, pcombo,par) {
  
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
  compartList    <-c(COMPART, "active","hospital","icu")
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
                    paste("sir_", compartList[j], "_" ,"hist",
                          place_hist,"size", place_size,
                          'contact',c_matrix, pcombo, 
                          "_T", paste(TTT[-1], collapse="-"),
                          "beta_0_", toString(par$beta_0),
                          "_beta_1_", toString(par$beta_1), "_", 
                          verTag, ".csv", sep=""))
    checkWrite(fn, cbind(out,outj), "SIR outputs") 
    
  }
}

# indicate fraction of people in this type that are actively working, 
# and also number of people in hospitalization  -----------------------------------------------------
tagActiveEmp <- function(sim1) {
  # add open status to sim1
  for (i in 1:length(types$ego)){
    sim1[[paste("active"  ,i,sep="")]]<-types$active_emp[i]
    sim1[[paste("hospital",i,sep="")]]<-sim1[[paste("Ia",i,sep="")]] * HOSPITAL[ageSickVec[i]]
    sim1[[paste("icu",i,sep="")]]     <-sim1[[paste("Ia",i,sep="")]] *      ICU[ageSickVec[i]]
  }
  return(sim1)
}  

# load contact matrix data -----------------------------------------------------
loadData <- function(place_hist, place_size, c_matrix, contact) {
  
  # load contact data
  CData <-loadContactData(place_hist, place_size, c_matrix, contact)
  
  # type composition and contact matrix
  types <<-CData$types
  Cmat <-CData$Cmat
  
  # some global variables for types
  globalTypeVectors(types)
  
  return(Cmat=Cmat)
}


### LEVEL 3: FUNCTIONS CALLED BY functions called by runSir 

# load contact matrix data -----------------------------------------------------
loadContactData <- function(place_hist, place_size, c_matrix, contact) {
  
  # load files
  fn <- file.path(contactMatrixPath, 
                  paste(ctMatData, place_hist, 'size', place_size,'contacts' ,c_matrix,contact, ".csv", sep=""))
  CData <-checkLoad(fn)
  ## define some global variables of types of agents
  types<- cbind(CData[,c("ego","naics","age","sick","wfh","shift","work_poi","active_emp","n")])
  
  #contact matrix
  Cmat <- as.matrix(CData[,grepl("rate", colnames(CData))])
  
  return(list(types=types, Cmat=Cmat))
}


# type vectors global variables -----------------------------------------------------
globalTypeVectors <- function(types) {
  
  #age X sick (age*10 + sick)
  ageSickVec <<-match(types$age*10 + types$sick, typeAgeSick)
  
  #healthcare worker
  healthVec <<-(types$naics==heathNAICS) * (types$work_poi==0)
  
  #old/sick people
  oldSickInd <<- (types$sick==1 | types$age>=4)
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
