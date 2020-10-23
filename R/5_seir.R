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
    parm<-setBeta(contact, par, j) 
    vpar["beta0"] <- parm$beta0
    vpar["beta1"] <- parm$beta1
    
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
    
    # organize outputs across three phasess
    if (j==np){
      pcombo<-paste(policy, collapse="" )
        
      # generate csv outputs
      if (outputSIR==1) exportSIR(sim0,place,contact,pcombo)
      packagePlot(sim0,place,pcombo,simRef)
    }
  }
  return(sim0)
}


#####################################
# run SIR
#####################################
start_time <- Sys.time()


### get max number of policy scenarios across MSAs
P <- checkLoad(policyParm) 
np <- length(grep("Scenario",colnames(P),perl=T))
policyVecName0<-rep("",np)

# foreach MSA
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
  
  #policy combo
  policyCombo<<-genPolicy(reopenPolicy, msaPD$refPolicy)
  
  # demo graphic accounting inputs
  Demo <- checkLoad(demoAcct) 
  Demo <- Demo[Demo$msa==m,]
  
  # aggregate key SEIR outputs
  outcomesofinterest<-c("MSA", 
                      "BaselineAdjDeaths", "BaselineAdjCases", "BaselineAdjEmpDays", "BaselineAdjHospDays",  "BaselineAdjICUDays", 
                      "Deaths", "Case", "EmpDaysLost", "HospitalizationDays", "ICUDays", "Population")
  df<-data.frame(matrix(ncol = np+length(outcomesofinterest)+length(ageNames)*2+length(raceincomeNames)*6, nrow = dim(policyCombo)[1], 
                        dimnames=list(NULL, 
                        c(paste("Policy",1:np, sep=""), outcomesofinterest,
                          paste("Death", ageNames, sep=""),
                          paste("Case" , ageNames, sep=""),
                          paste("DeathRate",      raceincomeNames, sep=""),
                          paste("HospitalDays",   raceincomeNames, sep=""),
                          paste("ICUDays",        raceincomeNames, sep=""),
                          paste("AvgEmpDaysLost", raceincomeNames, sep=""),
                          paste("AvgAge",         raceincomeNames, sep=""),
                          paste("Prob",           raceincomeNames, sep="")))), stringsAsFactors=FALSE)

  # for each policy combo
  simRef<-NA
  for (i in 1:dim(policyCombo)[1]){
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(paste("!! Running SIR for policy combo", i,"/", dim(policyCombo)[1], 
                ":", paste(policyCombo[i,], collapse=""),"at MSA", m, "......!!"))
    
    ## sir across phases
    sim_i<-runSir(m, as.vector(policyCombo[i,]), par, simRef)
    outstats_i<-calOutcome(sim_i, Demo)
    
    #the first is reference policies in each MSA
    if (i==1){
      simRef<-sim_i
      outstats_ref<-outstats_i
    }
    
    ## collect outputs
    policyVecName      <- policyVecName0
    policyVecName[1:np]<- gsub("_","",policyCombo[i,])
    df[i,]<-c(policyVecName, m, outstats_i[1:5]/outstats_ref[1:5]-1,outstats_i)
  }
  
  #export csv output for this msa
  checkWrite(file.path(outPath, paste("Agg_SIR_",m,".csv", sep="")), 
             df, "aggregate SIR outputs") 
}
rm(Cmat, par)




end_time <- Sys.time()
print(end_time - start_time)
