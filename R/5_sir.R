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
    vpar["beta"]<-setBeta(contact, par, j) 
    
    ### initial condition
    inits<-initialCondition(TTT[j],sim0)
    vt <- seq(0,diff(TTT)[j],1) 
    
    #run SIR
    sim_j = as.data.frame(lsoda(inits, vt, SEIIRRD_model, vpar))
    
    # add tag to keep track of which type is working
    sim_j <- tagActiveEmp(sim_j)

    #aggregate results
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

# aggregate key outputs
df<-data.frame(matrix(ncol = np+6, nrow = 0, 
                  dimnames=list(NULL, 
                        c(paste("Policy",1:np, sep=""), 
                          "MSA", "BaselineAdjDeaths", "BaselineAdjEmpHours", "Deaths", "EmpHoursLost", "Population"))
                  ), stringsAsFactors=FALSE)
dfRow<-0
policyVecName0<-rep("",np)

# foreach place
for (m in msaList){
  # load calibrated parameters for this location
  par<-calibratedPar(m)

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
  
  
  
  # for each policy combo
  simRef<-NA
  for (i in 1:dim(policyCombo)[1]){
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(paste("!! Running SIR for policy combo", i,"/", dim(policyCombo)[1], 
                ":", paste(policyCombo[i,], collapse=""),"at MSA", m, "......!!"))
    
    ## sir across phases
    sim_i<-runSir(m, as.vector(policyCombo[i,]), par, simRef)
    outstats_i<-calOutcome(sim_i)
    
    #the first is reference policies in each MSA
    if (i==1){
      simRef<-sim_i
      outstats_ref<-outstats_i
    }
    
    ## collect outputs
    dfRow<-dfRow+1
    policyVecName <- policyVecName0
    policyVecName[1:np]<-gsub("_","",policyCombo[i,])
    df[dfRow,]<-c(policyVecName, m, outstats_i[1:2]/outstats_ref[1:2]-1,outstats_i)
  }
}
rm(Cmat, par)


#export csv
checkWrite(file.path(outPath, paste("Agg_SIR", datv, ".csv", sep="")), 
           df, "aggregate SIR outputs") 

end_time <- Sys.time()
print(end_time - start_time)
