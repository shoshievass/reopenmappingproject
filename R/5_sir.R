##############################################################################
##############################################################################
######                                                                  ######
######            reopen mapping project -- run SEIR                    ######
######                                                                  ######
##############################################################################
##############################################################################



### wraper function for running SIR under one policy combination 
run_sir <- function(place, policy, par, sim_ref){
  print(paste("!!.... policy combo", i, ":", 
              policy[[1]], policy[[2]], policy[[3]], "at place", place, "......!!"))
  
  sim0<-NA
  # for each of three phases
  for (j in 1:3){
    contact <- as.vector(policy[[j]])
    # load contact matrix from synthetic population
    Cmat<<-loadData(place, contact)
    
    # load default parameters
    vpar<-vparameters0
    
    #transmission rate
    if (fixBETA==1){
      #fix at initial (high) level
      vpar["beta"]<-par$beta
    }else if (fixBETA==2){
      #fix at reduced level
      vpar["beta"]<-par$beta2        
    }else{
      # consider reduced transmission rate in phase 2 and 3
      vpar["beta"]<-ifelse(j>1,par$beta2,par$beta)
    }
    
    
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
    if (j==3){
      pcombo<-paste(policy[[1]], policy[[2]], policy[[3]], sep="")
      # generate plot/csv outputs
      if (outputSIR==1){
        exportSIR(sim0,place,contact,pcombo)
        packagePlot(sim0,place,pcombo,sim_ref)
      }else{
        # plot SIR output, not save
        par(mfrow=c(2,2))
        
        plotSIR(sim0,NA)
        plotSIRHealth(sim0,NA)
        plotIbyNaics(sim0)
      }
    }
  }
  
  print(paste("cumulative deaths:", round(max(extractState("D",sim0)*sum(types$n)/100))))
  
  return(sim0)
}


#####################################
# run SIR
#####################################
start_time <- Sys.time()

# foreach place
for (m in msaList){
  # load calibrated parameters for this location
  par<-calibratedPar(m)
  
  # initial condition
  initNumIperType<<-par$I0
  
  # for each policy combo
  for (i in 1:dim(policyCombo)[1]){
    if (i==1){
      sim_ref<-run_sir(m, as.vector(policyCombo[i,]), par, NA)
    }else{
      run_sir(m, as.vector(policyCombo[i,]), par, sim_ref)
    }
  }
}
rm(Cmat, par)
end_time <- Sys.time()
print(end_time - start_time)
