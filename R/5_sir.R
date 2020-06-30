##############################################################################
##############################################################################
######                                                                  ######
######            reopen mapping project -- run SEIR                    ######
######                                                                  ######
##############################################################################
##############################################################################



### wraper function for running SIR under one policy combination 
run_sir <- function(place, policy, par, sim_ref){
  sim0<-NA
  # for each of three phases
  for (j in 1:4){
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
    if (j==4){
      pcombo<-paste(policy[[1]], policy[[2]], policy[[3]], policy[[4]], sep="")
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
  
  return(sim0)
}


#####################################
# run SIR
#####################################
start_time <- Sys.time()

# aggregate key outputs
df <- data.frame(Policy1=character(),Policy2=character(),Policy3=character(),Policy4=character(),
                 MSA=character(), 
                 BaselineAdjDeaths=numeric(), 
                 BaselineAdjEmpHours=numeric(), 
                 Deaths=numeric(), 
                 EmpHoursLost=numeric(), 
                 Population=numeric(), 
                 stringsAsFactors=FALSE) 
dfRow<-0

# foreach place
for (m in msaList){
  # load calibrated parameters for this location
  par<-calibratedPar(m)
  
  # load location input parameters
  gsPar <- checkLoad(gsParm)
  gsPar <- as.vector(gsPar[gsPar$msa==m,])
  
  # timing for 4 phases
  TTT<<-c(0,
          gsPar$T2-gsPar$T1,
          gsPar$T3-gsPar$T1,
          gsPar$T4-gsPar$T1,
          gsPar$Tend)
  
  # initial condition
  initNumIperType<<-par$I0
  
  # for each policy combo
  for (i in 1:dim(policyCombo)[1]){
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(paste("!! Running SIR for policy combo", i,"/", dim(policyCombo)[1], ":", 
                policyCombo[i,1], 
                policyCombo[i,2], 
                policyCombo[i,3], 
                policyCombo[i,4], "at MSA", m, "......!!"))
    
    if (i==1){
      #the first is reference policies in each MSA
      sim_ref<-run_sir(m, as.vector(policyCombo[i,]), par, NA)
      outstats_ref<-calOutcome(sim_ref)
    }else{
      sim_i<-run_sir(m, as.vector(policyCombo[i,]), par, sim_ref)
      outstats_i<-calOutcome(sim_i)
      
      ## collect outputs
      dfRow<-dfRow+1
      df[dfRow,]<-c(gsub("_","",policyCombo[i,]), 
                   m, outstats_i[1:2]/outstats_ref[1:2]-1,outstats_i)
    }
  }
}
rm(Cmat, par)


#export csv
fn <- paste(outPath, 
            paste("Agg_SIR", datv, ".csv", sep=""), sep="/")
write.table(df, file=fn, sep=",",col.names=TRUE,row.names=FALSE)
print(paste("aggregate SIR outputs:",fn))

end_time <- Sys.time()
print(end_time - start_time)
