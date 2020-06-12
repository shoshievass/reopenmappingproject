##############################################################################
##############################################################################
######                                                                  ######
######            reopen mapping project -- run SEIR                    ######
######                                                                  ######
##############################################################################
##############################################################################



#####################################
# set up
#####################################


vpar<-vparameters0

#####################################
# run SIR
#####################################
start_time <- Sys.time()

# foreach place
for (m in msaList){
  place  <-m
  # load calibrated parameters for this location
  par<-calibratedPar(place)
  
  # initial condition
  initNumIperType<<-par$I0
  
  # for each policy combo
  for (i in 1:dim(policyCombo)[1]){
    print(paste("!!.... policy combo", i, ":", policyCombo[i,1],policyCombo[i,2],policyCombo[i,3], "at place", place, "......!!"))
    
    sim0<-NA
    # for each phase
    for (j in 1:3){
      contact <- as.vector(policyCombo[i,j])
      # load contact matrix from synthetic population
      Cmat<<-loadData(place, contact)
       
      # load inputed parameters
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
      inits<-initialCondition(nc,N,TTT[j],sim0)
      vt <- seq(0,diff(TTT)[j],1) 
      
      #run SIR
      sim_j = as.data.frame(lsoda(inits, vt, SEIIRRD_model, vpar))
      
      # add naics specific tag to keep track of opening sector
      sim_j <- tagOpenNaics(sim_j,contact)
      # add age specific tag for whether people in this age bin can work
      sim_j <- tagOpenAge(sim_j,contact)
      
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
        pcombo<-paste(policyCombo[i,1],policyCombo[i,2],policyCombo[i,3],sep="")
        print(pcombo)
        # generate plot/csv outputs
        if (outputSIR==1){
          exportSIR(sim0,place,contact,pcombo)
          if (i==1){
            #first run is reference
            sim2<-sim0
            packagePlot(sim0,place,pcombo,NA)
          }else{
            packagePlot(sim0,place,pcombo,sim2)
          }
        }else{
          # plot SIR output, not save
          par(mfrow=c(2,2))
          
          plotSIR(sim0,NA)
          plotSIRHealth(sim0,NA)
          plotIbyNaics(sim0)
        }
      }
    }
  }
}

end_time <- Sys.time()
print(end_time - start_time)
