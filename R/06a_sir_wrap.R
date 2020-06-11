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
# set up
#####################################


baserun<-"_socialdistance"
# baserun<-"_XXX"
vpar<-vparameters0

#####################################
# run SIR
#####################################
start_time <- Sys.time()

# foreach place
for (p in placeList){
  place  <-p
  # calibrated parameters for this location
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
      ### contact matrix from synthetic population
      Cmat<<-loadData(place, contact)
       
      vpar<-vparameters0
      if (j>1){
        vpar["beta"]<-par$beta2
      }else{
        vpar["beta"]<-par$beta
      }
      
      #fix beta to be the same
      if (fixBETA==1){
        vpar["beta"]<-par$beta
      }else if (fixBETA==2){
        vpar["beta"]<-par$beta2        
      }
      
      ### initial condition
      inits<-initialCondition2(nc,N,TTT[j],sim0)
      vt <- seq(0,diff(TTT)[j],1) 
      
      #run SIR
      sim_j = as.data.frame(lsoda(inits, vt, SEIIRRD_model, vpar))
      
      # add naics specific tag for opening sector
      sim_j <- tagOpenNaics(sim_j,contact,place)
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
        # save model results and produce plots
        fn<-paste(policyCombo[i,1],policyCombo[i,2],policyCombo[i,3],sep="")
        if (saveSIR==1){
          save(sim0, file = paste(outPath, "model", paste("ode_", place, fn, ver, ".rda", sep=""), sep="/"))
          print("  saved results")
        }
        # generate plot/csv outputs
        if (outpSIR==1){
          exportSIR(sim0,place,contact,fn)
          # packagePlot2(sim0,place,fn,NA)
          
          if (i==1){
            #first is reference
            sim2<-sim0
            packagePlot2(sim0,place,fn,NA)
          }else{
            packagePlot2(sim0,place,fn,sim2)
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



# copy for matteo:
# cd '/Users/hanyang/Dropbox/COVID Project/R/output/csv'
# cp sir_*_5600_*_combo.csv ../../../Employmentv2/rawdata/SIR/
