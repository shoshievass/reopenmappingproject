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
# SET UP
#####################################
policyCombo<<-expand.grid(c("_regular"),c("_socialdistance"),
                          c("_regular","_socialdistance","_alternateschoolwork","_shutold", "_wfhreopenschool","_cautious"))
npolicy<-dim(policyCombo)[1]

### MSA list
placeList<-c("5600", "1600", "6920")
placeList<-c("5600")
MSAstats<-matrix(0,length(placeList),8+npolicy)

### death statistics
COVID <-read.csv("covid_case_death_jh.csv",header=TRUE)



for (p in 1:length(placeList)){  
  place<-placeList[p]
  MSAstats[p,1]<-as.numeric(place)
  
  print(paste("running msa ",place,sep=""))
  
  #####################################
  # Grid Search
  #####################################

  
  ### load covid death and cases data
  covid<-COVID[COVID$msa==place,c("t","deathper100k","msaname")]
  nt<-dim(covid)[1]
  #msa name and cum deaths at t=90
  MSAstats[p,2]<-as.character(covid$msaname)[1]
  MSAstats[p,3]<-covid$deathper100k[covid$t==nt]
  
  #run this msa if there is meaningful deaths
  if (covid$deathper100k[60]>0){
    vpar<-vparameters0
    
    ### regular and shutdown contact matrix
    Cmat1<-loadData(place, "_regular")
    eig<-largestEigenvalue(Cmat1)
    Cmat2<-loadData(place, "_socialdistance")
    
    
    ### grid search
    # back out beta from the R0 and tau
    lb<-c(0,  -5)
    ub<-c(60,  4)
    step1<-c(10, 1)
    step<-rbind(step1, step1*0.1, step1*0.01)
    
    #change contact at t=15
    TTT <<-c(0,15,nt-1)
    #range of sample to fit deaths
    t0<-max(30,min(which(covid$deathper100k>0))+2)
    t1<-90
    tRange<-seq(t0,t1)
    
    # j rounds of grid search
    for (j in 1:3){
      if (j==1){
        #coarser grid
        rList<-seq(lb[1],ub[1],step[1,1])
        iList<-seq(lb[2],ub[2],step[1,2])
      }else{
        #finer grid
        rList<-seq(max(lb[1],g0[1]-step[j-1,1]),
                   min(ub[1],g0[1]+step[j-1,1]),
                   step[j-1,1])
        iList<-seq(max(lb[2],g0[2]-step[j-1,2]),
                   min(ub[2],g0[2]+step[j-1,2]),
                   step[j-1,2]) 
      }
      gList<-expand.grid(pmax(0.01,rList),iList)
      
      ng<-dim(gList)[1]
      fitDeath<-matrix(0,ng,nt)
      
      for (i in 1:ng){
        ### parameters
        beta <-gList[i,1]/(infectDuration*eig) 
        vpar["beta"] <-beta
        
        ### initial condition
        I0<-exp(gList[i,2])
        initNumIperType<<-I0
        
        ### change contact matrix after t=15
        sim0<-NA
        for (t in 1:2){
          inits<-initialCondition2(nc,N,TTT[t],sim0)
          vt <- seq(0,diff(TTT)[t],1) 
          
          if (t==1){
            Cmat <<-Cmat1
          }else{
            Cmat <<-Cmat2        
          }
          
          # RUN SIR
          sim_j = as.data.frame(lsoda(inits, vt, SEIIRRD_model, vpar))
          
          if (t>1){
            sim0<-rbind(sim0[1:TTT[t],], sim_j)
          }else{
            sim0<-sim_j
          }      
        }
        
        #store death  (per 100k)
        fitDeath[i,]<-extractState("D",sim0)*1e3
      }
    
      ## fit log death, min squared loss
      err<-log(fitDeath) - matrix(1,ng,1)%*%log(covid$deathper100k)
      sse<-rowSums(err[,tRange]^2)
      
      #best fit
      gstar<-which.min(sse)
      
      g0<-as.double(gList[gstar,])
      beta <-g0[1]/(infectDuration*eig) 
      I0<-exp(g0[2])
      
      if (min((g0<ub) * (g0>lb))==1){
        ubslack<-c(0,0)
        lbslack<-c(0,0)
      }else{
        ubslack<-ifelse(g0>=ub,1,0)
        lbslack<-ifelse(g0<=lb,1,0)
        lb <- lb -lbslack*step[j,]
        ub <- ub +ubslack*step[j,]
        print("hit grid boundary")
      }
    }
    
    plot(covid$t, log(covid$deathper100k), type="l")
    lines(covid$t, log(fitDeath[gstar,]),type="l", col="red")
    print(paste(" beta=", format(beta, digits=4)," I0=",   format(I0,digits=4),sep=""))
    
    #some stats to keep track
    MSAstats[p,4:8]<-c(fitDeath[gstar,nt],t0,t1,beta,I0)
  }else{
    #just use fixed default value if there is not much deaths
    beta<-0.0005
    I0  <-0.1
    MSAstats[p,6:7]<-c(beta,I0)
  }
  
  
  #####################################
  # SIR
  #####################################
  
  ## time period
  TTT <<-c(0,15,75,150)
  
  # initial condition
  initNumIperType<<-I0
  vpar<-vparameters0
  vpar["beta"]<-beta

  # for each policy combo
  for (i in 1:npolicy){
    print(paste("!!.... policy combo", i, ":", policyCombo[i,1],policyCombo[i,2],policyCombo[i,3], "at place", place, "......!!"))

    sim0<-NA
    # for each phase
    for (j in 1:3){
      contact <- as.vector(policyCombo[i,j])

      ### contact matrix from synthetic population
      Cmat<<-loadData(place, contact)

      ### initial condition
      inits<-initialCondition2(nc,N,TTT[j],sim0)
      vt <- seq(0,diff(TTT)[j],1)

      ### run SIR
      sim_j = as.data.frame(lsoda(inits, vt, SEIIRRD_model, vpar))

      ### add naics specific tag for opening sector
      sim_j <- tagOpenNaics(sim_j,contact,NA)
      ### add age specific tag for whether people in this age bin can work
      sim_j <- tagOpenAge(sim_j,contact)


      ### aggregate results
      if (j>1){
        sim0<-rbind(sim0[1:TTT[j],], sim_j)
      }else{
        sim0<-sim_j
      }

      # check population is static
      check<-rowSums(sim0[,c("S5","E5","Ia5","Ins5","Ihc5","Rq5","Rnq5","Rqd5","D5")])
      stopifnot(abs(min(check)-max(check))<0.1)

      # organize outputs across three phasess
      if (j==3){
        # save model results and produce plots
        fn<-paste(policyCombo[i,1],policyCombo[i,2],policyCombo[i,3],sep="")
        # generate plot/csv outputs
        # exportSIR(sim0,place,contact,fn)
      }
    }

    MSAstats[p,8+i]<-max(extractState("D",sim0))
  }
}
