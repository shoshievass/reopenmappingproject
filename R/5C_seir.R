##############################################################################
##############################################################################
######                                                                  ######
######            reopen mapping project -- run SEIR                    ######
######                   construct frontier                             ######
######                                                                  ######
##############################################################################
##############################################################################


### wrapper function for running SIR under one policy combination 
runSir <- function(place, policy, par){
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
    vpar["beta"] <- setBeta(contact, par) 
    
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
  }
  return(sim0)
}


#####################################
# run SIR
#####################################
start_time <- Sys.time()




###################################
## HY: stuff we want to switch
### different work policies
#do nothing, lock down, wfh/as
work<-c(0,1,2)

## number of periods and number of days in between policies
nb<-8

ndays<-14

## number of week in tottal
Nweek <<-round(nb*ndays/7)

## possible grouping of breaks
policies<-permutations(length(work), nb, v=work, repeats.allowed=TRUE)

## different work policies
policyCombo<-policies[apply(t(diff(t(policies))), 1, FUN=min) ==0,]

## wait and do nothing is never optimal
# initial condition does not seem to matter
policyCombo <- rbind(policyCombo[1,], policyCombo[policyCombo[,1]>0,])

## the rest of policies
# W(work), 
# S(school), N(neighbor), POI-resturants, POI-retail, POI-personal service, POI-entertainment, (E)eldery, mask=1
# HY: always set mask to 1, we can test different policies for school, etc and see what frontier looks like
# note due to large file size of contact matrices, I only commited S=1 or 3, N=2 and E=2
# so the possible non-work policies to explore are c({1,3},2,{1,2},{1,2},{1,2},{1,2},2,1)
# work policies has 1 (essential only), 2 (wfh), 3 (alternating schedule), 4(full reopen)
policyNonWork<-c(1,3,
                 2,2,2,2,
                 2,1)

policyNonWorkList<-rbind(c(1,3,2,2,2,2,2,1),
                         c(1,3,1,2,2,1,2,1),
                         c(1,3,1,1,1,1,2,1),
                         c(1,2,2,2,2,2,2,1),
                         c(1,2,1,2,2,1,2,1),
                         c(1,2,1,1,1,1,2,1),
                         c(1,1,2,2,2,2,2,1),
                         c(1,1,1,2,2,1,2,1),
                         c(1,1,1,1,1,1,2,1))

###################################
aggSimOutcome<-c()
for (j in 1:dim(policyNonWorkList)[1]){
  policyNonWork<-policyNonWorkList[j,]

  msaList=c("5600","1600","6920","3760")
  
  aggOut<-c()
  for (m in msaList){
    # load calibrated parameters for this location
    par<-calibratedPar(m, Generic)
      
    ### msa policy and date
    msaPD <- loadPolicyDates(m)
      
    # initial condition
    initNumIperType<<-par$I0
      
    out<-c()
    npolicy<-dim(policyCombo)[1]
    for (i in 1:npolicy){
      print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
      print(paste("!! Running SIR for policy combo", i,"/", npolicy, "......!!"))
      
      ## sir across phases
      policy_i <- policyCombo[i,]
      
      policy_i[policy_i==0]<-4
      uni_policy <- unique(policy_i)
      nt <- length(uni_policy)
      policy<-c()
      ttt<-c(0)
      for (t in 1:nt){
        # combine work policies and the rest of policies
        policy<-c(policy, policyTagString(matrix(c(uni_policy[t],policyNonWork),
                                                 1,length(policyNonWork)+1)))
        
        # define timing
        ttt<-c(ttt, ndays * max(which(policy_i==uni_policy[t])))
      }
      TTT<<-ttt
    
      sim_i<-runSir(m, policy, par)
      outstats_i<-calOutcome(sim_i, NA)
      eval1 <- largestEigenvalue(Cmat)
      print(paste(" largest eigenvalue of contact matrices:", 
                  format(eval1,digits=4) )) 
      
      ## collect outputs
      out<-rbind(out,c(outstats_i,eval1,policyNonWork))
    }
    
    ## job loss vs deaths relative to do nothing
    out2plot<-cbind(out[2:npolicy,3]/out[1,3], out[2:npolicy,1]/out[1,1])
    plot(out2plot[,1],out2plot[,2])
    
    ## round change in employment, and find minimized deaths
    # if (max(policyNonWork[1],policyNonWork[2])==1){
    #   out2plot<-cbind(out2plot,round(out2plot[,2]*10)/10)
    # }else{
    #   out2plot<-cbind(out2plot,round(out2plot[,2]*5)/5)
    # }
    # out2plot<-cbind(out2plot,round(out2plot[,2]*10)/10)
    # colnames(out2plot)<-c("deaths","emploss","emplossRound")
    # frontier <- as.data.frame(out2plot) %>% group_by(emplossRound) %>% summarise(deaths = min(deaths))
    # 
    # ### fit a curve for the frontier
    # X<-cbind(1, frontier$emplossRound, frontier$emplossRound^2)
    # y<-X%*%solve(t(X)%*%X, t(X)%*%frontier$deaths)
    # 
    # plot(frontier$emplossRound, frontier$deaths, type="l")
    # lines(frontier$emplossRound, y, type="l", col="red")
    # 
    ### aggregate results
    # we stack them vertically, because the results have different dimension across MSAs
    aggOut<-rbind(aggOut, cbind(as.numeric(m)*rep(1,npolicy-1), out2plot))
    # aggregate raw output
    aggSimOutcome<-rbind(aggSimOutcome, 
                         cbind(as.numeric(m)*rep(1,npolicy), 
                         out, policyCombo))
  }
  
  
  # compare frontier across MSAs -----------------------------------------------------
  plotFrontier <- function(fn, aggOut) {
    
    ### plot frontier for each msa - non-smoothed
    if (fn!="")  fn1 <- file.path(outPath, "figure", fn)
    if (fn!="")  png(fn1)
    msa<-aggOut[,1]
    
    # flip sign so higher number is more life/employment relative to first policy 
    aggOut[,c(2,3)]<--aggOut[,c(2,3)]+1
    plot( aggOut[,2], aggOut[,3], type="l", col="white",
          ylab=paste("Life saved relative to ",Nweek,"-week no policy",sep=""),
          xlab=paste("Employment saved relative to ",Nweek,"-week no policy",sep=""))
    lines(aggOut[msa==5600,2], aggOut[msa==5600,3], type="l", col="black")
    lines(aggOut[msa==1600,2], aggOut[msa==1600,3], type="l", col="red")
    lines(aggOut[msa==6920,2], aggOut[msa==6920,3], type="l", col="blue")
    lines(aggOut[msa==3760,2], aggOut[msa==3760,3], type="l", col="orange")
    
    legend("bottomleft",
           legend=c("NYC", "Chicago", "Sacramento", "Kansas City"),
           col=c("black", "red", "blue", "orange"),
           lty=c(1,1,1,1), 
           horiz=F,lwd=1.3,bty="n",cex=1)
    
    if (fn!=""){
      dev.off()
      print(paste("  saved plot:",fn1))
    }
  }
  
  
  ### plot fontier
  par(mfrow=c(1,2))
  # plotFrontier("", aggOut, 1) 
  plotFrontier("", aggOut) 
  
  
  # ## plot and save fontier
  fname <- paste("np",nb,"_ndays",ndays,"_lowI0_school",
                 policyNonWork[1],"_neighbor",policyNonWork[2],
                 "_poi",policyNonWork[3],policyNonWork[4],policyNonWork[5],policyNonWork[6],sep="")
  # plotFrontier(paste("frontier_smooth_",fname,".png",sep=""), aggOut, 1)
  plotFrontier(paste("frontier_raw_"   ,fname,".png",sep=""), aggOut)

}


colnames(aggSimOutcome)<-c("msa", 
                           "deaths","cases","emploss","hospital","icu","population","maxeval_lastt",
                           "S","N","B","R","P","F","E","M",
                           paste("W_t",seq(1,nb),sep=""))

fn <- file.path(outPath, "csv", 
                paste("frontier_np",nb,"_ndays",ndays,"_lowI0.csv", sep=""))
checkWrite(fn, aggSimOutcome, "SIR outputs") 

