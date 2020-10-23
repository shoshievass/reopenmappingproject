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
#do nothing, lock down , wfh
work<-c(0,1,2)

## number of periods and number of days in between policies
nb<-8

ndays<-7

## possible grouping of breaks
policies<-permutations(length(work), nb, v=work, repeats.allowed=TRUE)

## different work policies
policyCombo<-policies[apply(t(diff(t(policies))), 1, FUN=min) ==0,]

## the rest of policies
# W(work), 
# S(school), N(neighbor), POI-resturants, POI-retail, POI-personal service, POI-entertainment, eldery, mask=1
# HY: always set mask to 1, we can test different policies for school, etc and see what frontier looks like
policyNonWork<-c(1,2,1,1,1,1,2,1)
###################################


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
  
  for (i in 1:dim(policyCombo)[1]){
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(paste("!! Running SIR for policy combo", i,"/", dim(policyCombo)[1], "......!!"))
    
    ## sir across phases
    policy_i <- policyCombo[i,]
    
    policy_i[policy_i==0]<-4
    uni_policy <- unique(policy_i)
    nt <- length(uni_policy)
    policy<-c()
    ttt<-c(0)
    for (t in 1:nt){
      # combine work policies and the rest of policies
      policy<-c(policy, policyTagString(matrix(c(uni_policy[t],policyNonWork),1,length(policyNonWork)+1)))
      # define timing
      ttt<-c(ttt, ndays * max(which(policy_i==uni_policy[t])))
    }
    TTT<<-ttt
  
    sim_i<-runSir(m, policy, par)
    outstats_i<-calOutcome(sim_i, NA)
    
    ## collect outputs
    out<-rbind(out,outstats_i)
  }
  
  ## job loss vs deaths relative to do nothing
  out2plot<-cbind(out[,1]/out[1,1],out[,3]/out[1,3])
  plot(-out2plot[,1],-out2plot[,2])
  
  ## round change in employment, and find minimized deaths
  out2plot<-cbind(out2plot,round(out2plot[,2]*10)/10)
  colnames(out2plot)<-c("deaths","emploss","emplossRound")
  frontier <- as.data.frame(out2plot) %>% group_by(emplossRound) %>% summarise(deaths = min(deaths))
  
  ### fit a curve for the frontier
  X<-cbind(1, frontier$emplossRound, frontier$emplossRound^2)
  y<-X%*%solve(t(X)%*%X, t(X)%*%frontier$deaths)

  plot(frontier$emplossRound, frontier$deaths, type="l")
  lines(frontier$emplossRound, y, type="l", col="red")
  
  ### aggregate results
  # we stack them vertically, because the results have different dimension across MSAs
  aggOut<-rbind(aggOut, cbind(as.numeric(m)*rep(1,length(y)), frontier$deaths, y, X[,2]))
}



### plot frontier for each msa - smoothed
fn <- file.path(outPath, "figure", "frontier_v1.pdf")
pdf(fn)
msa<-aggOut[,1]
plot( -aggOut[,4]+1, -aggOut[,3]+1, type="l", col="white",
      ylab="Life saved relative to 8-week no policy",xlab="Employment saved relative to 8-week no policy")
lines(-aggOut[msa==5600,4]+1, -aggOut[msa==5600,3]+1, type="l", col="dimgray")
lines(-aggOut[msa==1600,4]+1, -aggOut[msa==1600,3]+1, type="l", col="red")
lines(-aggOut[msa==6920,4]+1, -aggOut[msa==6920,3]+1, type="l", col="blue")
lines(-aggOut[msa==3760,4]+1, -aggOut[msa==3760,3]+1, type="l", col="orange")

legend("bottomleft",
       legend=c("NYC", "Chicago", "Sacramento", "Kansas City"),
       col=c("dimgray", "red", "blue", "orange"),
       lty=c(1,1,1,1), 
       horiz=F,lwd=1.3,bty="n",cex=1)

dev.off()

### plot frontier for each msa - non-smoothed
fn <- file.path(outPath, "figure", "frontier_v2.pdf")
pdf(fn)
msa<-aggOut[,1]
plot( -aggOut[,4]+1, -aggOut[,2]+1, type="l", col="white",
      ylab="Life saved relative to 8-week no policy",xlab="Employment saved relative to 8-week no policy")
lines(-aggOut[msa==5600,4]+1, -aggOut[msa==5600,2]+1, type="l", col="dimgray")
lines(-aggOut[msa==1600,4]+1, -aggOut[msa==1600,2]+1, type="l", col="red")
lines(-aggOut[msa==6920,4]+1, -aggOut[msa==6920,2]+1, type="l", col="blue")
lines(-aggOut[msa==3760,4]+1, -aggOut[msa==3760,2]+1, type="l", col="orange")

legend("bottomleft",
       legend=c("NYC", "Chicago", "Sacramento", "Kansas City"),
       col=c("dimgray", "red", "blue", "orange"),
       lty=c(1,1,1,1), 
       horiz=F,lwd=1.3,bty="n",cex=1)

dev.off()

