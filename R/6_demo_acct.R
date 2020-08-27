

# load policy and dates for this msa
P <- checkLoad(policyParm) 
P <- P[P$MSA==m,]


# demo graphic accounting inputs
Demo <- checkLoad(demoAcct) 
Demo <- Demo[Demo$region==P$Division,]

# key outputs
t<-max(TTT)

######## sim_i is generated from running 5_seir.R ########
simRun <- sim_i
coln <- colnames(simRun)


### cumulative death rate
deathRate       <- simRun[,pickState("D",coln)]
types$deathRate <- 1e5*colSums(deathRate[t+1,])/types$n


## includes compartment that are working
s<-simRun[,pickState("S",coln)]
for (i in c("E","Ia","Ins","Rqd","Rnq")){
  s<-s+simRun[,pickState(i,coln)]
}
## consider whether these individuals are actively working
active<-s * simRun[,pickState("active",coln)]

## compute total employment loss in days
types$empLoss <- t*(types$naics>0) - colSums(active[2:(t+1),])/types$n


## first average across types
types_temp <- types %>% 
  group_by(age, naics, work_poi, sick, wfh, active_emp) %>% 
  summarise(deathRate = mean(deathRate), 
            empLoss   = mean(empLoss),
            n = sum(n)) 


## merge type specific outcome with demo distribution
## not all types exit in ACS data, 
Demo1 <- merge(x = Demo, y = types_temp)

stopifnot(sum(is.na(Demo1))==0)
Demo2 <- Demo1 %>% 
  group_by(race, income) %>% 
  summarise(totalW = sum(pr_type_ineq), 
            deathRate = sum(pr_type_ineq * deathRate)/sum(pr_type_ineq),
            empLoss   = sum(pr_type_ineq * empLoss)/sum(pr_type_ineq * active_emp),
            N = sum(N))

Demo2$pr_race_income <- Demo2$N*100/sum(Demo2$N)
as.matrix(Demo2)
