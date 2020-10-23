

# load policy and dates for this msa
P <- checkLoad(policyParm) 
P <- P[P$MSA==m,]


# key outputs
t<-max(TTT)

######## sim_i is generated from running 5_seir.R ########
simRun <- simRef
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


## first average outcomes across types we want to keep track
types_temp <- types %>% 
  group_by(age, naics, sick, wfh) %>% 
  summarise(deathRate  = sum(deathRate*n)/sum(n), 
            empLoss = sum(empLoss*n)/sum(n),
            n = sum(n)) 
## people with employed naics 
types_temp$employed <- 1*(types_temp$naics>0)

## merge type specific outcome with demo distribution
## not all types exit in contact matrix
Demo1 <- merge(x = Demo, y = types_temp)

stopifnot(sum(is.na(Demo1))==0)

## accounting for race/income based outcome
# pr_type_ineq = prob ( age, naics, sick, wfh | race income region)
Demo2 <- Demo1 %>% 
  group_by(race, income) %>% 
  summarise(totalW = sum(pr_type_ineq), 
            deathRate = sum(pr_type_ineq * deathRate)/sum(pr_type_ineq),
            empLoss   = sum(pr_type_ineq * empLoss * employed)/sum(pr_type_ineq * employed),
            age       = sum(pr_type_ineq * raw_age)/sum(pr_type_ineq),
            essential = sum(pr_type_ineq * essential)/sum(pr_type_ineq),
            N = sum(N))

Demo2$pr_race_income <- Demo2$N*100/sum(Demo2$N)
as.matrix(Demo2)
