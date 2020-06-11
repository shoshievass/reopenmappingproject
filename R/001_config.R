##############################################################################
##############################################################################
######                                                                  ######
######                 COVID Project Configuration                      ######
######                                                                  ######
##############################################################################
##############################################################################


#####################################
## load programs
require("deSolve")
require("latex2exp")


#####################################
### hardcode global variables
#####################################

#industry to plot
naics2plot<<-c(31,42,44,52,54,62,72)
naicsName<<-c("Manufacturing*","Wholesale*","Retail","Finance","Professional & IT","Healthcare*","Accommodation")

#essential naics2
ESSENTIAL<<-c(11,21,22,31,42,48,62,92)

#compartments
COMPART   <<-c("S","E","Ia","Ins","Ihc","Rq","Rqd","Rnq","D")

#### set of colors
mycol <<- c(t_col("red",     perc = 0),
            t_col("orange",  perc = 0),
            t_col("yellow2", perc = 0),
            t_col("green",   perc = 0),
            t_col("dimgray", perc = 0),
            t_col("blue",    perc = 0),
            t_col("purple",  perc = 0),
            t_col("red",     perc = 40),
            t_col("orange",  perc = 40),
            t_col("yellow2", perc = 40),
            t_col("green",   perc = 40),
            t_col("dimgray", perc = 40),
            t_col("blue",    perc = 40),
            t_col("purple",  perc = 40))
mycolLength<<-7


#####################################
# set up for run
#####################################
### scenarios/place

# NYC, Chicago, Sacramento
# msa  "5600", "1600", "6920"


## all scenarios
contactList<<-c("_regular","_socialdistance","_alternateschoolwork","_shutold", "_wfhreopenschool","_cautious")

## reference policies
refPolicy   <-expand.grid(c("_regular"),c("_socialdistance"),c("_regular"))

## combinations
policyComboFull<- expand.grid(c("_regular"),contactList,contactList)


adhocPolicy   <-expand.grid(c("_regular"),c("_socialdistance"),contactList)


policyCombo<<-rbind(refPolicy,adhocPolicy)

## time period
TTT <<-c(0,15,75,150)
placeList<-c("5600", "1600", "6920")
# placeList<-c( "1600")
# placeList<-c("6920")

#####################################
# Versions
#####################################

# ver<<-"_combo"
# ver<<-"_comboReplica"
ver<<-"_comboReplicav0"

# ver<<-"_comboReplicav0_beta0"
# ver<<-"_comboReplicav0_beta1"
# ver<<-"_comboReplicav0_beta09"
# ver<<-"_comboReplicav0_beta075"
# ver<<-"_comboReplicav0_beta05"
# ver<<-"_comboReplicav0_beta025"
# # 
# ver<<-"_comboReplica_psi08"
# ver<<-"_comboReplica_psi09"
# ver<<-"_comboReplica_psi06_8"
# ver<<-"_comboReplica_psi05_9"





# fix beta as counterfactual, 0 for varying beta, 1 for beta0 and 2 for beta1
fixBETA<<-0

#beta scaler
scalBETA<<-1
#90% 75% 50% 25%

#psi scaler
scalPSI<<-1
#scaled to 0.8/0.7, 0.9/0.7


#####################################
# SIR parameters
#####################################


beta    <- 0.015
betaH   <<- 0.1         # transmission from infected to healthcare workers
epsilon <- 1/4         # incubation period
EPSILON <<- c(0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.33,0.33,0.35,0.35,0.35,0.35,0.35)         # incubation period
gamma   <- 1/22        # 1/symptom to recovery
gammaD  <- 1/20        # 1/symptom to death
tauAvg  <- 1/5         # infected to symptom or never symptom
TAU     <<-c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.25,0.25,0.33,0.33,0.33,0.33,0.33)
psiAvg  <- 0.4         # prob of showing symptom (average). 
# Simon: 50 to 75 percent of infected cases are asymptomatic based on a study in Italy
# PSI     <<- c(0, 0.01, 0.1, 0.3, 0.3, 0.4, 0.45, 0.55, 0.5, 0.6, 0.6, 0.6, 0.6, 0.6) # prob of showing symptom (age X sick )
PSI <<- rep(0.7, 14)*scalPSI
# PSI     <<- c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8) 
# PSI     <<- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.7, 0.7, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9)

eta     <<- 1/14        # time spent in quarantine
mort    <- c(0.0001,0.00015,0.0001,0.00015,0.001,0.002,0.005,0.01,0.025,0.04,0.03,0.05,0.06,0.07) # prob death | infected (age X sick)
DELTAhc <<-       mort*gammaD/pmax(1e-10,PSI)         # death rate (age X sick)
# acemoglu et al: death rate 0.1%, 1% and 6% for <49, 50-64, 65+
#DELTAa  <<- 0.001 *mort*   TAU/pmax(1e-10,PSI)        # adjustment of death rate of asymoptomatic people (at this level, similar to influenza https://www.cdc.gov/nchs/data/nvsr/nvsr68/nvsr68_09-508.pdf)
GAMMA <<- gammaD - DELTAhc # recovery rate to account for variation in death
vparameters0 <<- c(gamma=gamma,beta=beta,betaH=betaH,eta=eta)

## unique age X sick type: age*10 + sick
typeAgeSick <<-as.matrix(expand.grid(0:1,0:6)) %*% c(1,10)

#wtd average duration of infected
infectDuration<<-psiAvg * (1/tauAvg) + (1-psiAvg) * (1/tauAvg + 1/gamma)

# initial condition: number of people in I^A per type
initNumIperType<<-1

#



#data version
datv<<-"_replica"


### save results/plots?
saveSIR<<-1
outpSIR<<-1


### for plots
psize<<-1.25


