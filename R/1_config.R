##############################################################################
##############################################################################
######                                                                  ######
######        reopen mapping project -- configuration                   ######
######                                                                  ######
##############################################################################
##############################################################################


#####################################
## check pakges
#####################################
packages <- c("deSolve")
newPackages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)

library(deSolve)



#####################################
## check files
#####################################
setwd(parmPath)

## input parameters
if (max(list.files()=="params.csv")==0){
  stop(paste("missing input parameters files: ", parmPath, "/params.csv", sep=""))
}

## industry reopen status under each policies
if (max(list.files()=="naics2essentialpolicy.csv")==0){
  stop(paste("missing industry open status files: ", parmPath, "/naics2essentialpolicy.csv", sep=""))
}

setwd(dataPath)

## death counts used for calibration
if (max(list.files()=="covid_case_death_jh.csv")==0){
  stop(paste("missing deaths count from Johns Hopkins: ", parmPath, "/covid_case_death_jh.csv", sep=""))
}




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

### setting for plots
psize<<-1.25

#####################################
# set up for run
#####################################
### scenarios/place

## all scenarios
contactList<<-c("_regular","_socialdistance","_alternateschoolwork","_shutold", "_wfhreopenschool","_cautious")

## reference policies
refPolicy   <-expand.grid(c("_regular"),c("_socialdistance"),c("_cautious"))

## combinations
policyFull <- expand.grid(c("_regular"),contactList,contactList)


policyFull <- expand.grid(c("_regular"),c("_socialdistance"),contactList)


policyCombo<<-rbind(refPolicy)

## time period
TTT <<-c(0,15,75,150)

## locations
# NYC, Chicago, Sacramento
# msa  "5600", "1600", "6920"
placeList<-c("5600", "1600", "6920")
placeList<-c("1600")

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


# fix beta as counterfactual, 0 for varying beta, 1 for beta0 and 2 for beta1
fixBETA<<-0

#beta scaler
scalBETA<<-1
#90% 75% 50% 25%

#####################################
# load SIR parameters
#####################################
PAR <-read.csv(paste(parmPath,"params.csv",sep="/"),header=TRUE)

#constant parameters
beta  <-min(PAR$beta) #beta will be calibrated later on
betaH <-min(PAR$betaH)
gamma <-1/min(PAR$gamma_inv)
gammaD<-1/min(PAR$gammaD_inv)
psi   <-min(PAR$psi)
eta   <-1/min(PAR$eta_inv)
vparameters0 <<- c(gamma=gamma,beta=beta,betaH=betaH,eta=eta,psi=psi)

#type specific parameters
EPSILON<<-PAR$epsilon
TAU    <<-PAR$tau

#input mortality conditional on infected, transform into transition rate
mort     <-PAR$mort
DELTAhc <<-mort*gammaD/psi     # death rate 
GAMMA   <<-gammaD - DELTAhc    # recovery rate to account for variation in death, so total transition rate out of infected is kept at gammaD

## unique age X sick type: age*10 + sick
typeAgeSick <<-as.matrix(PAR$age*10+PAR$sick)

#wtd average duration of infected
infectDuration<<-psi * (1/mean(TAU)) + (1-psi) * (1/mean(TAU) + 1/gamma)

# initial condition: number of people in I^A per type
initNumIperType<<-1


#data version
datv<<-"_replica"


### save results/plots?
outputSIR<<-1





