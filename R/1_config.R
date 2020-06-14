##############################################################################
##############################################################################
######                                                                  ######
######        reopen mapping project -- configuration                   ######
######                                                                  ######
##############################################################################
##############################################################################


#####################################
## check packages
#####################################
packages <- c("deSolve","plyr","dplyr","tidyr")
newPackages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)

library(deSolve)
library(plyr)
library(dplyr)
library(tidyr)


#####################################
## check input files
#####################################

## user input file names

# data for fitting death
deathData <<- "covid_case_death_jh.csv"

# industry open/close definition in each policy
indPlData <<- "naics2essentialpolicy.csv"

# SEIR model parameters
seirParm <<- "params.csv"

# contact matrix file names start with
ctMatData <<- "C_msa"

# calibrated parameter name
caliParm  <<- "calibrated_parm_msa"


## input parameters
setwd(parmPath)
if (max(list.files()==seirParm)==0){
  stop(paste("missing input parameters files: ", parmPath, seirParm, sep="/"))
}

## industry reopen status under each policies
if (max(list.files()==indPlData)==0){
  stop(paste("missing industry open status files: ", parmPath, indPlData, sep="/"))
}

## death counts used for calibration
setwd(dataPath)
if (max(list.files()==deathData)==0){
  stop(paste("missing deaths count from Johns Hopkins: ", parmPath, deathData, sep="/"))
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



#####################################
# policy and locations
#####################################
### scenarios/place

## all scenarios
policyList<<-c("_NP","_EO","_AS","_I60","_WFH","_CR")

## reference policies
refPolicy   <-expand.grid(c("_NP"),c("_EO"),c("_CR"))

## combinations
policyFull <- expand.grid(c("_NP"),c("_EO"),policyList)


policyCombo<<-rbind(refPolicy)

## time period of each phases with different policies
TTT <<-c(0,15,75,150)

## locations
# NYC, Chicago, Sacramento
# msa  "5600", "1600", "6920"
msaList<<-c("5600")

#####################################
# Versions
#####################################

#output version
ver<<-"_combo"

#input (data/contact matrix/parameter) version
datv<<-""

### save results/plots?
outputSIR<<-1


# fix beta as counterfactual, 0 for varying beta, 1 for beta1 and 2 for beta2
fixBETA<<-0

#beta scale factor to test sensitivity
scalBETA<<-1





#####################################
# load SIR parameters
#####################################
PAR <-read.csv(paste(parmPath,seirParm,sep="/"),header=TRUE)

#constant parameters
betaH <-min(PAR$betaH)
gamma <-1/min(PAR$gamma_inv)
gammaD<-1/min(PAR$gammaD_inv)
psi   <-min(PAR$psi)
eta   <-1/min(PAR$eta_inv)
vparameters0 <<- c(gamma=gamma,betaH=betaH,eta=eta,psi=psi)

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
infectDuration<<-min(PAR$infectionDuration)
# infectDuration<<-psi * (1/mean(TAU)) + (1-psi) * (1/mean(TAU) + 1/gamma)

  
# initial condition: number of people in I^A per type
initNumIperType<<-1

# setting for plots
gen_col()


# clear variables
rm(PAR, mort, gamma, gammaD, betaH, eta, psi)

