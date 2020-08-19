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

# demo accounting
demoAcct  <<- file.path(dataPath, "demo_accounting.csv")

# data for fitting death
deathData <<- paste(dataPath,"covid_case_death_jh.csv",sep="/")

# SEIR model parameters
seirParm  <<- paste(parmPath,"params.csv",sep="/")

# user input contact matrix file names start with
ctMatRaw <<- "contact_poi_msa"
# ctMatRaw <<- "contact_msa"

# contact matrix file names start with
ctMatData <<- "C_msa"

# msa specific policy scenarios and dates
# policyParm <<- file.path(parmPath, "msa_policy_scenarios_dates.csv")
policyParm <<- file.path(parmPath, "msa_poi_policy_scenarios_dates.csv")

# inputs for grid search calibration
gsParm <<- file.path(parmPath, "gridsearch.csv")

# calibrated parameter name
caliParm  <<- "calibrated_parm_msa"



#####################################
# policy and locations
#####################################
### scenarios/place

## contact definition for 
# W(work), S(school), N(neighbor) contacts, whether we quarantine E(elderly), M(mask), although M no impact on contact
policyLetterCode <<- c("W", "S", "N", "E", "M")
policyPOILetterCode <<- c("W", "S", "N", "B", "R", "P", "F", "E", "M")



## reference policies
refPhase1 <<-"_W4-S3-N3-E2-M3"
refPhase2 <<-"_W1-S1-N1-E2-M2"
refPhase3 <<-"_W4-S3-N1-E2-M2"
refPolicy   <-c(refPhase1,refPhase2,refPhase3)


# all combinations of reopening policies including POI policies
contactPolicyPOI<<-expand.grid(1:4,1:3,1:3,1:2,1:2,1:2,1:2,1:2,1:4)
poiLevel <<- c("restaurants", "retail", "services", "entertainment")


### key policies considered
#NP, EO, CR, AS, WFH, 60+
reopenPolicy<<-rbind(c(4,3,3,2,2,2,2,2,1),c(1,1,1,1,1,1,1,2,1),c(4,3,1,2,2,2,2,2,1),c(3,2,1,2,2,2,2,2,1),c(2,3,1,2,2,2,2,2,1),c(4,3,1,2,2,2,2,1,1))

### all possible combo for reopen policy in the final phase
# reopenPolicy<<-contactPolicyPOI



## MSAs
# NYC, Chicago, Sacramento, Houston, Kansas City
msaList<<-c("5600","1600","6920","3360","3760")
msaList<<-c("5600")
#####################################
# key global variables
#####################################


## age number of 60
age60<<-4

## school age (age>=school age is definitely not attending schools)
ageSchool<<-3

## naics code for healthcare
heathNAICS<<-62

## source of input for contact matrix
Csource<<-list(msa5600="fred",msa1600="replica",msa6920="replica")

#####################################
# input/output versions
#####################################

#output version
ver <<-"_combo"

#input (data,contact matrix) version
datv<<-""

### save detailed seir compartment X type X time level results and plots for internal checking?
outputSIR<<-1

# fix beta as counterfactual, 0 for varying beta, 1 for beta1 and 2 for beta2
fixBETA  <<-0

#beta scale factor to test sensitivity
scalBETA <<-1





#####################################
# load and define SIR parameters
#####################################
#compartments
COMPART   <<-c("S","E","Ia","Ins","Ihc","Rq","Rqd","Rnq","D")

# load user defined SIR parameters
PAR <-checkLoad(seirParm)

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
DELTAhc <<-(mort/psi)   * gammaD     # death rate 
GAMMA   <<-gammaD - DELTAhc         # recovery rate to account for variation in death, so total transition rate out of infected is kept at gammaD

## unique age X sick type: age*10 + sick
typeAgeSick <<-as.matrix(PAR$age*10+PAR$sick)

#wtd average duration of infected (not in the SIR model, a scaling factor in calibration exercise)
# infectDuration<<-min(PAR$infectionDuration)
infectDuration<<-psi * (1/mean(TAU)) + (1-psi) * (1/mean(TAU) + 1/gamma)
tEtoD <<-round(1/mean(EPSILON) + 1/mean(TAU) +  1/gammaD)
  
  
# initial condition: number of people in I^A per type
initNumIperType<<-1


#####################################
### global variables for plots
#####################################

#industry to plot
naics2plot<<-c(31,42,44,52,54,62,72)
naicsName<<-c("Manufacturing*","Wholesale*","Retail","Finance","Professional & IT","Healthcare*","Accommodation")


# setting for plots
gen_col()


# clear variables
rm(PAR, mort, gamma, gammaD, betaH, eta, psi)

