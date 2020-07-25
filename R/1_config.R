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
## check folders exist
#####################################

# check if directories exist
dir.create(tempPath, showWarnings = FALSE)
dir.create(contactMatrixPath, showWarnings = FALSE)
dir.create(calibratedParPath, showWarnings = FALSE)
dir.create(parmPath, showWarnings = FALSE)
dir.create(outPath, showWarnings = FALSE)
dir.create(figPath, showWarnings = FALSE)



#####################################
## input file names
#####################################


# msa type adjustment
msaType  <<- file.path(dataPath, "msa_type.csv")

# data for fitting death
deathData <<- file.path(dataPath,"covid_case_death_jh.csv")

# SEIR model parameters
seirParm  <<- file.path(parmPath,"seir_parameters.csv")

# contact matrix file names start with
ctMatData <<- "C_msa"

# msa specific policy scenarios and dates
policyParm <<- file.path(parmPath, "msa_policy_scenarios_dates.csv")

# inputs for grid search calibration
gsParm <<- file.path(parmPath, "gridsearch.csv")

# calibrated parameter name
caliParm  <<- "calibrated_parm_msa"



#####################################
# policy and locations
#####################################

## contact definition for 
# W(work), S(school), N(neighbor) contacts, whether we quarantine E(elderly), M(mask), although M no impact on contact

policyLetterCode <<- c("W", "S", "N", "E", "M")

#default E and M for initial vs subsequent phases in calibration
policyCaliEM <<- c("-E2-M4","-E2-M1")


# all combinations of reopening policies
contactPolicy<<-expand.grid(1:4,1:3,1:3,1:2,1:4)

### key policies considered
#NP, EO, CR, AS, WFH, 60+
#reduced beta, normal beta, even lower beta
reopenPolicy<<-rbind(c(4,3,3,2,1),c(1,1,1,2,1),c(4,3,1,2,1),c(3,2,1,2,1),c(2,3,1,2,1),c(4,3,1,1,1))


### all possible combo for reopen policy in the final phase
reopenPolicy<<-contactPolicy

# generate results for multiple reference policies, 
# this is for area plot
genRef4Area<<-0

## MSAs
# NYC, Chicago, Sacramento, Houston, Kansas City
# msa  "5600", "1600", "6920","3360", "3760"
 msaList<<-c("5600","6920","3760")
#msaList<<-c("1600")

#####################################
# key global variables
#####################################

## age number of 60
age60<<-4

## naics code for healthcare
heathNAICS<<-62

## fix beta as counterfactual, 0 for varying beta, 1 for beta1 and 2 for beta2
fixBETA  <<-0

## beta scale factor to test sensitivity
scalBETA <<-1

## t0, starting time for SIR model
TNAUGHT <<- as.Date(unique("3/5/2020"), "%m/%d/%Y")

#####################################
# input/output versions
#####################################

#output version
verTag <<-"_combo"


### save detailed sir compartment X type X time level results and plots for internal checking?
outputSIR<<-0


#####################################
# load and define SEIR parameters
#####################################
#compartments
COMPART   <<-c("S","E","Ia","Ins","Ihc","Rq","Rqd","Rnq","D")

# load user defined SIR parameters
PAR <-checkLoad(seirParm)

#constant parameters
betaH <-min(PAR$beta_healthcare)
gamma <-1/min(PAR$gamma_inv)
gammaD<-1/min(PAR$gammaD_inv)
psi   <-min(PAR$psi)
eta   <-1/min(PAR$eta_inv)
vparameters0 <<- c(gamma=gamma,betaH=betaH,eta=eta,psi=psi)

#type specific parameters
EPSILON<<-PAR$epsilon
TAU    <<-PAR$tau

#input mortality conditional on infected, transform into transition rate conditional on symptomatic
mort     <-PAR$mort
DELTAhc <<-mort*gammaD/psi     # death rate 
GAMMA   <<-gammaD - DELTAhc    # recovery rate to account for variation in death, so total transition rate out of infected is kept at gammaD

## unique age X sick type: age*10 + sick
typeAgeSick <<-as.matrix(PAR$age*10+PAR$sick)

#wtd average duration of infected (not in the SIR model, a scaling factor in calibration exercise)
infectDuration<<-psi * (1/mean(TAU)) + (1-psi) * (1/mean(TAU) + 1/gamma)
  
# initial condition: number of people in I^A per type
initNumIperType<<-1


#########################################################
### global variables for plots for internal testing
#########################################################

#industry to plot
naics2plot<<-c(62,31,42,44,52,54,72)
naicsName<<-c("Healthcare*","Manufacturing*","Wholesale*","Retail","Finance","Professional & IT","Accommodation")


# setting for plot color
gen_col()


# clear variables
rm(PAR, mort, gamma, gammaD, betaH, eta, psi)

