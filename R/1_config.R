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
packages <- c("deSolve","plyr","dplyr","tidyr","pracma")
newPackages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)

library(deSolve)
library(plyr)
library(dplyr)
library(tidyr)
library(pracma)

#####################################
## check input files
#####################################

## user input file names

# data for fitting death
deathData <<- file.path(dataPath,"covid_case_death_jh.csv")

# SEIR model parameters
seirParm  <<- file.path(parmPath,"params.csv")

# contact matrix file names start with
ctMatData <<- "C_msa"

# inputs for grid search calibration
gsParm <<- file.path(parmPath, "gridsearch.csv")

# calibrated parameter name
caliParm  <<- "calibrated_parm_msa"



#####################################
# policy and locations
#####################################
### scenarios/place


## reference policies: NP, EO, NP(M2), NP(M2)
refPhase1 <<-"_W4-S3-N3-E2-M3"
refPhase2 <<-"_W1-S1-N1-E2-M2"
refPhase3 <<-"_W4-S3-N3-E2-M2"
refPhase4 <<-"_W4-S3-N3-E2-M2"


## contact definition for 
# W(work), S(school), N(neighbor) contacts, whether we quarantine E(elderly), M(mask), although M no impact on contact

# all reopening policies
contactPolicy<<-expand.grid(1:4,1:3,1:3,1:2,1:3)

#NP, EO, CR, AS, WFH, 60+
#reduced beta, normal beta, even lower beta
reopenPolicy<<-rbind(c(4,3,3,2,2),c(1,1,1,2,2),c(4,3,1,2,2),c(3,2,1,2,2),c(2,3,1,2,2),c(4,3,1,1,2))
policyCombo<<-genPolicy(reopenPolicy,refPhase1,refPhase2,refPhase3,refPhase4)


#####################################
# hard code parameters
#####################################


## MSAs
# NYC, Chicago, Sacramento, Houston, Kansas City
# msa  "5600", "1600", "6920","3360", "3760"
# msaList<<-c("1600","6920","3760")
msaList<<-c("5600")


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

#input (data,contact matrix) version
datv<<-""

### save results/plots?
outputSIR<<-1

## source of input for contact matrix
Csource<<-list(msa5600="fred",msa7240="fred", msa3360="fred", msa1920="fred",
               msa1600="replica",msa6920="replica", msa3760="replica")



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

