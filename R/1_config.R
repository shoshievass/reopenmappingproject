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
packages <- c("deSolve","plyr","dplyr","tidyr","gtools","latex2exp")
newPackages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)

library(deSolve)
library(plyr)
library(dplyr)
library(tidyr)
library(gtools)
library(latex2exp)


#####################################
## check folders exist
#####################################

# check if directory exist
dir.create(tempPath, showWarnings = FALSE)
dir.create(contactMatrixPath, showWarnings = FALSE)
dir.create(calibratedParPath, showWarnings = FALSE)
dir.create(parmPath, showWarnings = FALSE)
dir.create(outPath, showWarnings = FALSE)


#####################################
## input file names
#####################################



# msa type adjustment
msaType  <<- file.path(dataPath, "msa_type.csv")

# demo accounting
demoAcct  <<- file.path(dataPath, "demo_accounting.csv")

# data for fitting death
deathData <<- paste(dataPath,"covid_case_death_jh.csv",sep="/")

# SEIR model parameters
seirParm  <<- file.path(parmPath,"seir_parameters.csv")

# user input contact matrix file names start with
# version with POI
ctMatRaw <<- "contact_poi_msa"
# version without POI
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
policyLetterCode    <<- c("W", "S", "N", "E", "M")
policyPOILetterCode <<- c("W", "S", "N", "B", "R", "P", "F", "E", "M")


#default E and M for initial vs subsequent phases in calibration
# by default, we do not isolate elderly,
# M1-4 corresponds to beta estimates in each period
policyCaliEM <<- c("-E2-M1","-E2-M2","-E2-M3","-E2-M4")

# all combinations of reopening policies including POI policies
# the final entry is mask, we only run mask in phase I and mask in phase IV
contactPolicyPOI<<-expand.grid(1:4,1:3,1:3,1:2,1:2,1:2,1:2,1:2,c(1,4))
poiLevel <<- c("restaurants", "retail", "services", "entertainment")

### key policies considered (not specifying the mask)
#NP, EO, CR, AS, WFH, 60+ (set of policies in NBER working paper)
# reopenPolicy<<-rbind(c(4,3,3,2,2,2,2,2),
#                      c(1,1,1,1,1,1,1,2),
#                      c(4,3,1,2,2,2,2,2),
#                      c(3,2,1,2,2,2,2,2),
#                      c(2,3,1,2,2,2,2,2),
#                      c(4,3,1,2,2,2,2,1))

#NP, EO, CR, NR, AS, WFH, 60+ (new set of policies, to be discussed)
reopenPolicy<<-rbind(c(4,3,3,2,2,2,2,2),
                     c(1,1,1,1,1,1,1,2),
                     c(4,1,2,2,2,2,1,2),
                     c(4,1,3,2,2,2,2,2),
                     c(3,2,1,2,2,2,2,2),
                     c(2,1,1,1,2,1,1,2),
                     c(4,3,1,2,2,2,2,1))

reopenPolicy<<-rbind(c(4,3,3,2,2,2,2,2))

# generate results for multiple reference policies, 
# this is for area plot in the NBER working paper
genRef4Area<<-0

# do we run interaction of 'reopenPolicy' across all phases (=1) 
# or keep first 3 phases as the reference policy and run all possible reopen policies 'contactPolicyPOI' in the last phase (=0)?
AllPhases<<-1


## MSAs
# NYC, Chicago, Sacramento, Kansas City
msaList<<-c("5600","1600","6920","3760")
msaList<<-c("5600")


#####################################
# key global variables
#####################################


## age number of 60
age60<<-6

## school age (age>=school age is definitely not attending schools)
ageSchool<<-3

## age group names
ageNames <<-c("5_15", "16_29", "30_39", "40_49", "50_59", "60_69", "70_79", "80")
employNames <<-c("Unemploy", "CanNotWFH", "CanWFH")
raceincomeNames <<-c("BlackLowInc", "BlackMidInc", "BlackHigInc",
                     "HispnLowInc", "HispnMidInc", "HispnHigInc",
                     "WhiteLowInc", "WhiteMidInc", "WhiteHigInc")

## naics code for healthcare
heathNAICS<<-62

## t0, starting time for SIR model
TNAUGHT <<- as.Date(unique("3/5/2020"), "%m/%d/%Y")



#####################################
# input/output versions
#####################################

# output version
verTag <<-"_combo"

# save detailed seir compartment X type X time level results and plots for internal checking?
outputSIR<<-0

# use estimated or generic beta (default 0, 1 and 2 for the compare MSA exercise)
# 0: use MSA specific estimated parameter
# 1: use MSA specific beta and generic initial condition
# 2: use generic parameters (used in our current draft)
Generic<<-0
if (Generic>0){
  verTag <<-paste(verTag,'_par',Generic,sep="")
}


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

# rate at which Ia flow into hospital/icu
# hospital/icu are  (hospital|symptomatic) and  (icu|symptomatic)
# to compute flow from Ia we need to multiply by the prob of symptom (psi) and the rate (TAU)
HOSPITAL <<-PAR$hospital * TAU * psi
ICU      <<-PAR$hospital * TAU * psi

# hospitalization/icu duration
HOSPDAYS<<-min(PAR$hospital_days)
ICUDAYS <<-min(PAR$icu_days)


# input mortality conditional on (hospita|symptomatic) from https://healthcare-in-europe.com/en/news/covid-19-high-mortality-in-hospital-patients.html
# transform into (death|symptomatic) * rate (gammaD)
DELTAhc <<-PAR$death_in_hosp * PAR$hospital * gammaD         # death rate 
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
rm(PAR, gamma, gammaD, betaH, eta, psi)

