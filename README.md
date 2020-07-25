# Read Me for the Reopen Mapping Project

This repository contains code to run estimation and simulation for [Socioeconomic Network Heterogeneity and Pandemic Policy Response](https://www.nber.org/papers/w27374). Directory structure and instructions for use below. Code for the [website](reopenmappingproject.com) is [here](https://github.com/dalek2point3/reopen-website/).


## Directory Contents

### R/
Folder for R code

#### 0_run_all.R
Set up directories, load functions, and run the entire process. 

#### 1_config.R
Users need to make changes here for the specific run. 
This script includes some hardcodes for policies, locations, etc. It also loads SEIR disease parameters.  

It includes the following key sections:

1. Input file names: specify file names of input data and parameters, and intermediate outputs

2. Policy and locations: specify policies and MSAs to run. 

Policies code names have the following convention:

W: work contact

	W1: essential work only;
	W2: remote if possible;
	W3: alternate schedules;
	W4: all in-person
	
S: school contact

	S1: online;
	S2: alternate schedules;
	S3: all in-person
	
N: neighbor contact

	N1: necessities only;
	N2: moderate interaction;
	N3: all in-person
	
E: protect elderly

	E1: isolate 60+; 
	E2: no specific policies
	
M: mask usage and social distancing measures 

	M1: widespread mask usage and social distancing;
	M2: mostly consistent;
	M3: not consistent;
	M4: rare

3. Key global variables: specify age cohort corresponding to 60. Naics code for healthcare (to account for transmission from patients to healthcare workers),

4. Load and define SEIR parameters: load SEIR parameters from /parameter/seir_parameters.csv, and set up some global variables to store these parameters for simulations. 


#### 2_programs.R
Library of functions used across R scripts

#### 3_gen_contact.R
Load contact matrices from /data/contact_msa\<msa code\>_\<source\>, expand types according to /data/msa_type.csv, and define contacts under different policies.

#### 4_grid_search.R
We calibrate three parameters, using grid search

	beta1: initial transmission rates assuming no mask and social distancing measures before lockdown policies. This corresponds to M4;
	beta2: reduced transmission rates assuming adoption of mask and distancing measures. This corresponds to M1;
	I0: initial condition in terms of the fraction of people in each type that are infected but have not yet developed symptoms

#### 5_seir.R
Simulate SEIR models


### data/
This folder includes the following input data

#### Contact matrices
A csv file for each MSA of the contact rates by age and industry, named as "contact_msa<msa code>_<source>".
We have two sources for contact matrices: "replica" or "fred". 
Replica contact matrices is generated based on data from Replica. 
FRED contact matrices are generated from public available synthetic population (https://fred.publichealth.pitt.edu/)

Contact rates mean for a focal individual, the expected number of contact with a target individual at each contact level. 
Each observation is for a focal age (age_i), focal industry (naics_i), target age (age_j), target industry (naics_j), and contact level (contactlvl), with the following variables

	num_people: number of people in each focal age and industry type;
	num_contact_per_person_day: the expected number of contact with people in each target age and industry type;
	contact_time_min_per_person_day: the expected number of contact weighted by duration in minutes, which we use to define contact rates in our main analysis.

Age cohort largely follows CDC age group definition with more cohorts for the elder population. 
Industry cohort is based on NAICS 2 digit code. 
Please refer to age_doc.csv and naics_doc.csv in supplement/
Contact level includes household(hh), school, work, and neighborhood.


#### Deaths and cases in the US
Time series of cumulative deaths and cases by MSA in the US: covid_case_death_jh.csv 
We download deaths and cases in the US from COVID-19 Data Repository by the Center for Systems Science and Engineering (CSSE) at Johns Hopkins University, 
and aggregate county level information to MSA level. 



####  Adjustment for work and health types
msa_type.csv contains the distribution of work and health types, conditional on previously defined age and industry types in each MSA. 

	sick: binary health type indicating whether the individual is high risk (obese or diabetic) based on MEPS (Medical Expenditure Panel Survey);
	sick_w: probability an individual in a MSA X age X industry is high risk;
	wfh: binary type for whether the individual can work from home, computed from O*NET following Dingel and Neiman (2020);
	wfh_w: probability an individual in a MSA X age X industry can work from home;
	shift: a binary split representing alternating schedules for non-essential workers and students;
	shift_w: probability an individual in a MSA X age X industry is in each of the two alternating schedules


### parameter/
This folder includes user input parameters

####  Lockdown and reopening policies and start dates in each MSA: msa_policy_scenarios_dates
Each row corresponds to a MSA. 

	Scenario\<X\> : correspond to the policy code for work, school and neighbor contacts in phase X;
	Date\<X\> : start dates for each phase X;
	T\<X\> : translate the start dates to integer as the number of days from March 5th, 2020;
	Tend: the number of periods the model simulates results for. 


#### Settings for calibrating parameters: gridsearch.csv
Each row corresponds to a MSA. We specify the lower bound (lb), upper bound (ub) and step size (step) for the first round of grid search for the three parameters we calibrate (beta1, beta2, I0). T_range_start and T_range_end specify the start and end dates of the time series we use to measure the goodness of fit. 


#### SEIR model parameters: seir_parameters.csv
Disease specific parameters are taken from related literature. 
We allow these parameters to vary by age and health type.
Please refer to the paper for further descriptions. 

### output/
/csv/: SEIR model simulation results for each compartment by type and period

/figure/: plot time series of selected compartments for internal checking purposes

### temp/
Temporary folders to store intermediate outputs:
1. Expanded contact matrices for each MSA and policy in /contactmatrix/, generated from 3_gen_contact.R
2. Calibrated transmission and initial condition parameter in /calibratedparameter/, generated from 4_grid_search.R 

### supplement/
This folder contains information on our age and industry cohorts


## Instructions for Replication
Specify MSAs (msaList) and reopening policies (reopenPolicy) to run in 1_config.R 

Generate contact matrices with 3_gen_contact.R

Calibrate parameters with 4_grid_search.R. Check reference policies in each phases specified in parameter/msa_policy_scenarios_dates.csv
If grid search hits the boundary of some parameters or is running slow, please revise set up in parameter/gridsearch

Run SEIR to generate simulation results under different reopening policies. To save compartment by type by period outputs, set outputSIR<<-1 in 1_config.R


