# Readme
This code is intended to launch a large number of simulations at the same time and submit the batched workload to some SLURM cluster. 

## Definition of Simulation
Each simulations is described as a tuple containing the following information:
1. MSA composition (i.e. a distribution of the population over types)
2. MSA size (i.e. total number of inhabitants)
3. Contact Matrix
4. A reopening policy. For the complete and formal definition refer to the old report nr 13 in our shared Dropbox folder
5. A vector $(\beta_0, \beta_1)$ of parameters for the SEIRD model. One for each of the two phases.
6. $T_p$ time to change policy 
7. $T_{\beta}$ time to change beta

## Other parameters
All simulations share the same $I_0 = 0.15$ and simulate the first $T=100$ days of the epidemic. All remaining disease parameters are as specified in the file `parameter/seir_parameters.csv`. We use three different iterations of the file:
- one with COVID parameters (pre loaded)
- one with COVID parameters but we swap the hospitalization and icu probabilities with numbers corresponding to the flu season (`alternative_parameters`)
- same as above but divide death in hospital by 10 (in `alternative_parameters`)

## What simulations do we run?
The file `SLURM_workload.csv` contains all the vectors identifying the simulations we plan to run. A different workload file can easily be created using the `workload.R` script in the `auxiliary` folder.


## Output readme
A batch produces 6 files:
- `aggregate` aggregates all outcomes 
- `byage` disaggregate outcomes by agebins
	- Legend: (1 is age<16, 2 is  16<age<29, 3 is 30<age<39, 4 is 40<age<49 [...] and 8 is age>80 ) 
- `bynaics` disaggregate outcomes by  2-digits naics. 
	- Legend: For a description of the naics codes consult [this](https://www.census.gov/programs-surveys/economic-census/year/2022/guidance/understanding-naics.html) link. Naics 0 indicates unemployment, retirement or more generally individuals not in the workforce.
- `byraceinc` for this exercise we use the race and income composition of the "histogram" msa used. For more info refer to report nr. 14 in our Dropbox folder. The breakdown by race and income is in the variable `p_race_inc` in the file `parameter\demog_accounting.csv`.
	- Income legend: 1 is less than 125% of federal poverty line 2 is between 125 and 400% of federal poverty line and 3 is higher than 400% of federal poverty line 
	- Race legend: 1 is black 2 is non-black hispanic and 3 is all else
- `bywfh` disaggregate outcomes based on work from home status (legend: 0 is does not WFH  and 1 does wfh)
- `bysick` disaggregate outcomes based on whether the possibility pre-existing comorbidities.  (legend: 0 means no pre-existing conditions)

Outcomes:
- Deaths: total deaths in the bin. Not standardized.
- Employment lost: total in the entire bin (i.e. not standardized). Measured in millions of days.
- ICU days: total in the entire population factoring 5 days for each person entering the icu. The number of days for each person flowing into icu can be modified in the `seir_parameters.csv` file.
- Hosp days: total in the entire population factoring 5 days for each person entering hospitalization. The number of days for each person flowing into hospitals can be modified in the `seir_parameters.csv` file.
- N: number of individuals in bin
NB all outcomes--except N--are reported at $T=T_p$ (i.e. time of policy change) and at $T=100$.


## Code overview
The code is composed of 3 main files:
1. `main.R` main wrapper to launch the run in SLURM batches. Main output are the output in the command that launches the slurm smulations (i.e. number of nodes and cpu per node) and the vector of batches id to run
2. `slurm_func.R` main file takes as input a batch id, runs all the simulation, aggregates results and prints the 6 csv output files corresponding to the batch. See output readme for a more complete descriptions of outputs
3. `functions.R` contains all the subfunctions called by `slurm_func.R` 
NB: remember to manually set the project folder in `main` and `slurm_func` before you start the run!!!


## Auxiliary files
-  `build_cmatrixes_swap.m`: using the 16 contact matrices originally produced by Cody in the `original_cmatrices` subfolder, it creates recombines contacts with compositions and size to create the new set of contact matrices. Also it gets rid of extra type so we have 796 types in all msa's.
- `calibrated_values.xlsx` contains vector of calibrated disease parameters and detailed instructions to replicate the calibration using the code from the master branch.
- `workload.R` generates the `SLURM_workload.csv` file. Refer to the paragraphs above for its function.
- `alternative_parameters` contains parameters for the flu and flu with low mortality runs.