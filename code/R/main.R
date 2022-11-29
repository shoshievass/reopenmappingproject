

rm(list = ls(all.names = TRUE)) #clear all objects includes hidden objects.

proj <- "/accounts/grad/msaccarola/covid_proj/for_slurmv2"
codePath <- file.path(proj, "R")


#####################################
## check packages and load custom functaions functions
#####################################
packages <- c("deSolve","plyr","dplyr","tidyr","gtools","rslurm")
newPackages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)

library(deSolve)
library(plyr)
library(dplyr)
library(tidyr)
library(gtools)
library(rslurm)

source(file.path(codePath,"slurm_func.R"))

to_run = c(1,2)

SLURM_ids = data.frame(task_id=to_run)

sjob <- slurm_apply(SLURM_wrapper,
                    SLURM_ids,
                    jobname = 'simulation_test',
                    nodes = 20, cpus_per_node = 1, submit = TRUE)


