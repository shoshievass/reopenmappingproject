##############################################################################
##############################################################################
######                                                                  ######
######           reopen mapping project -- run all                      ######
######                                                                  ######
##############################################################################
##############################################################################



rm(list = ls())

## set up directory
dir<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)
setwd('..')
proj<-getwd()
dataPath <- file.path(proj, "data")
tempPath <- file.path(proj, "temp")
contactMatrixPath <- file.path(tempPath, "contactmatrix")
calibratedParPath <- file.path(tempPath, "calibratedparameter")
codePath <- file.path(proj, "R")
outPath  <- file.path(proj, "output")
parmPath <- file.path(proj, "parameter")

#check if project directory is properly set up
stopifnot(endsWith(proj, "reopenmappingproject"))


## set up global variables and functions
source(file.path(codePath,"2_programs.R"))
source(getCodePath("1_config.R"))


## generate contact matrix
source(getCodePath("3b_gen_contact.R"))

## run grid search 
source(getCodePath("4_grid_search.R"))

## run SIR
source(getCodePath("5_seir.R"))

## run SIR and compare MSAs
source(getCodePath("5B_seir.R"))




