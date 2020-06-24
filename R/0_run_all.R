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
dataPath <- paste(proj, "data", sep="/")
codePath <- paste(proj, "R",    sep="/")
outPath  <- paste(proj, "output", sep="/")
parmPath <- paste(proj, "parameter", sep="/")

#check if project directory is properly set up
stopifnot(endsWith(proj, "reopenmappingproject"))


## set up global variables and functions
source(paste(codePath,"2_programs.R",sep='/'))
source(getCodePath("1_config.R"))


## generate contact matrix
source(getCodePath("3_gen_contact.R"))

## run grid search 
source(getCodePath("4_grid_search.R"))

## run SIR
source(getCodePath("5_sir.R"))
