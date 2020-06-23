##############################################################################
##############################################################################
######                                                                  ######
######           reopen mapping project -- run all                      ######
######                                                                  ######
##############################################################################
##############################################################################



rm(list = ls())

## set up directory
dir<-Sys.getenv("HOME")

if (dir=="/Users/hanyang") {
  proj <- paste(dir, "Documents","GitHub","reopenmappingproject", sep="/")
  stopifnot(dir.exists(proj))
  dataPath <- paste(proj, "data", sep="/")
  codePath <- paste(proj, "R",    sep="/")
  outPath  <- paste(proj, "output", sep="/")
  parmPath <- paste(proj, "parameter", sep="/")
} 
#check if project directory is properly set up
stopifnot(endsWith(proj, "reopenmappingproject"))
setwd(proj)


## set up global varibales and functions
source(paste(codePath,"2_programs.R",sep='/'))
source(getCodePath("1_config.R"))


## generate contact matrix
source(getCodePath("3_gen_contact.R"))

## run grid search 
source(getCodePath("4_grid_search_v2.R"))

## run SIR
source(getCodePath("5_sir.R"))




