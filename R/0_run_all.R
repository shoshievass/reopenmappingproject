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
  dataPath <- paste(proj, "data", sep="/")
  codePath <- paste(proj, "R", sep="/")
  outPath  <- paste(proj, "output", sep="/")
  parmPath <- paste(proj, "parameter", sep="/")
} 
setwd(dataPath)


## set up global varibales and functions
source(paste(codePath,"2_programs.R"  ,sep='/'))
source(paste(codePath,"1_config.R"  ,sep='/'))


## generate contact matrix
source(paste(codePath,"5_gen_contact.R",sep='/'))



## run grid search 
source(paste(codePath,"3_grid_search.R",sep='/'))

## run SIR
source(paste(codePath,"4_sir.R",sep='/'))




