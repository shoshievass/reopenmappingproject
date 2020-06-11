

rm(list = ls())
dir<-Sys.getenv("HOME")

if (dir=="/Users/hanyang") {
  proj <- paste(dir, "Dropbox", "COVID Project", sep="/")
  dataPath <- paste(proj, "stata", "data", sep="/")
  codePath <- paste(proj, "R", sep="/")
  outPath  <- paste(proj, "R", "output", sep="/")
} 
setwd(dataPath)

## set up global varibales and functions
source(paste(codePath,"002_programs.R"  ,sep='/'))
source(paste(codePath,"001_config.R"  ,sep='/'))

### eigenvalue


##
place<-"6920"
contactList<-c("_regular","_socialdistance","_alternateschoolwork","_shutold", "_wfhreopenschool")
contact<-contactList[2]

##
Cmat<-loadData(place, contact)


eig<-largestEigenvalue(Cmat)


#duration
duration<-PSI * (1/TAU) + (1-PSI) * (1/TAU + 1/gamma)
d<-duration[ageSickVec]
# scale matrix
C<-Cmat * (matrix(1,length(d),1) %*% d)

eig2<-largestEigenvalue(C)

eig
eig2
