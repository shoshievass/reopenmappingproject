##############################################################################
##############################################################################
######                                                                  ######
######                 COVID Project SIR                                ######
######                                                                  ######
##############################################################################
##############################################################################



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


#####################################
# export detailed SIR outputs
#####################################

# foreach place/scenarios
for (p in placeList){
  for (c in contactList){
    Ttt<<-setTiming(c)
    
    place  <-p
    contact<-c
    
    ## load data for dimensions
    CData <<- read.csv(paste("C2_msa", place, contact,".csv",sep=""),header=TRUE)
    
    ## load SIR simulation result
    load(paste(outPath, "model", paste("ode_", place, contact, rver, ".rda", sep=""), sep="/"))
    
    ### export 
    exportSIRbyNAICS(sim1)
  }
}

#####################################
# compute eigen values
#####################################

len_c<-length(contactList)
len_p<-length(placeList)

EIGEN<-matrix(0,len_c,len_p)
EIGENbeta<-matrix(0,len_c,len_p)
AvgDegree<-matrix(0,len_c,len_p)

# foreach place/scenarios
for (p in 1:len_p){
  par<-calibratedPar(placeList[p])
  for (c in 1:len_c){
    
    # unweighted and weighted eigen values
    Cmat<-loadData(placeList[p], contactList[c])
    EIGEN[c,p]<-largestEigenvalue(Cmat, N)
    
    EIGENbeta[c,p]<-EIGEN[c,p]*par$beta
    
    w<-N/sum(pop)
    AvgDegree[c,p]<-sum(w * rowSums(Cmat))
  }
}


out<-cbind(expand.grid(contactList,placeList),as.vector(EIGEN),as.vector(EIGENbeta), as.vector(AvgDegree))

#export csv
fn <- paste(outPath, "csv", 
            paste("max_eigen_values.csv", sep=""), sep="/")
write.table(out, file=fn,sep=",",col.names=FALSE,row.names=FALSE)
