##############################################################################
##############################################################################
######                                                                  ######
######       reopen mapping project -- generate contact matrix          ######
######                                                                  ######
##############################################################################
##############################################################################



### contact at work  
workContact <- function(Th, poli) {
  
  #essential worker  
  essentialWorker<- (Th$essential_i==1 & Th$essential_j==1)
  
  if (poli==1){
    # essential only
    contactW <-  essentialWorker
  }else if (poli==2){
    # WFH: essential only + cannot work from home
    contactW <- essentialWorker | (Th$wfh_i==0 & Th$wfh_j==0)
  }else if (poli==3){
    # AS
    contactW <- essentialWorker + 0.5*((essentialWorker==0) & (Th$shift_i==Th$shift_j))
  }else if (poli==4){
    contactW<-1
  }else{
    stop("Only support W 1-4")
  }
  return(contactW)
}

### contact at school  
schoolContact <- function(Th, poli) {
  if (poli==1){
    # no school
    contactS <-  0
  }else if (poli==2){
    # AS
    contactS <- 0.5*(Th$shift_i==Th$shift_j)
  }else if (poli==3){
    contactS<-1
  }else{
    stop("Only support S 1-3")
  }
  return(contactS)  
}  

### contact at neighborhood
neighborContact <- function(Th, poli){
  if (poli==1){
    contactN<-1/3
  }else if (poli==2){
    contactN<-2/3
  }else if (poli==3){
    contactN<-1
  }else{
    stop("Only support N 1-3")
  }
  return(contactN)
}

### whether we quarantee elderly at work
elderQuarantine <- function(Th, poli){
  if (poli==1){
    # keep work contact only if both <60 
    contactE<- (Th$age_i<age60 & Th$age_j<age60)
  }else if (poli==2){
    contactE<-1    
  }else{
    stop("Only support E 1-2")
  }
  return(contactE)
}  


### POI contacts
POIcontact <- function(Th, poli){
  if (poli==1){
    contactPOI<-0
  }else if (poli==2){
    contactPOI<-1    
  }else{
    stop("Only support 1-2")
  }
  return(contactPOI)
}  



### indicate the fraction that each type is actively working
genActiveEmp <- function(Th, W, E, POI) {
  
  if (W==1){
    # essential only, non-essential workers can work from home
    activeW <-pmax(Th$essential_i,Th$wfh_i)
  }else if (W==3){
    # AS
    activeW <-Th$essential_i + (Th$essential_i==0) * (0.5+0.5*Th$wfh_i) 
  }else if (W==4|W==2){
    # WFH or regular, not impact on employment
    activeW <-1
  }else{
    stop("Only support W 1-4")
  }
  
  # work at POI
  if (any(is.na(POI))==FALSE){
    activeW <- activeW * (Th$work_poi_i==0)
    for (i in 1:length(POI)){
      activeW <- activeW + (Th$work_poi_i==i) * (POI[i]==2)
    }
  }
  
  if (E==1){
    # 60+ only essential or WFH
    # POI workers are not essential and cannot WFH so no conditions here
    activeW <- (Th$age_i< age60) * activeW + 
      (Th$age_i>=age60) * pmax(Th$essential_i,Th$wfh_i)
  }
  
  return(activeW)
}


### get variable for type j from type i
typeJfromI <- function(Cp) {

  Cp$idx <- seq.int(nrow(Cp))
  
  # extract status for active employment
  Cp_1 <- Cp %>% 
    group_by(age_i, naics_i, work_poi_i, sick_i, wfh_i, shift_i) %>% 
    summarise( activeEmp_i = mean(activeEmp))

  # switch i and j updated offsite and update status for j
  colnames(Cp_1) <- gsub("_i", "_j", colnames(Cp_1))
  
  #update onsite/offsite status of worker i
  Cp2 <- merge(x = Cp, y = Cp_1, by = paste(typeVec, "j", sep="_"), all = TRUE)
  
  # sort and re-arrange
  Cp2 <-Cp2[order(Cp2$idx),]
  Cp2 <-Cp2[,!(colnames(Cp2) %in% "idx")]
  
  stopifnot(nrow(Cp2)==nrow(Cp)) 
  
  return(Cp2)
}

# conditional on working, whether this individual is working onsite or offsite
OnOrOffWork <- function(Cp, Th, W) {
  
  # no longer employed
  unemploy_i <- (Cp$activeEmp  ==0) * Cp$naics_i>0
  unemploy_j <- (Cp$activeEmp_j==0) * Cp$naics_j>0
  
  if (W==3){
    #AS: AS workers are not at risk, only POI workers are
    Cp$atRiskW   <- Th$atRiskW * (Th$work_poi_i>0)
    
    # essential and POI workers always onsite, those participating in AS onsite half the time
    Cp$offsite_i <- (Cp$activeEmp>0)   * (Th$essential_i==0 & Th$work_poi_i==0) * 0.5 + unemploy_i
    Cp$offsite_j <- (Cp$activeEmp_j>0) * (Th$essential_j==0 & Th$work_poi_j==0) * 0.5 + unemploy_j
    
  }else{
    Cp$atRiskW   <- Th$atRiskW
    
    Cp$offsite_i <- (Cp$activeEmp>0)   * (Th$wfh_i==1) + unemploy_i
    Cp$offsite_j <- (Cp$activeEmp_j>0) * (Th$wfh_j==1) + unemploy_j     
  }
  
  #at risk contact (workers for other industry or POI)
  Cp$atRiskContact <- Cp$contactPolicy *  ((Th$work_poi_i!=Th$work_poi_j) | (Th$naics_i!=Th$naics_j))
  
  return(Cp)
}





##########################################################################
######## main script foreach MSA
##########################################################################


msaList<<-c("5600","1600","6920","3760")
# wfh-close POI, {EO, WFH, AS} X open POI
pList = c(610,1149,1150,1151)
offsiteIterMat <- matrix(0, nrow=50, ncol=length(msaList)*length(pList))
i_col=1
for (m in msaList){
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(paste("!! Generating contact for MSA", m))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  
  #####################################
  ### load data
  #####################################
  # load contact
  C <-loadCmat(paste(ctMatRaw, m, sep=""))
  
  # we might no have duration weight in all sources of contact matrics
  if (("contact_time_min_per_person_day" %in% colnames(C))==FALSE){
    C$contact_time_min_per_person_day<-C$num_contact_per_person_day * 60 * 16
    print("duration is missing from contact matrix. Weight all contacts equally and scale by 16 hours")
  }
  
  # do we have POI types?
  hasPOI <- ("work_poi_i" %in% colnames(C))
  
  # load types
  TYPE <-checkLoad(msaType)
  TYPE <-TYPE[TYPE$msa==m,!(colnames(TYPE) %in% c("msaname"))]
  
  # types for individual i and j
  TYPE_i<-TYPE
  TYPE_j<-TYPE
  colnames(TYPE_i) <- paste(colnames(TYPE),"_i",sep="")
  colnames(TYPE_j) <- paste(colnames(TYPE),"_j",sep="")
  
  # extend types for individual i and j
  if (hasPOI){
    C <- merge(x = C, y = TYPE_i, by = c("age_i", "naics_i", "work_poi_i"), all = FALSE)
    C <- merge(x = C, y = TYPE_j, by = c("age_j", "naics_j", "work_poi_j"), all = FALSE)
  }else{
    C <- merge(x = C, y = TYPE_i, by = c("age_i", "naics_i"), all = FALSE)
    C <- merge(x = C, y = TYPE_j, by = c("age_j", "naics_j"), all = FALSE)
  }
  
  rm(TYPE, TYPE_i, TYPE_j)
  
  
  #####################################
  ### policy definition
  #####################################
  
  ## full type vector
  typeVec <-c("age", "naics", "work_poi", "sick", "wfh", "shift")
  
  ## adjust number of people and contact by weight for health, WFH, and shift types 
  # we use contact weighted by overlap duration (data is in minute and we translate into hours)
  C$num_people<-C$num_people * (C$sick_w_i * C$shift_w_i * C$wfh_w_i)
  C$contact<-(C$contact_time_min_per_person_day / 60 ) * (C$sick_w_j * C$shift_w_j * C$wfh_w_j) 
  C$contact_unwtd<-(C$num_contact_per_person_day ) * (C$sick_w_j * C$shift_w_j * C$wfh_w_j) 
  
  # type space for individual i and j
  # define contact between type i and j, and work status for each type i
  Th <-C[,c(paste(typeVec,"i",sep="_"), 
            paste(typeVec,"j",sep="_"),
            "essential_i","essential_j","contactlvl")]
  
  # at risk workers: non-essential workers who cannot work from home
  atRiskW <- (C$contactlvl=="work" & C$naics_i>0 & C$wfh_i==0 & C$essential_i==0)
  # and POI workers
  if (hasPOI) {
    for (i in 1:length(poiLevel)){
      poi_i <- paste("work",poiLevel[i],sep="-")
      atRiskW <- atRiskW + 
        (C$contactlvl==poi_i & C$naics_i>0 & C$wfh_i==0 & C$essential_i==0)
    }
  }
  Th$atRiskW <- atRiskW
  
  ## contact definition for 
  # W(work), S(school), N(neighbor), B(bar/restuarant), R(retail), P(personal services), F(entertainment)
  # and whether we quarantine E(elderly), 
  # we do not consider mask here since that does not affect contacts
  if (hasPOI){
    contactList <- unique(contactPolicyPOI[,1:8])
  }else{
    contactList <- unique(contactPolicy[,1:4])
  }
  nPoli <- dim(contactList)[2]
  
  #foreach policy under consideration
  for (p in pList){
    #policy tag
    policy <- gsub("_", "", policyTagString(contactList[p,]))
    Cp <- C
    print(paste("contact ", p, "/",dim(contactList)[1], " : ", policy, sep=""))
    
    ################################################
    ### define contact under policy intervention
    ################################################
    
    #work contact
    contactW <- workContact(Th, contactList[p,1])   
    #school contact
    contactS <- schoolContact(Th, contactList[p,2]) 
    #neighbor contact
    contactN <- neighborContact(Th, contactList[p,3])
    #elderly  quarantine (=1 if not quarantined)
    contactE <- elderQuarantine(Th, contactList[p,nPoli]) 
    
    
    #activate contact at each level
    Cp$contactActive <-  
      (Th$contactlvl=="hh") + 
      (Th$contactlvl=="school")   * contactS + 
      (Th$contactlvl=="neighbor") * contactN + 
      (Th$contactlvl=="work") * (contactW * contactE + (1-contactE) * workContact(Th, 1)) #if isolate elderly, allow essential works
    
    #POI contacts, isolate ederly overwrites POI. Other work policies (WFH/AS) no impact on POIs
    # treat those who visit POI the same as those who work at POI, for simple on/off policies
    if (hasPOI) {
      for (i in 1:length(poiLevel)){
        Cp$contactActive <- Cp$contactActive + 
          (Th$contactlvl %in% paste(c("work-", ""), poiLevel[i],sep="")) * POIcontact(Th, contactList[p,3+i]) * contactE
      }
      poiContact <- unlist(contactList[p,3+c(1:length(poiLevel))])
    }else{
      poiContact <- NA
    }
    
    # included contact X contact rates
    Cp$contactPolicy <- Cp$contactActive * C$contact
    Cp$contactPolicy_unwtd <- Cp$contactActive * C$contact_unwtd
    
    #indicate whether these individuals are actively working or work from home
    Cp$activeEmp <- genActiveEmp(Th, 
                              contactList[p,1], contactList[p,nPoli], poiContact) * (Th$naics_i>0)
    stopifnot(max(Cp$activeEmp)<=1)  
    
    #indicate work status for individual j
    Cp <- typeJfromI(Cp)
    stopifnot(max(Cp$activeEmp_j)<=1)


    #####################################
    ### TBD
    #####################################
    
    # conditional on working, whether working offsite
    Cp <- OnOrOffWork(Cp, Th, contactList[p,1])
    
    # only focus on those employed but at risk
    CpFp <- Cp[(Cp$atRiskW * Cp$activeEmp)>0,]
    CpFp$idx <- seq.int(nrow(CpFp))
    
    offsite_old <- CpFp$offsite_i
    
    # fix point
    offsite_iter <-c()
    i <-1
    err<-999
    while (err>1e-5 & i<50){
      
      # what fraction of all contacts at work for individual i is lost because some individual j is not comming on site
      Cp_1 =  CpFp %>% 
        group_by(age_i, naics_i, work_poi_i, sick_i, wfh_i, shift_i) %>% 
        summarise(offsetContact = mean(offsite_j * atRiskContact), 
                  allContact    = mean(atRiskContact),
                  popOffsite    = mean(num_people*offsite_i),
                  pop           = mean(num_people)) %>%
        mutate(offsite_i_prime = offsetContact/allContact,
               offsite_frac = popOffsite/pop) %>%
        select(-(offsetContact:popOffsite))
      
      #show wtd avg offside fraction
      frac_off <- sum(Cp_1$offsite_frac*Cp_1$pop)/sum(Cp_1$pop)
      offsite_iter<-c(offsite_iter, frac_off)
      print(paste("iteration", i, ": off-site fraction", format(frac_off,digits=4), "... error =",format(err,digits=4)))
      
      # sort and re-arrange
      Cp_1 <-Cp_1[, !(colnames(Cp_1) %in% c("offsite_frac", "pop"))]      
      
      #update offsite status of worker i
      CpFp <- merge(x = CpFp, y = Cp_1, by = paste(typeVec, "i", sep="_"), all = TRUE)
      CpFp$offsite_i <- pmax(CpFp$offsite_i,  CpFp$offsite_i_prime, na.rm=TRUE)
      if (min(CpFp$offsite_i)<0 | max(CpFp$offsite_i)>1){
        stop("i's offset indicator out of [0,1] bounds")
      }
      
      # switch i and j updated offsite and update offiste status of worker j
      colnames(Cp_1) <- gsub("_i", "_j", colnames(Cp_1))
      CpFp <- merge(x = CpFp, y = Cp_1, by = paste(typeVec, "j", sep="_"), all = TRUE)
      CpFp$offsite_j <- pmax(CpFp$offsite_j,  CpFp$offsite_j_prime, na.rm=TRUE)
      if (min(CpFp$offsite_j)<0 | max(CpFp$offsite_j)>1){
        stop("j's offset indicator out of [0,1] bounds")
      }
      
      # sort and re-arrange
      CpFp <-CpFp[order(CpFp$idx), 
                  !(colnames(CpFp) %in% c("offsite_i_prime", "offsite_j_prime"))]
      
      offsite_new <- CpFp$offsite_i
      
      err <- max(abs(offsite_new-offsite_old))
      offsite_old <- offsite_new
      i <- i+1
    }

    
    #store onsite fraction 
    offsiteIterMat[1:length(offsite_iter),i_col] <- offsite_iter
    i_col=i_col+1 
  }
}  


offsiteIterMatCut<-offsiteIterMat[1:5,]

for (i in 2:nrow(offsiteIterMatCut)){
  offsiteIterMatCut[i,]<-pmax(offsiteIterMatCut[i,],offsiteIterMatCut[i-1,])
}



idxm<-0

for (m in msaList){
  fn <- file.path(outPath, "figure", paste("frac_non_wfh_still_emp_", m, ".pdf", sep=""))  
  if (fn!="") pdf(fn)
  plot(1:nrow(offsiteIterMatCut), 1:nrow(offsiteIterMatCut), ylim=c(0,25), col="white",
       type="l", ylab="Job lost among at risk workers (%)", xlab="iteration")
  
  
  colorList <- c(t_col("red",perc = 0),t_col("blue", perc = 0),t_col("green",perc = 0), t_col("orange", perc = 0))
  
  for (i in 1:4){
    lines(1:nrow(offsiteIterMatCut), 100*offsiteIterMatCut[,idxm+i], type="l", col=colorList[i])  
  }
  legend("topleft",
         legend=c("WFH-Close POI", "EO-Open POI",
                  "WFH-Open POI", "AS-Open POI"),
         col=colorList, horiz=F,lwd=1.3,bty="n",cex=1)
  if (fn!="") dev.off()
  
  idxm<-idxm+4
  
}


