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

### indicate the fraction that each type is actively working
activeEmp <- function(Th, W, E) {
  
  if (W==1){
    # essential only
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
  
  if (E==1){
    # 60+ only essential or WFH
    activeW <- (Th$age_i< age60) * activeW + 
               (Th$age_i>=age60) * pmax(Th$essential_i,Th$wfh_i)
  }
  
  return(activeW)
}


##########################################################################
######## main script foreach MSA
##########################################################################

for (m in msaList){
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(paste("!! Generating contact for MSA", m))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  
  #####################################
  ### load data
  #####################################
  # load contact
  C <-loadCmat(paste("contact_msa", m, sep=""))
  
  # we might no have duration weight in all sources of contact matrics
  if (("contact_time_min_per_person_day" %in% colnames(C))==FALSE){
    C$contact_time_min_per_person_day<-C$num_contact_per_person_day * 60 * 16
    print("duration is missing from contact matrix. Weight all contacts equally and scale by 16 hours")
  }
  
  # load types
  TYPE <-checkLoad(msaType)
  TYPE <-TYPE[TYPE$msa==m,!(colnames(TYPE) %in% c("msaname"))]
  
  # types for individual i and j
  TYPE_i<-TYPE
  TYPE_j<-TYPE
  colnames(TYPE_i)<- paste(colnames(TYPE),"_i",sep="")
  colnames(TYPE_j)<- paste(colnames(TYPE),"_j",sep="")
  
  # extend types for individual i and j
  C <- merge(x = C, y = TYPE_i, by = c("age_i"), all = FALSE)
  C <- merge(x = C, y = TYPE_j, by = c("age_j"), all = FALSE)
  rm(TYPE, TYPE_i, TYPE_j)
  
  #####################################
  ### policy definition
  #####################################
  
  ## full type vector
  #typeVec <-c("age", "naics", "sick", "wfh", "shift")
  typeVec <-c("age",  "sick", "wfh", "shift")
  ## adjust number of people and contact by weight for health, WFH, and shift types 
  # we use contact weighted by overlap duration (data is in minute and we translate into hours)
  C$num_people<-C$num_people * (C$sick_w_i * C$shift_w_i * C$wfh_w_i)
  C$contact<-(C$contact_time_min_per_person_day / 60 ) * (C$sick_w_j * C$shift_w_j * C$wfh_w_j) 
  
  # type space for individual i and j
  # define contact between type i and j, and work status for each type i
  Th <-C[,c(paste(typeVec,"i",sep="_"), 
            paste(typeVec,"j",sep="_"),
            "essential_i","essential_j","contactlvl",'working_i','working_j')]
  
  ## contact definition for 
  # W(work), S(school), N(neighbor) contacts, whether we quarantine E(elderly), 
  # we do not consider mask here since that does not affect contacts
  contactList <- unique(contactPolicy[,1:4])
  
  # and M(mask), although M no impact on contact
  for (p in 1:dim(contactList)[1]){
    
    #policy tag
    policy <- gsub("_", "", policyTagString(contactList[p,]))
    Cp <- C
    
    ################################################
    ### define contact under policy intervention
    ################################################
 
    #work contact
    contactW <- workContact(Th, contactList[p,1])   
    #school contact
    contactS <- schoolContact(Th, contactList[p,2]) 
    #neighbor contact
    contactN <- neighborContact(Th, contactList[p,3])
    #elderly  quarantine
    contactE <- elderQuarantine(Th, contactList[p,4]) 
    #contact at each level
    Cp$contactPolicy <- C$contact * (
      (Th$contactlvl=="hh") + 
      (Th$contactlvl=="school")   * contactS + 
      (Th$contactlvl=="neighbor") * contactN + 
      (Th$contactlvl=="work") * (contactW * contactE + (1-contactE) * workContact(Th, 1)) ) # for elderly, only essential
    
    #indicate whether these individuals are actively working or work from home
    Cp$activeEmp <- activeEmp(Th, contactList[p,1], contactList[p,4]) * (Th$working_i>0)
    stopifnot(max(Cp$activeEmp)<=1)  
    #####################################
    ### contact matrix for each policy
    #####################################
    
    ## aggregate across contact levels
    Cmat<-Cp %>% 
      group_by(age_i,  sick_i, wfh_i, shift_i,working_i,
               age_j,  sick_j, wfh_j, shift_j,working_j) %>%    # group by type i and j
      summarise(n          = mean(num_people), 
                active_emp = mean(activeEmp), 
                contact    = sum(contactPolicy))  # number of people, active work status of i and contact
    
    #unique indicator of types
    Cmat$ego  <- as.integer(interaction(
      Cmat$age_i, Cmat$sick_i, Cmat$wfh_i, Cmat$shift_i, drop = TRUE, lex.order = TRUE))
    Cmat$rate <- as.integer(interaction(
      Cmat$age_j,  Cmat$sick_j, Cmat$wfh_j, Cmat$shift_j, drop = TRUE, lex.order = TRUE))
    if(length(unique(Cmat$ego)) != length(unique(Cmat$rate))){
      stop("focal and target types in the contact matrix need to be the same!")
    }

    #reshape long to wide
    Cmat2 <- Cmat[,!(colnames(Cmat) %in% paste(typeVec,"j",sep="_"))] %>% 
      spread(key=rate, value=contact,sep="", fill = 0) %>% 
      rename(age=age_i,  sick=sick_i, wfh=wfh_i, shift=shift_i)
    print(colnames(Cmat))
    #check dimensions, 
    #there are 8 additional columns for number of people, active work status, 
    #5 types columns and 1 indicator for all posible types
    if ((dim(Cmat2)[2]-dim(Cmat2)[1])!=(length(typeVec)+2)){
      stop("contact matrix dimensions do not match")
    }
    
    #export contact matrix csv
    checkWrite(file.path(contactMatrixPath, 
                         paste(ctMatData, m, "_", policy, ".csv", sep="")), 
               Cmat2, "contact matrix")

    #collapse to lower dimension contact matrix
    Cmat1 <- Cmat %>% 
      group_by(ego, age_i,  sick_i, wfh_i, shift_i, age_j) %>%    # group by target type 
      summarise(n=mean(n), rate=sum(contact)) %>% # total contact across type j
      group_by(age_i, age_j) %>%                  # group by focal type
      summarise(Value = weighted.mean(rate, n)) %>% #weighted average contact across type i
      rename(FocalAge = age_i, TargetAge = age_j)
    
    Cmat1$PolicyID <- policy
    
    if (p==1){
      CmatAll <- Cmat1
    }else{
      CmatAll <- rbind(CmatAll,Cmat1)
    }
  }
  
  
  ## for each mask policy, duplicate aggregated contact outputs
  mList <- unique(contactPolicy[,5])
  for (i in 1:length(mList)){
    cmat_i <- CmatAll
    cmat_i$PolicyID <- paste(cmat_i$PolicyID, "-M", i, sep="")
    if (i==1){
      CmatAllOut <- cmat_i
    }else{
      CmatAllOut <- rbind(CmatAllOut, cmat_i)
    }
  }
  #export aggregate contact matrix by age
  checkWrite(file.path(outPath, 
                       paste(ctMatData, m, "_allPolicy", ".csv", sep="")), 
             CmatAllOut, "aggregate contact matrix")  
  
  rm(C, Th, Cmat, Cmat2, Cp, policy, typeVec, contactList, mList)
}

