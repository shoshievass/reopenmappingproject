##############################################################################
##############################################################################
######                                                                  ######
######       reopen mapping project -- generate contact matrix          ######
######                                                                  ######
##############################################################################
##############################################################################



#### foreach MSA
for (m in msaList){
  
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(paste("!! Generating contact for MSA", m))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  
  #####################################
  ### load data
  #####################################
  # load contact
  fn <- paste(dataPath, paste("contact_msa", m, "_", getCmatSource(m), ".csv",sep=""), sep="/")
  C <-checkLoad(fn)
  
  # load types
  TYPE <-checkLoad(paste(dataPath, "msa_type.csv", sep="/"))
  TYPE <-TYPE[TYPE$msa==m,!(colnames(TYPE) %in% c("msaname"))]
  
  # types for individual i and j
  TYPE_i<-TYPE
  TYPE_j<-TYPE
  colnames(TYPE_i)<- paste(colnames(TYPE),"_i",sep="")
  colnames(TYPE_j)<- paste(colnames(TYPE),"_j",sep="")
  
  # extend types for individual i and j
  C <- merge(x = C, y = TYPE_i, by = c("age_i", "naics_i"), all = FALSE)
  C <- merge(x = C, y = TYPE_j, by = c("age_j", "naics_j"), all = FALSE)
  
  rm(TYPE, TYPE_i, TYPE_j)
  
  #####################################
  ### policy definition
  #####################################
  
  ## full type vector
  typeVec <-c("age", "naics", "sick", "wfh", "shift")
  
  
  ## adjust number of people and contact by weight for health, WFH, and shift types 
  C$num_people<-C$num_people * (C$sick_w_i * C$shift_w_i * C$wfh_w_i)
  C$num_contact_per_person_day<-C$num_contact_per_person_day * (C$sick_w_j * C$shift_w_j * C$wfh_w_j) 
  
  
  # type space for individual i and j
  # define contact between type i and j, and work status for each type i
  Th<-C[,c(paste(typeVec,"i",sep="_"),paste(typeVec,"j",sep="_"),"essential_i","essential_j","contactlvl")]
  ThW<-Th
  
  # no policy, keep all contacts and everyone worker is employed
  Th$NP <-1
  ThW$NP<-1
  
  # essential only: 
  #   contact: 10% neighbor, only essential work
  #   employ:  only essential and work from home
  Th$EO <- (Th$contactlvl=="hh") + 
              0.1*(Th$contactlvl=="neighbor") + 
              (Th$contactlvl=="work" & Th$essential_i==1 & Th$essential_j==1)
  ThW$EO<- pmax(ThW$essential_i,ThW$wfh_i)
    
  # cautious reopening:
  #   contact: 10% neighbor, normal school/work
  #   employ:  full employment
  Th$CR <- (Th$contactlvl=="hh") + 
              0.1*(Th$contactlvl=="neighbor") + 
              (Th$contactlvl=="work") + (Th$contactlvl=="school")
  ThW$CR<- 1
  
  # isolate 60+:
  #   contact: 10% neighbor, 60+ only essential work
  #   employ:  <60 full employment, 60+ only essential and work from home
  Th$I60<- (Th$contactlvl=="hh") + 
    0.1*(Th$contactlvl=="neighbor") + 
    (Th$contactlvl=="work" &  ((Th$essential_i==1 & Th$essential_j==1) | (Th$age_i<4 & Th$age_j<4)))  +
    (Th$contactlvl=="school")
  ThW$I60<- pmax(ThW$essential_i,ThW$wfh_i,ThW$age_i<4)

  
  # work from home:
  #   contact: 10% neighbor, school, essential and non-work from home at work
  #   employ: full employment
  Th$WFH<- (Th$contactlvl=="hh") + 
    0.1*(Th$contactlvl=="neighbor") + 
    (Th$contactlvl=="work" &  ((Th$essential_i==1 & Th$essential_j==1) | (Th$wfh_i==0 & Th$wfh_j==0)))  +
    (Th$contactlvl=="school")
  ThW$WFH<- 1
  
  # alternate schedule
  #   contact: 10% neighbor, alternate schedule at school and non-essential work
  #   employ: alternate schedule people stay home half the time
  Th$AS <- (Th$contactlvl=="hh") + 
    0.1*(Th$contactlvl=="neighbor") + 
    (Th$contactlvl=="work" &  (Th$essential_i==1 & Th$essential_j==1)) +
    0.5*(Th$contactlvl=="work"   &  (Th$essential_i==0 | Th$essential_j==0) & Th$shift_i==Th$shift_j) +
    0.5*(Th$contactlvl=="school" & Th$shift_i==Th$shift_j)
  ThW$AS<- (ThW$essential_i==1) + (ThW$essential_i==0) * (0.5+0.5*ThW$wfh_i) 
  stopifnot(max(ThW$AS)<=1)
  
  #####################################
  ### contact matrix for each policy
  #####################################
  
  for (p in 1:length(policyList)){
    
    policy <- gsub("_", "", policyList[p])
    Cp <- C
    
    ## adjust number of contact by policy vector 
    Cp$num_contact_per_person_day<-C$num_contact_per_person_day * Th[[policy]]
    
    ## indicate whether these individuals are actively working or work from home
    Cp$active_emp <- ThW[[policy]] * (ThW$naics_i>0)
    
    
    ## aggregate across contact levels
    Cmat<-Cp %>% 
      group_by(age_i, naics_i, sick_i, wfh_i, shift_i, 
               age_j, naics_j, sick_j, wfh_j, shift_j) %>%    # group by type i and j
      summarise(n          = mean(num_people), 
                active_emp = mean(active_emp), 
                contact    = sum(num_contact_per_person_day))  # numebr of people, active work status of i and contact
    
    #unique indicator of types
    Cmat$ego  <- as.integer(interaction(
                  Cmat$age_i, Cmat$naics_i, Cmat$sick_i, Cmat$wfh_i, Cmat$shift_i, drop = TRUE, lex.order = TRUE))
    Cmat$rate <- as.integer(interaction(
                  Cmat$age_j, Cmat$naics_j, Cmat$sick_j, Cmat$wfh_j, Cmat$shift_j, drop = TRUE, lex.order = TRUE))
    
    #reshape long to wide
    Cmat2 <- Cmat[,!(colnames(Cmat) %in% paste(typeVec,"j",sep="_"))] %>% 
            spread(key=rate, value=contact,sep="", fill = 0) %>% 
            rename(age=age_i, naics=naics_i, sick=sick_i, wfh=wfh_i, shift=shift_i)
    
    #check dimensions, 
    #there are 8 additional columns for number of people, active work status, 
    #5 types columns and 1 indicator for all posible types
    stopifnot((dim(Cmat2)[2]-dim(Cmat2)[1])==(length(typeVec)+3))
    
    #export csv
    fn <- paste(dataPath, 
                paste(ctMatData, m, "_", policy, datv, ".csv", sep=""), sep="/")
    write.table(Cmat2, file=fn, sep=",",col.names=TRUE,row.names=FALSE)
    print(paste("export contact matrix:",fn))
  }
  
  rm(C, Th, ThW, Cmat, Cmat2, Cp, policy, typeVec)
}

