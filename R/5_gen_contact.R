##############################################################################
##############################################################################
######                                                                  ######
######       reopen mapping project -- generate contact matrix          ######
######                                                                  ######
##############################################################################
##############################################################################



#### foreach MSA
for (m in 1:length(msaList)){
  
  
  #####################################
  ### load data
  #####################################
  #load contact
  C <-checkLoad(paste("contact_msa", msaList[m], datv, ".csv",sep=""))
  
  # load types
  TYPE <-checkLoad("msa_type.csv")
  TYPE<-TYPE[TYPE$msa==msaList[m],!(colnames(TYPE) %in% c("msaname"))]
  
  # types for individual i and j
  TYPE_i<-TYPE
  TYPE_j<-TYPE
  colnames(TYPE_i)<- paste(colnames(TYPE),"_i",sep="")
  colnames(TYPE_j)<- paste(colnames(TYPE),"_j",sep="")
  
  # extend types for individual i and j
  C <- merge(x = C, y = TYPE_i, by = c("age_i", "naics_i"), all = FALSE)
  C <- merge(x = C, y = TYPE_j, by = c("age_j", "naics_j"), all = FALSE)
  
  
  #####################################
  ### policy definition
  #####################################
  ##
  
  # type space for individual i and j
  Th<-C[,c("age_i","naics_i","sick_i","shift_i","wfh_i","essential_i",
           "age_j","naics_j","sick_j","shift_j","wfh_j","essential_j","contactlvl")]
  
  
  # normal
  Th$NP<-1
  
  # essential only
  Th$EO<- (Th$contactlvl=="hh") + 
              0.1*(Th$contactlvl=="neighbor") + 
              (Th$contactlvl=="work" & Th$essential_i==1 & Th$essential_j==1)
  
  # cautious reopening
  Th$CR<- (Th$contactlvl=="hh") + 
              0.1*(Th$contactlvl=="neighbor") + 
              (Th$contactlvl=="work") + (Th$contactlvl=="school")
  
  # isolate 60+
  Th$I60<- (Th$contactlvl=="hh") + 
    0.1*(Th$contactlvl=="neighbor") + 
    (Th$contactlvl=="work" &  ((Th$essential_i==1 & Th$essential_j==1) | (Th$age_i<4 & Th$age_j<4)))  +
    (Th$contactlvl=="school")
  
  # work from home
  Th$WFH<- (Th$contactlvl=="hh") + 
    0.1*(Th$contactlvl=="neighbor") + 
    (Th$contactlvl=="work" &  ((Th$essential_i==1 & Th$essential_j==1) | (Th$wfh_i==0 & Th$wfh_j==0)))  +
    (Th$contactlvl=="school")
  
  # alternate schedule
  Th$AS<- (Th$contactlvl=="hh") + 
    0.1*(Th$contactlvl=="neighbor") + 
    (Th$contactlvl=="work" &  (Th$essential_i==1 & Th$essential_j==1)) +
    0.5*(Th$contactlvl=="work"   &  (Th$essential_i==0 | Th$essential_j==0) & Th$shift_i==Th$shift_j) +
    0.5*(Th$contactlvl=="school" & Th$shift_i==Th$shift_j)
  
  
  
  #####################################
  ### contact matrix for each policy
  #####################################
  
  for (p in 1:length(policyList)){
    
    policy <- gsub("_", "", policyList[p])
    Cp <- C
    
    # weight
    # note we currently have two types shift (for alternating schedule) and wfh. 
    # They are the same but with different weight: shift_w and wfh_w
    # So we need to define weights conditional on policy
    if (policy=="WFH"){
      w_i<-C$sick_w_i * C$wfh_w_i
      w_j<-C$sick_w_j * C$wfh_w_j
    }else{
      w_i<-C$sick_w_i * C$shift_w_i
      w_j<-C$sick_w_j * C$shift_w_j    
    }
    
    ## adjust number of people and contact by weight and policy vector 
    Cp$num_people<-C$num_people * w_i
    Cp$num_contact_per_person_day<-C$num_contact_per_person_day * w_j * Th[[policy]]
    
    
    ## aggregate across contact levels
    Cmat<-Cp %>% 
      group_by(age_i, naics_i, sick_i, shift_i, age_j, naics_j, sick_j, shift_j) %>%    # multiple group columns
      summarise(n = mean(num_people), contact = sum(num_contact_per_person_day))  # multiple summary columns
    
    #unique indicator of types
    Cmat$ego  <- as.integer(interaction(Cmat$age_i, Cmat$naics_i, Cmat$sick_i, Cmat$shift_i, drop = TRUE,lex.order = TRUE))
    Cmat$rate <- as.integer(interaction(Cmat$age_j, Cmat$naics_j, Cmat$sick_j, Cmat$shift_j, drop = TRUE,lex.order = TRUE))
    
    #reshape long to wide
    Cmat2 <- Cmat[,!(colnames(Cmat) %in% c("age_j","naics_j","sick_j","shift_j"))] %>% 
            spread(key=rate,value=contact,sep="", fill = 0) %>% 
            rename(age=age_i, naics=naics_i, sick=sick_i, shift=shift_i)
    
    #check dimensions, there are 6 additional columns for number of people, 4 types columns and 1 indicator for all posible types
    stopifnot((dim(Cmat2)[2]-dim(Cmat2)[1])==6)
    
    #export csv
    fn <- paste(dataPath, 
                paste("C_msa", msaList[m], "_" , policy, CmatVer, ".csv", sep=""), sep="/")
    write.table(Cmat2, file=fn, sep=",",col.names=TRUE,row.names=FALSE)
    print(paste("export contact matrix :",fn))
  
  }
}