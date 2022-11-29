SLURM_wrapper <- function( task_id){

    outputSIR<<-1                    #=1 if you want output to be saved
  verTag <<-"SLURMv2"            #tag for output files 
  
  proj <- "/accounts/grad/msaccarola/covid_proj/for_slurmv2"
  setwd(proj)
  
  
  contactMatrixPath <<- file.path(proj, "contactmatrix")
  codePath <<- file.path(proj, "R")
  outPath  <<- file.path(proj, "output")
  parmPath <<- file.path(proj, "parameter")
  
  
  source(file.path(codePath,"functions.R"))
  
  
  dir.create(file.path(outPath), showWarnings = TRUE)
  dir.create(file.path(outPath,"csv"), showWarnings = TRUE)
  
  #####################################
  # Keep useful parm combos and gen matrix of results
  #####################################
  toimport  <<- file.path(proj,"SLURM_workload.csv")
  workload <-checkLoad(toimport)
  
  toimport  <<- file.path(proj,"parameter","demog_accounting.csv")
  demog_accounting <-checkLoad(toimport)
 
  
  
  simulations_to_run = workload %>% filter(slurm_id == task_id)
  

  #uncomment this line for test runs
  #simulations_to_run = simulations_to_run[c(1,2),]
  
  outcomes_names = c( "msa_composition",
                      "msa_size",
                      "c_matrix",
                      "N_population",
                      "policy_1_name",
                      "policy_1_code",
                      "beta_0",
                      "beta_1",
                      "T_change_policy",
                      "T_change_beta",
                      "death_T",
                      "emplost_T",
                      "ICUdays_T",
                      "HOSPdays_T",
                      "death_100",
                      "emplost_100",
                      "ICUdays_100",
                      "HOSPdays_100")
  
  
  outcomes <-data.frame(matrix(nrow = nrow(simulations_to_run),
                               ncol = length(outcomes_names) ))
  
  colnames(outcomes) = outcomes_names
  
  
  # GEN TABLE OF OUTCOMES BY AGEBIN
  colnames_agebin <-c ("age",
      "msa_composition",
      "msa_size",
      "c_matrix",
      "policy_1_name",
      "policy_1_code",
      "beta_0",
      "beta_1",         
      "T_change_policy",
      "T_change_beta",
      "n",
      "Death_T",
      "Death_100",
      "Emplost_T",
      "Emplost_100",
      "HOSPdays_T",
      "HOSPdays_100",
      "ICUdays_T",      
      "ICUdays_100"
      )            
  outcomes_agebin_master <-data.frame(matrix(nrow = 0 ,
                                            ncol = length(colnames_agebin)  ))
  colnames(outcomes_agebin_master) <- colnames_agebin  
  
  
  
  # GEN TABLE OF OUTCOMES BY RACEINC
  colnames_raceinc <-c ("race",
                        "income",
                       "msa_composition",
                       "msa_size",
                       "c_matrix",
                       "policy_1_name",
                       "policy_1_code",
                       "beta_0",
                       "beta_1",         
                       "T_change_policy",
                       "T_change_beta",
                       "n",
                       "Death_T",
                       "Death_100",
                       "Emplost_T",
                       "Emplost_100",
                       "HOSPdays_T",
                       "HOSPdays_100",
                       "ICUdays_T",      
                       "ICUdays_100"
                       )            
  outcomes_raceinc_master <-data.frame(matrix(nrow = 0 ,
                               ncol = length(colnames_raceinc)  ))
  colnames(outcomes_raceinc_master) <- colnames_raceinc  
  
  
  #GEN TABLE OF OUTCOMES BY NAICS
  colnames_naics <-c ("naics",
                        "msa_composition",
                        "msa_size",
                        "c_matrix",
                        "policy_1_name",
                        "policy_1_code",
                        "beta_0",
                        "beta_1",         
                        "T_change_policy",
                        "T_change_beta",
                        "n",
                        "Death_T",
                        "Death_100",
                        "Emplost_T",
                        "Emplost_100",
                        "HOSPdays_T",
                        "HOSPdays_100",
                        "ICUdays_T",      
                        "ICUdays_100"
                        )            
  outcomes_naics_master <-data.frame(matrix(nrow =0,
                               ncol = length(colnames_naics)  ))
  colnames(outcomes_naics_master) <- outcomes_naics_master  
  
  #GEN TABLE OF OUTCOMES BY WFH
  colnames_wfh <-c ("wfh",
                      "msa_composition",
                      "msa_size",
                      "c_matrix",
                      "policy_1_name",
                      "policy_1_code",
                      "beta_0",
                      "beta_1",         
                      "T_change_policy",
                      "T_change_beta",
                      "n",
                      "Death_T",
                      "Death_100",
                      "Emplost_T",
                      "Emplost_100",
                      "HOSPdays_T",
                      "HOSPdays_100",
                      "ICUdays_T",      
                      "ICUdays_100"
                      )            
  outcomes_wfh_master <-data.frame(matrix(nrow = 0,
                                            ncol = length(colnames_wfh)  ))
  colnames(outcomes_wfh_master) <- colnames_wfh
  
  #GEN TABLE OF OUTCOMES BY SICK
  colnames_sick <-c ("sick",
                    "msa_composition",
                    "msa_size",
                    "c_matrix",
                    "policy_1_name",
                    "policy_1_code",
                    "beta_0",
                    "beta_1",         
                    "T_change_policy",
                    "T_change_beta",
                    "n",
                    "Death_T",
                    "Death_100",
                    "Emplost_T",
                    "Emplost_100",
                    "HOSPdays_T",
                    "HOSPdays_100",
                    "ICUdays_T",      
                    "ICUdays_100")
  
  outcomes_sick_master <-data.frame(matrix(nrow = 0,
                                          ncol = length(colnames_sick)  ))
  colnames(outcomes_sick_master) <- colnames_sick
  
                         
  
  #####################################
  # load and define SEIR parameters
  #####################################
  # SEIR model parameters
  seirParm  <<- file.path(parmPath,"seir_parameters.csv")
  
  # contact matrix file names start with
  ctMatData <<- "C_msa"
  
  #compartments
  COMPART   <<-c("S","E","Ia","Ins","Ihc","Rq","Rqd","Rnq","D")
  
  # load user defined SIR parameters
  PAR <-checkLoad(seirParm)
  
  #constant parameters
  betaH <-min(PAR$beta_healthcare)
  gamma <-1/min(PAR$gamma_inv)
  gammaD<-1/min(PAR$gammaD_inv)
  psi   <-min(PAR$psi)
  eta   <-1/min(PAR$eta_inv)
  vparameters0 <<- c(gamma=gamma,betaH=betaH,eta=eta,psi=psi)
  
  #type specific parameters
  EPSILON<<-PAR$epsilon
  TAU    <<-PAR$tau
  
  # rate at which Ia flow into hospital/icu
  # hospital/icu are  (hospital|symptomatic) and  (icu|symptomatic)
  # to compute flow from Ia we need to multiply by the prob of symptom (psi) and the rate (TAU)
  HOSPITAL <<-PAR$hospital * TAU * psi
  ICU      <<-PAR$icu * TAU * psi
  
  # hospitalization/icu duration
  HOSPDAYS<<-min(PAR$hospital_days)
  ICUDAYS <<-min(PAR$icu_days)
  
  
  # input mortality conditional on (hospita|symptomatic) from https://healthcare-in-europe.com/en/news/covid-19-high-mortality-in-hospital-patients.html
  # transform into (death|symptomatic) * rate (gammaD)
  DELTAhc <<-PAR$death_in_hosp * PAR$hospital * gammaD         # death rate 
  GAMMA   <<-gammaD - DELTAhc    # recovery rate to account for variation in death, so total transition rate out of infected is kept at gammaD
  
  ## unique age X sick type: age*10 + sick
  typeAgeSick <<-as.matrix(PAR$age*10+PAR$sick)
  
  #wtd average duration of infected (not in the SIR model, a scaling factor in calibration exercise)
  infectDuration<<-psi * (1/mean(TAU)) + (1-psi) * (1/mean(TAU) + 1/gamma)
  
  ## naics code for healthcare
  heathNAICS<<-62
  
  
  
  
  #  main part: loop over policies and parameters
  policy_0 = "_W4-S3-N3-B2-R2-P2-F2-E2"
  
  
  
  
  
  
  
  #####################################################
  ### RUN SEIIIRD SIMULATIONS FOR EACH COMBO, SAVE CSVS
  #############################
  ########################
  
  for (iter in 1:nrow(simulations_to_run)) {
    histogram    = simulations_to_run$msa_composition[iter]
    size         = simulations_to_run$msa_size[iter]
    c_matrix     = simulations_to_run$c_matrix[iter]
    policy_1     = simulations_to_run$policy_1_code[iter]
    T            = simulations_to_run$T_change_policy[iter]
    T_beta       = simulations_to_run$T_change_beta[iter]
    beta_0       = simulations_to_run$beta_0[iter]
    beta_1       = simulations_to_run$beta_1[iter]
    
    # Case 1 T=T_beta (2 phases)
    if (T==T_beta){
    policy = c(policy_0,policy_1)
    par    = list(0.15, beta_0, beta_1)
    names(par) <- c("I0", "beta_0", "beta_1")
    TTT    <<- c(0,T,100)
    }
    
    # Case 2: T>T_beta: change betas before you change policy
    if (T>T_beta){
      policy = c(policy_0,policy_0,policy_1)
      par    = list(0.15, beta_0,beta_1 ,beta_1)
      names(par) <- c("I0", "beta_0", "beta_mid", "beta_1")
      TTT    <<- c(0,T_beta,T,100)
    }
    
    
    # Case 3: T<T_beta: change betas after you change policy
    if (T<T_beta){
      policy = c(policy_0,policy_1,policy_1)
      par    = list(0.15, beta_0,beta_0 ,beta_1)
      names(par) <- c("I0", "beta_0", "beta_mid", "beta_1")
      TTT    <<- c(0,T,T_beta,100)
    }
    
    # initial condition: number of people in I^A per type
    initNumIperType<<-par$I0
    
    # run main SIR program
    temp_results<-runSir(histogram, size, c_matrix, as.vector(policy), par)
    
    ################################
    # compute relevant aggregates
    ##############################
    outcomes$msa_composition[iter] = histogram
    outcomes$msa_size[iter] = size
    outcomes$c_matrix[iter] = c_matrix
    outcomes$policy_1_name[iter] = simulations_to_run$policy_1_name[iter]
    outcomes$policy_1_code[iter] = policy_1
    outcomes$beta_0[iter] = beta_0
    outcomes$beta_1[iter] = beta_1
    outcomes$T_change_policy[iter] = T
    outcomes$T_change_beta[iter] = T_beta
    
    ## importCSVs
    pcombo<-paste(policy, collapse="" )
    compartList    <-c(COMPART, "active")
    nj<-length(compartList)
    
    for(j in 1:nj){
      #gen file name for compartment
      filename = file.path(outPath, "csv", 
                           paste("sir_", compartList[j], "_" ,"hist",
                                 histogram,"size", size,
                                 'contact',c_matrix, pcombo, 
                                 "_T", paste(TTT[-1], collapse="-"),
                                 "beta_0_", toString(par$beta_0),
                                 "_beta_1_", toString(par$beta_1), "_", 
                                 verTag, ".csv", sep=""))
      
      
      #import compartment
      assign(compartList[j],checkLoad(filename))
      
    }
    
    
    
    
    
    # deaths
    outcomes$death_T[iter] = sum(D[,T+9])
    outcomes$death_100[iter] = sum(D[,100+9])
    
    # employment
    actual_workdays   = active[,1:100+9] * (E[,1:100+9]+S[,1:100+9]+Ia[,1:100+9]+Ins[,1:100+9]+Rnq[,1:100+9]+Rqd[,1:100+9])
    
    tot_workdays_T    = sum(active$T0 * active$n * T)
    true_workdays_T   = sum(actual_workdays[,1:T])
    
    outcomes$emplost_T[iter]        = (tot_workdays_T - true_workdays_T)*10^(-6)
    
    tot_workdays_100     = sum(active$T0 * active$n * 100)
    true_workdays_100    = sum(actual_workdays[,1:100])
    outcomes$emplost_100[iter] = (tot_workdays_100 - true_workdays_100)*10^(-6)
    
    # hospitalization
    ageSickVec         = match(Ia$age*10 + Ia$sick, typeAgeSick)
    epsilonV          = EPSILON[ageSickVec]
    flow_in_Ia        = epsilonV*as.matrix(E[,1:100+9])
    
    # flow_in_Ia = as.matrix(Ia[,1:TTT[3]+9])
    
    flow_in_hosp_days = matrix(nrow=length(D$ego), ncol=100)
    flow_in_icu_days  = flow_in_hosp_days
    
    for (type in 1:length(D$ego)){
      flow_in_hosp_days[type,] = flow_in_Ia[type,] * HOSPITAL[ageSickVec[type]] * HOSPDAYS
      flow_in_icu_days[type,]  = flow_in_Ia[type,] * ICU[ageSickVec[type]] * ICUDAYS
    }
    

    outcomes$ICUdays_T[iter]  =  sum(flow_in_icu_days[,1:T])
    outcomes$HOSPdays_T[iter] =  sum(flow_in_hosp_days[,1:T])
    
    outcomes$ICUdays_100[iter]  =  sum(flow_in_icu_days[,1:100])
    outcomes$HOSPdays_100[iter] =  sum(flow_in_hosp_days[,1:100])
    
    
    
    # replace flow in Ia with Ia itself and uncomment this line to compare against Hanbin's result
    # print(outcomes$HOSPdays_100[iter])
    # print(sum(as.matrix(hospital[,1:TTT[3]+9])) * HOSPDAYS)
    
    
    # population
    outcomes$N_population[iter] = sum(D$n)
    
    
    #########################################
    # Compute relevant aggregates by all bins
    ##########################################
    outcomes_bytypes = types
    
    # Death
    outcomes_bytypes$Death_100 = D[,100+9]
    outcomes_bytypes$Death_T   = D[,T+9]
    
    
    # Employment
    outcomes_bytypes$tot_workdays_T   = active$T0 * active$n * T
    outcomes_bytypes$tot_workdays_100 = active$T0 * active$n * 100
    
    outcomes_bytypes$true_workdays_T     = rowSums(actual_workdays[,1:T])
    outcomes_bytypes$true_workdays_100   = rowSums(actual_workdays[,1:100])
    
    outcomes_bytypes$emplost_T      = (outcomes_bytypes$tot_workdays_T - outcomes_bytypes$true_workdays_T)*10^(-6)
    outcomes_bytypes$emplost_100    = (outcomes_bytypes$tot_workdays_100 - outcomes_bytypes$true_workdays_100)*10^(-6)
    
    
    # Hospitalization
    outcomes_bytypes$ICUdays_100  =  rowSums(flow_in_icu_days[,1:100])
    outcomes_bytypes$HOSPdays_100 =  rowSums(flow_in_hosp_days[,1:100])
    
    outcomes_bytypes$ICUdays_T          =  rowSums(flow_in_icu_days[,1:T])
    outcomes_bytypes$HOSPdays_T         =  rowSums(flow_in_hosp_days[,1:T])
    
    
    outcomes_temp <- outcomes_bytypes %>% 
        group_by(age, naics, sick, wfh) %>% 
        summarise(  Death_T      = sum(Death_T),
                    Death_100    = sum(Death_100),
                    Emplost_T    = sum(emplost_T), 
                    Emplost_100  = sum(emplost_100),
                    HOSPdays_T   = sum(HOSPdays_T) ,
                    HOSPdays_100 = sum(HOSPdays_100),
                    ICUdays_T    = sum(ICUdays_T),
                    ICUdays_100  = sum(ICUdays_100),
                    n            = sum(n)) 
    
    demog_accounting_touse <-  demog_accounting %>% filter( msa == histogram)
    
    merged_outcomes <- merge(outcomes_temp, demog_accounting_touse, by=c("age", "naics", "sick", "wfh"), all=TRUE)
    
    merged_outcomes <- merged_outcomes %>% mutate(   
                                  n            = n            * p_race_inc,
                                  Death_T      = Death_T      * p_race_inc,
                                  Emplost_T    = Emplost_T    * p_race_inc,
                                  Death_100    = Death_100    * p_race_inc,
                                  Emplost_100  = Emplost_100  * p_race_inc,            
                                  ICUdays_T    = ICUdays_T    * p_race_inc,
                                  HOSPdays_T   = HOSPdays_T   * p_race_inc,
                                  ICUdays_100  = ICUdays_100  * p_race_inc,  
                                  HOSPdays_100 = HOSPdays_100 * p_race_inc )
    
    
    ### BY RACE AND INCOME
    outcomes_receincome <- merged_outcomes %>% 
      group_by(race, income) %>% 
      summarise(  Death_T      = sum(Death_T),
                  Death_100    = sum(Death_100),
                  Emplost_T    = sum(Emplost_T), 
                  Emplost_100  = sum(Emplost_100),
                  HOSPdays_T   = sum(HOSPdays_T) ,
                  HOSPdays_100 = sum(HOSPdays_100),
                  ICUdays_T    = sum(ICUdays_T),
                  ICUdays_100  = sum(ICUdays_100),
                  n            = sum(n))
    
    outcomes_receincome$msa_composition = histogram
    outcomes_receincome$msa_size = size
    outcomes_receincome$c_matrix = c_matrix
    outcomes_receincome$policy_1_name = simulations_to_run$policy_1_name[iter]
    outcomes_receincome$policy_1_code = policy_1
    outcomes_receincome$beta_0 = beta_0
    outcomes_receincome$beta_1  = beta_1
    outcomes_receincome$T_change_policy  = T
    outcomes_receincome$T_change_beta  = T_beta
    
    outcomes_raceinc_master <- rbind(outcomes_raceinc_master, outcomes_receincome)
    

    
    ## BY AGEBIN 
    outcomes_byage <- merged_outcomes %>% 
      group_by(age) %>% 
      summarise(  Death_T      = sum(Death_T),
                  Death_100    = sum(Death_100),
                  Emplost_T    = sum(Emplost_T), 
                  Emplost_100  = sum(Emplost_100),
                  HOSPdays_T   = sum(HOSPdays_T) ,
                  HOSPdays_100 = sum(HOSPdays_100),
                  ICUdays_T    = sum(ICUdays_T),
                  ICUdays_100  = sum(ICUdays_100),
                  n            = sum(n))
    
    outcomes_byage$msa_composition = histogram
    outcomes_byage$msa_size = size
    outcomes_byage$c_matrix = c_matrix
    outcomes_byage$policy_1_name = simulations_to_run$policy_1_name[iter]
    outcomes_byage$policy_1_code = policy_1
    outcomes_byage$beta_0 = beta_0
    outcomes_byage$beta_1  = beta_1
    outcomes_byage$T_change_policy  = T
    outcomes_byage$T_change_beta  = T_beta
    
    outcomes_agebin_master <- rbind(outcomes_agebin_master, outcomes_byage)
    
    
    ## BY NAICS
    outcomes_bynaics <- merged_outcomes %>% 
      group_by(naics) %>% 
      summarise(  Death_T      = sum(Death_T),
                  Death_100    = sum(Death_100),
                  Emplost_T    = sum(Emplost_T), 
                  Emplost_100  = sum(Emplost_100),
                  HOSPdays_T   = sum(HOSPdays_T) ,
                  HOSPdays_100 = sum(HOSPdays_100),
                  ICUdays_T    = sum(ICUdays_T),
                  ICUdays_100  = sum(ICUdays_100),
                  n            = sum(n))
    
    outcomes_bynaics$msa_composition = histogram
    outcomes_bynaics$msa_size = size
    outcomes_bynaics$c_matrix = c_matrix
    outcomes_bynaics$policy_1_name = simulations_to_run$policy_1_name[iter]
    outcomes_bynaics$policy_1_code = policy_1
    outcomes_bynaics$beta_0 = beta_0
    outcomes_bynaics$beta_1  = beta_1
    outcomes_bynaics$T_change_policy  = T
    outcomes_bynaics$T_change_beta  = T_beta
    
    outcomes_naics_master <- rbind(outcomes_naics_master,outcomes_bynaics)
    
    
    ## BY WFH
    outcomes_bywfh <- merged_outcomes %>% 
      group_by(wfh) %>% 
      summarise(  n            = sum(n),
                  Death_T      = sum(Death_T),
                  Death_100    = sum(Death_100),
                  Emplost_T    = sum(Emplost_T), 
                  Emplost_100  = sum(Emplost_100),
                  HOSPdays_T   = sum(HOSPdays_T) ,
                  HOSPdays_100 = sum(HOSPdays_100),
                  ICUdays_T    = sum(ICUdays_T),
                  ICUdays_100  = sum(ICUdays_100),
                  )
    
    outcomes_bywfh$msa_composition = histogram
    outcomes_bywfh$msa_size = size
    outcomes_bywfh$c_matrix = c_matrix
    outcomes_bywfh$policy_1_name = simulations_to_run$policy_1_name[iter]
    outcomes_bywfh$policy_1_code = policy_1
    outcomes_bywfh$beta_0 = beta_0
    outcomes_bywfh$beta_1  = beta_1
    outcomes_bywfh$T_change_policy  = T
    outcomes_bywfh$T_change_beta  = T_beta
    
    outcomes_wfh_master <- rbind(outcomes_wfh_master,outcomes_bywfh)
    
    # BYSICK
    outcomes_bysick <- merged_outcomes %>% 
      group_by(sick) %>% 
      summarise(  n            = sum(n),
                  Death_T      = sum(Death_T),
                  Death_100    = sum(Death_100),
                  Emplost_T    = sum(Emplost_T), 
                  Emplost_100  = sum(Emplost_100),
                  HOSPdays_T   = sum(HOSPdays_T) ,
                  HOSPdays_100 = sum(HOSPdays_100),
                  ICUdays_T    = sum(ICUdays_T),
                  ICUdays_100  = sum(ICUdays_100),
      )
    
    outcomes_bysick$msa_composition = histogram
    outcomes_bysick$msa_size = size
    outcomes_bysick$c_matrix = c_matrix
    outcomes_bysick$policy_1_name = simulations_to_run$policy_1_name[iter]
    outcomes_bysick$policy_1_code = policy_1
    outcomes_bysick$beta_0 = beta_0
    outcomes_bysick$beta_1  = beta_1
    outcomes_bysick$T_change_policy  = T
    outcomes_bysick$T_change_beta  = T_beta
    
    outcomes_sick_master <- rbind(outcomes_sick_master,outcomes_bysick)

    
    ################################
    # Finalize
    ##############################
    
    
    # empty csv folder to save space in server
    compartList    <-c(COMPART, "active","hospital", "icu")
    nj<-length(compartList)
    
    
    for(j in 1:nj){
      #gen file name for compartment
      filename = file.path(outPath, "csv",
                           paste("sir_", compartList[j], "_" ,"hist",
                                 histogram,"size", size,
                                 'contact',c_matrix, pcombo,
                                 "_T", paste(TTT[-1], collapse="-"),
                                 "beta_0_", toString(par$beta_0),
                                 "_beta_1_", toString(par$beta_1), "_",
                                 verTag, ".csv", sep=""))


      #import compartment
      file.remove(filename)}

  }
  
  
  
  # save final table of aggregate outcomes
  fn <- file.path(outPath,
                  paste("aggregate_data",
                        verTag, as.character(task_id) ,".csv", sep=""))
  checkWrite(fn,outcomes, "SIR outputs")
  
  # Outcomes disaggregated by type
  fn_naics <- file.path(outPath,
                      paste("bynaics_data",
                            verTag, as.character(task_id) ,".csv", sep=""))
  checkWrite(fn_naics,outcomes_naics_master, "SIR outputs")

  
  fn_raceinc <- file.path(outPath,
                        paste("byraceinc_data",
                              verTag, as.character(task_id) ,".csv", sep=""))
  checkWrite(fn_raceinc,outcomes_raceinc_master, "SIR outputs")
  
  
  fn_wfh <- file.path(outPath,
                          paste("bywfh_data",
                                verTag, as.character(task_id) ,".csv", sep=""))
  checkWrite(fn_wfh,outcomes_wfh_master, "SIR outputs")
  
  
  fn_age <- file.path(outPath,
                  paste("byage_data",
                        verTag, as.character(task_id) ,".csv", sep=""))
  checkWrite(fn_age,outcomes_agebin_master, "SIR outputs")

  fn_sick <- file.path(outPath,
                      paste("bysick_data",
                            verTag, as.character(task_id) ,".csv", sep=""))
  checkWrite(fn_sick,outcomes_sick_master, "SIR outputs")
  }
  
