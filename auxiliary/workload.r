rm(list = ls(all.names = TRUE)) #clear all objects includes hidden objects.


#####################################
## Paths
#####################################

dir<-"/accounts/grad/msaccarola/covid_proj/for_slurm"
setwd(dir)
setwd('..')
proj<-"/accounts/grad/msaccarola/covid_proj/for_slurm"




library(deSolve)
library(plyr)
library(dplyr)
library(tidyr)
library(gtools)
library(latex2exp)



#parameters to iterate over



# Define policy and name
policies_1_names    = c("WFH", "60+", "EO", "NP")
policies_1          = c("_W2-S1-N1-B1-R1-P1-F1-E2","_W4-S3-N1-B2-R2-P2-F2-E1", "_W1-S1-N1-B1-R1-P1-F1-E2",  "_W4-S3-N3-B2-R2-P2-F2-E2")

# Size and composition couplets
msa_composition            = c(c("1600","1600"),  c("5600","5600"),  c("6920") ,c("3760"))
msa_size                   = c(c("1600","5600"),  c("5600","1600"),  c("6920") ,c("3760"))
msa_comp_size = cbind(msa_composition, msa_size)

# Contacts
msa_contact         = c("1600","5600", "6920" ,"3760")

# define rest of parameters
beta_0_vector       = c(0.0004, 0.0008, 0.0016)
beta_1_ratio        = seq(0.1,1,.1)
T_chage_phase       = seq(2,50,2)
T_change_beta       = c(10,20)

# Iterate over the 4 msa compositions
tot_iterations =  length(msa_contact) * length(policies_1) *length(T_chage_phase) * 
                      length(beta_0_vector)*length(beta_1_ratio)* (length(T_change_beta)+1) * 
                      length(msa_size)
                  
                    

outcomes_names = c( "msa_composition",
                    "msa_size",
                    "c_matrix",
                    "policy_1_name",
                    "policy_1_code",
                    "beta_0",
                    "beta_1",
                    "T_change_policy",
                    "T_change_beta",
                    "slurm_id")


outcomes <-data.frame(matrix(nrow = tot_iterations,
                             ncol = length(outcomes_names) ))

colnames(outcomes) = outcomes_names

#####################################################
### RUN SEIIIRD SIMULATIONS FOR EACH COMBO, SAVE CSVS
#####################################################
iter = 0
slurm_id=0


for (msa in 1:length(msa_size)) {
  composition = msa_comp_size[msa,1]
  size = msa_comp_size[msa,2]
  
  for (cmat in 1:length(msa_contact)) {
    c_matrix = msa_contact[cmat]
  
  for (pol in 1:length(policies_1)){
    policy_1 = policies_1[pol]
    
    for (change in 1:length(T_chage_phase)){
      T = T_chage_phase[change]
      slurm_id = slurm_id+1  
      
      for (b0 in 1:length(beta_0_vector) ){
        beta_0 = beta_0_vector[b0]
        
        for (b1 in 1:length(beta_1_ratio)){
          beta_1 = beta_1_ratio[b1]*beta_0
          
          for (T_b in 1:length(T_change_beta) ) {
          iter = iter +1 
          outcomes$msa_composition[iter]  = composition
          outcomes$msa_size[iter]         = size
          outcomes$c_matrix[iter]         = c_matrix
          outcomes$policy_1_name[iter]    = policies_1_names[pol]
          outcomes$policy_1_code[iter]    = policies_1[pol]
          outcomes$beta_0[iter]           = beta_0
          outcomes$beta_1[iter]           = beta_1
          outcomes$T_change_policy[iter]  = T
          outcomes$T_change_beta[iter]    = T_change_beta[T_b]
          outcomes$slurm_id[iter]         = slurm_id
          
          }
          
          iter = iter +1 
          outcomes$msa_composition[iter]  = composition
          outcomes$msa_size[iter]         = size
          outcomes$c_matrix[iter]         = c_matrix
          outcomes$policy_1_name[iter]    = policies_1_names[pol]
          outcomes$policy_1_code[iter]    = policies_1[pol]
          outcomes$beta_0[iter]           = beta_0
          outcomes$beta_1[iter]           = beta_1
          outcomes$T_change_policy[iter]  = T
          outcomes$T_change_beta[iter]    = T
          outcomes$slurm_id[iter]         = slurm_id
          
              }
            }
          }
        }
      }
    }

# save final table
write.table( outcomes,file="workload_SLURM", sep=",",col.names=TRUE,row.names=FALSE)