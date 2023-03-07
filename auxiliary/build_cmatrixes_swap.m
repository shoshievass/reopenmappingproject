places = ["5600", "1600", "3760","6920" ]; %histograms
sizes  = ["5600", "1600", "3760","6920" ]; %sizes
contacts = ["5600", "1600", "3760","6920" ]; %contacts


policies=[ "_W4-S3-N3-B2-R2-P2-F2-E2","_W2-S1-N1-B1-R1-P1-F1-E2","_W4-S3-N1-B2-R2-P2-F2-E1", "_W1-S1-N1-B1-R1-P1-F1-E2"];
policies=[ "_W2-S1-N1-B2-R2-P2-F2-E2"];

base_cmat = "/Users/matteo/Desktop/original_cmatrices/C_msa";


 


colnames = strcat("rate", string(1:796));

%% Types to keep 
filename_types =  strcat(base_cmat,"3760", policies(1), ".csv");
C_mat_types = readtable(filename_types);
ego_id = C_mat_types(:,9);
C_mat_types = C_mat_types(:,1:6);

%% Do swaps
for msa = 1:length(places)
    for s = 1: length(sizes)
        for contacts  = 1:length(places)
            for p=1:length(policies)
        
        %%  import size mat
        filename_size = strcat(base_cmat,sizes(s), policies(p), ".csv");
        C_mat_size = readtable(filename_size);

        types_to_keep_size = ismember( C_mat_size(:,1:6), C_mat_types ,'row');
        
        % only keep size and appropriate types
        C_mat_size = C_mat_size(types_to_keep_size,1:9);
        
        current_size= sum(C_mat_size.n);

                
        %%  import main mat
        filename_msa = strcat(base_cmat,places(msa), policies(p), ".csv");
        C_mat_to_modify = readtable(filename_msa);

        types_to_keep_main = ismember(C_mat_to_modify(:,1:6) , C_mat_types  ,'row');

       
        
        % Keep population data from main mat
        C_mat_to_modify = C_mat_to_modify(types_to_keep_main,1:9);
        
        % Standardize histogram
        C_mat_to_modify.n = C_mat_to_modify.n/sum(C_mat_to_modify.n); 
                
        % re-normalize histogram
        C_mat_to_modify.n = C_mat_to_modify.n * current_size;
        
        % check that new histogram is of correct size
        check = abs(sum(C_mat_to_modify.n) - current_size)<1 ;
        assert(check)

        
        %% Import contact matrix
        filename_contact = strcat(base_cmat,places(contacts), policies(p), ".csv");
        C_mat_contacts   = readtable(filename_contact);
        types_to_keep_contacts = ismember(C_mat_contacts(:,1:6) , C_mat_types  ,'row');

        % parse contact rates into main matrix
        C_mat_to_modify(:,10:805) = C_mat_contacts(types_to_keep_contacts,10:805);
        C_mat_to_modify.Properties.VariableNames(10:805) = colnames;
        
        % replace ego ids
        C_mat_to_modify(:,9)= ego_id; 
        
        % write file
        filename_out = strcat(base_cmat,places(msa),"size", sizes(s) , "contacts", places(contacts), policies(p), ".csv");
        writetable(C_mat_to_modify,filename_out)
            
            end
        end
    end
end

