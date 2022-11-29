places = ["5600", "1600", "3760","6920" ]; %histograms
sizes  = ["5600", "1600", "3760","6920" ]; %sizes


policies=[ "_W4-S3-N3-B2-R2-P2-F2-E2","_W2-S1-N1-B1-R1-P1-F1-E2","_W4-S3-N1-B2-R2-P2-F2-E1", "_W1-S1-N1-B1-R1-P1-F1-E2"];
policies=[ "_W4-S3-N3-B2-R2-P2-F2-E2"];

base_cmat = "/Users/matteo/Desktop/Tebaldi Project/matrices/C_msa";

colnames = strcat("rate", string(1:796));

for msa = 1:length(places)
    for s = 1: length(sizes)
        for contacts  = 1:length(places)
            for p=1:length(policies)
        
        %%  import size mat
        filename_size = strcat(base_cmat,sizes(s),"contacts",sizes(s), policies(p), ".csv");
        C_mat_size = readtable(filename_size);
        C_mat_size = C_mat_size(:,1:9);
        
        current_size= sum(C_mat_size.n);

                
        %%  import main mat
        filename_msa = strcat(base_cmat,places(msa),"contacts",places(msa), policies(p), ".csv");
        C_mat_to_modify = readtable(filename_msa);
        
        
        
        % Keep population data from main mat
        C_mat_to_modify = C_mat_to_modify(:,1:9);
        
        % Standardize histogram
        C_mat_to_modify.n = C_mat_to_modify.n/sum(C_mat_to_modify.n); 
                
        % re-normalize histogram
        C_mat_to_modify.n = C_mat_to_modify.n * current_size;
        
        % check that new histogram is of correct size
        check = abs(sum(C_mat_to_modify.n) - current_size)<1 ;
        assert(check)

        
        %% Import contact matrix
        filename_contact = strcat(base_cmat,places(contacts),"contacts",places(contacts), policies(p), ".csv");
        C_mat_contacts   = readtable(filename_contact);
        
        % parse contact rates into main matrix
        C_mat_to_modify(:,10:805) = C_mat_contacts(:,10:805);
        C_mat_to_modify.Properties.VariableNames(10:805) = colnames;
        
        filename_out = strcat(base_cmat,places(msa),"size", sizes(s) , "contacts", places(contacts), policies(p), ".csv");
        writetable(C_mat_to_modify,filename_out)
            
            end
        end
    end
end

