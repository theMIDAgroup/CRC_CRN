clc
clear
close all

addpath(fullfile('../funcs'))

%% Step 1. Define input data
% 1.1. Files and folders
target_folder = '../data/';
folder_results = './results_paper';
file_mim = fullfile(target_folder, 'CRC_CRN_nodrug.mat');
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');

% 1.2. Load MIM
load(file_mim, 'new_CMIM'); CMIM = new_CMIM; clear new_CMIM

% 1.3. Define which mutation have been implemented
lof_mutations = {'PTEN'};
gof_mutations = {'Ras'};
lof_mutations_type2 = {};%{'TP53'};
% NOTE: Raf <-> Braf, Ras <-> Kras

% 1.4. Define which mutations to simulate
all_proteins = [gof_mutations, lof_mutations, lof_mutations_type2];
perc_all = [0, 0.3, 0.6];

max_t = 2.5*10^7;
rate_constants = CMIM.rates.std_values;

%% Step 2. Run NLPC on the physiological network
script_physiological_NLPC

%% Step 3. Run NLPC on the mutated network
% 3.a. Load results from physiological condition
load(file_ris_phys, 'nlpc_phys')
x_eq = nlpc_phys(1).x;
x_0 = CMIM.species.std_initial_values;

for ij=1:numel(perc_all)
% 3.b. Run each of the required mutation
    for ip = 1:numel(all_proteins)
        
        protein = all_proteins{ip};
        
        switch protein 
             case gof_mutations
        
              perc = perc_all(ij); 
              script_mutation_NLPC;
              
            case [lof_mutations, lof_mutations_type2]
                
              perc = perc_all(ij);
              script_mutation_NLPC;

        end
        
        clear protein
        
    end
end   


