clc
clear
close all

addpath(fullfile('../funcs'))

%% Step 1. Define input data
% 1.1. Files and folders
target_folder = '../data/ci_servono';
folder_results = './results_paper';
file_mim = fullfile(target_folder, 'CRC_CRN_nodrug.mat');
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');

% 1.2. Load MIM
load(file_mim, 'CMIM'); 

% 1.3. Define which mutation have been implemented
lof_mutations = {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT', 'ARF'};
gof_mutations = {'Raf', 'Ras', 'PI3K', 'BetaCatenin'};
lof_mutations = {'PTEN'};
gof_mutations = {'Ras'};
lof_mutations_type2 = {};%{'TP53'};
% NOTE: Raf <-> Braf, Ras <-> Kras

% 1.4. Define which mutations to simulate
all_proteins = [gof_mutations, lof_mutations, lof_mutations_type2];
perc_all = [0, 0.3, 0.6];

% 1.5. Define which part of code to run
do_physiological = 0;
do_mutation = 1;

max_t = 2.5*10^7;
rate_constants = CMIM.rates.std_values;

%% Step 2. If required simulate MIM dynamic in physiological condition 
if do_physiological

    script_physiological_NLPC;
    
end

%% Step 3. If required simulate MIM dynamic in mutated condition
if do_mutation
    
% 3.a. Load results from physiological condition
load(file_ris_phys, 'nlpc_phys')
x_eq = nlpc_phys(1).x;
x_0 = nlpc_phys(1).x0;

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
end

