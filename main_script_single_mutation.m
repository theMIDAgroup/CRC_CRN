clc
clear
close all

%% Step 1. Define input data
% 1.1. Data 
target_folder = './data';
file_mim = fullfile(target_folder, 'CRC_CRN.mat');

% 1.2. Other files and folders
folder_results = fullfile('.', 'results');
file_ris_phys = fullfile(folder_results, 'results_physiological.mat');

% 1.3. Add necessary folders to path
addpath(fullfile('./funcs'))

% 1.4. Load 
load(file_mim, 'CMIM');

% 1.5. Define which mutation have been implemented
lof_mutations = {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT', 'ARF'};
gof_mutations = {'Raf', 'Ras', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
% NOTE: Raf <-> Braf, Ras <-> Kras

% 1.6. Define which mutations to simulate
all_proteins = [gof_mutations, lof_mutations, lof_mutations_type2];

% 1.7. Define which part of code to run
do_physiological = 1;
do_mutation = 1;
do_drugs = 0;

max_t = 2.5*10^7;
rate_constants = CMIM.rates.std_values;

%% Step 2. If required simulate MIM dynamic in physiological condition 
if do_physiological

    script_physiological;
    
end

%% Step 3. If required simulate MIM dynamic in mutated condition
if do_mutation
    
% 3.a. Load results from physiological condition
load(file_ris_phys, 'ris_phys')
x_eq = ris_phys.x_eq;
x_0 = ris_phys.x_0;
    
% 3.b. Run each of the required mutation
    for ip = 1:numel(all_proteins)
        
        protein = all_proteins{ip};
        
        switch protein 
             case gof_mutations
        
              perc = 0; 
              script_mutation;
              
            case [lof_mutations, lof_mutations_type2]
                
              perc = 0;
              script_mutation;

        end
        
        clear protein
        
    end
    
end

