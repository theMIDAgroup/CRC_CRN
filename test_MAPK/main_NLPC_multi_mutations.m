clc
clear
close all

% NOTE: Run main_script_single_mutation.m first to compute the values of the 
% protein concentrations in the physiological network.

addpath('./funcs')

%% Step 1. Define input data
% 1.1. Files and folders
target_folder = '../data/ci_servono';
folder_results = './results_paper';
file_mim = fullfile(target_folder, 'CRC_CRN_nodrug.mat');
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');
% folder_results_old = './results_NLPC';

% 1.2. Define wich mutations to simulate
all_combo = { {'PI3K', 'Ras'}, {'PTEN', 'Ras'}, {'Raf', 'PTEN'}, {'Raf', 'PI3K'}, {'Ras', 'PI3K'}, ....
              {'Ras', 'PTEN'},{'AKT', 'PTEN'},  {'AKT', 'PI3K'}};
all_perc_m1 = 0;
all_perc_m2 = 0;

max_counter = 500;

%% Step 2. Load data
load(file_mim, 'CMIM');
load(file_ris_phys, 'nlpc_phys')

initial_condition = nlpc_phys(1).x; 
initial_rate_constants = CMIM.rates.std_values;

for ic = 1:numel(all_combo)
    
all_proteins = all_combo{ic};
name_proteins = all_proteins;
% for i=1:length(name_proteins)
%     name_proteins{i} = strrep(name_proteins{i}, 'Ras', 'KRAS');
% end

for ip = 1:numel(all_perc_m1)
    
perc_m1 = all_perc_m1(ip);
perc_m2 = all_perc_m2(ip);
        
%% Step 3 Implement the two mutations together
[x_0_mut_comb, rate_constants_mut_comb] = ...
    f_define_mutated_condition(all_proteins{1}, initial_condition, ...
    initial_rate_constants, CMIM, perc_m1);
[x_0_mut_comb, rate_constants_mut_comb] = ...
    f_define_mutated_condition(all_proteins{2}, x_0_mut_comb, ...
    rate_constants_mut_comb, CMIM, perc_m2);
rho_mut = CMIM.matrix.Nl * x_0_mut_comb;
disp('*****     Combined mutation   *****')
disp('Solving through NLPC... ')
ind_one = CMIM.matrix.ind_one;
v = CMIM.matrix.v;
idx_basic_species = find(CMIM.species.std_initial_values>0);

time_init = tic;
aux_mut = f_NLPC_restart(nlpc_phys(1).x0, rate_constants_mut_comb, CMIM.matrix.S,...
    CMIM.matrix.Nl, rho_mut, idx_basic_species, v, ind_one, max_counter, 0);
aux_mut.elapse_time = toc(time_init);
x_eq_mut_comb = aux_mut.x;
nlpc_mut = aux_mut;
clear aux_mut tim_init

%% Step 5. Save
nlpc_combo.proteins = all_proteins;
nlpc_combo.x_eq_mut_comb = x_eq_mut_comb;
nlpc_combo.max_counter = max_counter;

aux_save = sprintf('nlpc_mut_%s_%s.mat', ...
    name_proteins{1}, name_proteins{2});
save(fullfile(folder_results, 'mutations', aux_save), 'nlpc_combo')

clear nlpc_combo

end
end

