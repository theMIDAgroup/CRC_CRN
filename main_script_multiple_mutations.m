clc
clear
close all

% NOTE: Run main_script_single_mutation.m first to compute the values of the 
% protein concentrations in the physiological network.

%% Step 1. Define input data
% 1.1. Data 
target_folder = './data';
file_mim = fullfile(target_folder, 'CRC_CRN.mat');

% 1.2. Others results
folder_results = fullfile('.', './results');
file_ris_phys = fullfile(folder_results, 'results_physiological.mat');

% 1.3. Add necessary folders to path
addpath('./funcs')

% 1.4. Define wich mutations to simulate
all_combo = { {'Raf', 'PTEN'}, {'Raf', 'PI3K'}, {'Ras', 'PI3K'}, ....
              {'Ras', 'PTEN'},{'AKT', 'PTEN'},  {'AKT', 'PI3K'}};
all_perc_m1 = [0];
all_perc_m2 = [0];

max_t = 5*10^7;

%% Step 2. Load data
load(file_mim, 'CMIM');
load(file_ris_phys, 'ris_phys')

initial_condition = ris_phys.x_eq; % Starting from the equilibrium 
                                   % in the physiological cell
initial_rate_constants = CMIM.rates.std_values;

for ic = 1:numel(all_combo)
    
all_proteins = all_combo{ic};

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

disp('*****     Combined mutation   *****')
disp('Solving ODE system... ')
[time_mut_comb, x_t_mut_comb] = ode15s(@(t_, x_) f_odefun_MIM(...
    t_, x_, rate_constants_mut_comb, CMIM, 'Sv'), [0 max_t], x_0_mut_comb);
disp('Done!')

x_t_mut_comb = x_t_mut_comb';
x_eq_mut_comb = x_t_mut_comb(:, end);

%% Step 5. Save
ris_combo.proteins = all_proteins;
ris_combo.x_eq_mut_comb = x_eq_mut_comb;
ris_combo.max_t = max_t;

aux_save = sprintf('results_multi_mutations_%s_%s_pm1_%1.2f_pm2_%1.2f.mat', ...
    all_proteins{1}, all_proteins{2}, perc_m1, perc_m2);
save(fullfile(folder_results, aux_save), 'ris_combo')

clear ris_combo

end
end

