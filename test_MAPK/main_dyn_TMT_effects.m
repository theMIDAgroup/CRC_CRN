clc
clear
close all

% This code computes the dynamic of CR-CRNs affected by a GoF of k-Ras
% after adding the drug Dabrafenib (with different dosages) to the network

addpath(fullfile('../funcs'))

%% Step 1. Define general parameters
% 1.1. Data
target_folder = '../data/ci_servono';
file_mim = fullfile(target_folder, 'CRC_CRN_nodrug.mat');

% 1.2. Folders and files
folder_results = './results_paper';
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');

% 1.3. Starting mutation 
mut_prot = 'Ras';
perc_all = 0;
drug = 'TMT';

max_t = 5*10^7;

%% Step 2. Load and store data
load(file_mim, 'CMIM');
load(file_ris_phys, 'nlpc_phys')

x_eq_phys = nlpc_phys(1).x;
rate_constants_phys = CMIM.rates.std_values;

%% Step 3. Add drug and define drug parameters
% 3.1. Add DBF to CRN
[CMIM_drug, n_new_species] = f_add_drug_Raf_from_file(CMIM, drug);

% 3.2. Rate constants
k1_drug = 0.106 * 1e-3; k2_drug = 1.2296 * 1e-3; k3_drug = k1_drug;
k4_drug = k2_drug; k5_drug = k1_drug; k6_drug = 1.184 * 1e-2;
k = [k1_drug k2_drug k3_drug k4_drug k5_drug k6_drug];
    
[~, idx_k1] = ismember('cd_3', CMIM_drug.rates.names);
[~, idx_k2] = ismember('cd_4', CMIM_drug.rates.names);
[~, idx_k3] = ismember('cd_5', CMIM_drug.rates.names);
[~, idx_k4] = ismember('cd_6', CMIM_drug.rates.names);
[~, idx_k5] = ismember('cd_7', CMIM_drug.rates.names);
[~, idx_k6] = ismember('cd_8', CMIM_drug.rates.names);
idx_k = [idx_k1 idx_k2 idx_k3 idx_k4 idx_k5 idx_k6];
    
%   3.2. Drug initial value
[~, idx_d] = ismember(drug, CMIM_drug.species.names);
init_drug_all = [200, 150, 100, 50];

for ip = 1:numel(perc_all)

    perc = perc_all(ip);

%% Step 4. Simulate mutation and drug one after the other
    % 4.1. Mutation
    x_0_mut = x_eq_phys; rate_constants_mut = rate_constants_phys;
    [x_0_mut, rate_constants_mut] = f_define_mutated_condition(mut_prot, x_0_mut,...
        rate_constants_mut, CMIM, perc);

    disp('*****     Mutation    *****')
    disp('Solving ODE system... ')
    [time_mut, x_t_mut] = ode15s(@(t_, x_) f_odefun_MIM(t_, x_, rate_constants_mut,...
        CMIM, 'Sv'), [0 max_t], x_0_mut);
    disp('Done!')

    x_t_mut = x_t_mut';
    x_eq_mut = x_t_mut(:, end);

for id = 1:numel(init_drug_all)

    init_drug = init_drug_all(id);
    fprintf('*******  Analysing perc = %2.2f, init_drug = %2.2f \n ********', ...
            perc, init_drug)

    % 4.2. Drug
    x_0_mut_drug = padarray(x_eq_mut, [n_new_species, 0], 0, 'post'); 
    x_0_mut_drug(idx_d) = init_drug;

    rate_constants_mut_drug = rate_constants_mut;
    rate_constants_mut_drug(idx_k) = k;
    
    disp('*****     Drug    *****')
    disp('Solving ODE system... ')
    [time_mut_drug, x_t_mut_drug] = ode15s(@(t_, x_) f_odefun_MIM(...
        t_, x_, rate_constants_mut_drug, CMIM_drug, 'Sv'), [0 max_t], x_0_mut_drug);
    disp('Done!')

    x_t_mut_drug = x_t_mut_drug';
    x_eq_mut_drug = x_t_mut_drug(:, end);
    
    %% Step 5. Simulate combined effect of mutation and drug
    x_0_combo = padarray(x_eq_phys, [n_new_species, 0], 0, 'post');
    rate_constants_combo = rate_constants_phys;

    x_0_combo(idx_d) = init_drug;
    rate_constants_combo(idx_k) = k;

    [x_0_combo, rate_constants_combo] = ...
        f_define_mutated_condition(mut_prot, x_0_combo, ...
        rate_constants_combo, CMIM_drug, perc);

    disp('*****     Combined effect of drug and mutation    *****')
    disp('Solving ODE system... ')
    [time_combo, x_t_combo] = ode15s(@(t_, x_) f_odefun_MIM(...
        t_, x_, rate_constants_combo, CMIM_drug, 'Sv'), [0 max_t], x_0_combo);
    disp('Done!')

    x_t_combo = x_t_combo';
    x_eq_combo = x_t_combo(:, end);

    %% Step 6. Save
    ris_drug.mut_prot = mut_prot;
    ris_drug.k = k;
    ris_drug.init_drug = init_drug;
    ris_drug.max_t = max_t;
    
    ris_drug.rate_constants_mut = rate_constants_mut;
    ris_drug.rate_constants_mut_drug = rate_constants_mut_drug;
    ris_drug.rate_constants_combo = rate_constants_combo;

    ris_drug.x_eq_mut = x_eq_mut;
    ris_drug.x_eq_mut_drug = x_eq_mut_drug;
    ris_drug.x_eq_combo = x_eq_combo;

    ris_drug.time_mut = time_mut;
    ris_drug.x_t_mut = x_t_mut;

    ris_drug.time_mut_drug = time_mut_drug;
    ris_drug.x_t_mut_drug = x_t_mut_drug;

    ris_drug.time_combo = time_combo;
    ris_drug.x_t_combo = x_t_combo;

    aux_save = sprintf('dyn_%s_on_mut_%s_%2.2f.mat', drug, mut_prot, ...
        init_drug_all(id));
    save(fullfile(folder_results, 'drugs', aux_save), 'ris_drug')
    
    clearvars -except target_folder folder_results max_t perc mut_prot ...
        perc_all init_drug_all k k1_drug k2_drug k3_drug idx_k idx_k1 idx_k2 idx_k3 idx_d ...
        CMIM_drug CMIM ris_phys x_eq_phys rate_constants_phys drug...
        time_mut x_t_mut x_eq_mut perc_clt rate_constants_mut ...
        cons_law_target folder_results_dyn_deg max_t n_new_species
    
end

    clear time_mut x_t_mut x_eq_mut rate_constants_mut
end



