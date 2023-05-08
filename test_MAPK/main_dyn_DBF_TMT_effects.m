clc
clear
close all

% This code computes the dynamic of CR-CRNs affected by a GoF of k-Ras
% after adding the combination of drugs Dabrafenib and Trametinib (with
% different dosages) to the network

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
drug1 = 'DBF'; drug2 = 'TMT';
drug = strcat(drug1, '_', drug2);

max_t = 5*10^7;

%% Step 2. Load and store data
load(file_mim, 'CMIM');
load(file_ris_phys, 'nlpc_phys')

idx_basic_species = find(CMIM.species.std_initial_values>0);
x_eq_phys = nlpc_phys(1).x;
rate_constants_phys = CMIM.rates.std_values;

%% Step 3. Add drug and define drug parameters
% 3.1. Add DBF to CRN
[CMIM_drug, n_new_species1] = f_add_drug_Raf_from_file(CMIM, drug1);
[CMIM_drug, n_new_species2] = f_add_drug_Raf_from_file(CMIM_drug, drug2);
n_new_species = n_new_species1 + n_new_species2;

% 3.2. Rate constants
k1_drug = 0.106 * 1e-3; k2_drug = 0.593 * 1e-4; k3_drug = k1_drug; k4_drug = 1.2296 * 1e-3;
k5_drug = k1_drug; k6_drug = k4_drug; k7_drug = 0.1 * 1e-1; k8_drug = 0.33 * 1e-2;
k = [k1_drug k2_drug k3_drug k4_drug k5_drug k6_drug k7_drug k8_drug];
   
[~, idx_k1] = ismember('cd_1', CMIM_drug.rates.names); [~, idx_k2] = ismember('cd_2', CMIM_drug.rates.names);
[~, idx_k3] = ismember('cd_3', CMIM_drug.rates.names); [~, idx_k4] = ismember('cd_4', CMIM_drug.rates.names);
[~, idx_k5] = ismember('cd_5', CMIM_drug.rates.names); [~, idx_k6] = ismember('cd_6', CMIM_drug.rates.names);
[~, idx_k7] = ismember('cd_7', CMIM_drug.rates.names); [~, idx_k8] = ismember('cd_8', CMIM_drug.rates.names);
idx_k = [idx_k1 idx_k2 idx_k3 idx_k4 idx_k5 idx_k6 idx_k7 idx_k8];

%   3.2. Drug initial value
[~, idx_d1] = ismember(drug1, CMIM_drug.species.names);
[~, idx_d2] = ismember(drug2, CMIM_drug.species.names);

idx_basic_species_drug = [idx_basic_species; idx_d1; idx_d2];
init_drug_all1 = [50, 37.5, 25, 12.5];
init_drug_all2 = [200, 150, 100, 50];
 
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

    for id1 = 1:numel(init_drug_all1)

        for id2 = 1:numel(init_drug_all2)

            init_drug1 = init_drug_all1(id1);
            init_drug2 = init_drug_all2(id2);

            fprintf('*******  Analysing perc = %2.2f, init_drug1 = %2.2f, init_drug2 = %2.2f \n ********', ...
                            perc, init_drug1, init_drug2)
            % 4.2. Drug
            x_0_mut_drug = padarray(x_eq_mut, [n_new_species, 0], 0, 'post'); 
            x_0_mut_drug(idx_d1) = init_drug1; x_0_mut_drug(idx_d2) = init_drug2;

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

            x_0_combo(idx_d1) = init_drug1; x_0_combo(idx_d2) = init_drug2;
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
            ris_drug.init_drug1 = init_drug1;
            ris_drug.init_drug2 = init_drug2;
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

            ris_drug.CMIM_drug = CMIM_drug;
            
            aux_save = sprintf('prova_dyn_%s_on_mut_%s_%2.2f_%2.2f.mat', drug, mut_prot, ...
                init_drug1, init_drug2);
            save(fullfile(folder_results, 'drugs', aux_save), 'ris_drug')

            clearvars -except target_folder folder_results max_t perc mut_prot ...
                perc_all idx_k k idx_d CMIM_drug CMIM x_eq_phys rate_constants_phys ...
                rate_constants_mut drug1 drug2 drug time_mut x_t_mut x_eq_mut ...
                cons_law_target n_new_species init_drug_all1 init_drug_all2 id1 id2 ...
                idx_d1 idx_d2

        end
    end

    clear time_mut x_t_mut x_eq_mut rate_constants_mut
end

