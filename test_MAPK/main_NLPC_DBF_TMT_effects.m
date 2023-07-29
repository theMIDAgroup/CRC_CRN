clc
clear
close all


% This code computes the equilibrium points of CR-CRNs affected by a GoF of k-Ras
% after adding the combination of drugs Dabrafenib and Trametinib (with different
% dosages) to the network

addpath(fullfile('../funcs'))

%% Step 1. Define general parameters
% 1.1. Data
target_folder = '../data/';
file_mim_clean = fullfile(target_folder, 'CRC_CRN_nodrug_complete.mat');

% 1.2. Folders and files
folder_results = './results_paper';
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');

% 1.3. Starting mutation 
mut_prot = 'Ras';
drug1 = 'DBF'; drug2 = 'TMT';
drug = strcat(drug1, '_', drug2);

perc_all = 0;

%% Step 2. Load and store data
load(file_mim_clean, 'CMIM'); 
load(file_ris_phys, 'nlpc_phys')

x_eq_phys = nlpc_phys(1).x; 
rate_constants_phys = CMIM.rates.std_values;

%% 

for ip = 1:numel(perc_all)

    perc = perc_all(ip);

%% Step 4. Simulate mutation and drug one after the other
    % 4.1. Mutation
    x_0_mut = x_eq_phys; rate_constants_mut = rate_constants_phys;

    [x_0_mut, rate_constants_mut] = ...
        f_define_mutated_condition(mut_prot, x_0_mut, ...
        rate_constants_mut, CMIM, perc); 
    
    S_mut = CMIM.matrix.S;
    Nl = CMIM.matrix.Nl;
    rho_mut = Nl * x_0_mut;
    ind_one = CMIM.matrix.ind_one;
    idx_basic_species = find(CMIM.species.std_initial_values>0);
    v = CMIM.matrix.v;
    max_counter = 500;
    
    disp('*****     Mutation    *****')
    disp('Solving with NLPC... ')
    time_init = tic;
    nlpc_mut = f_NLPC_restart(x_0_mut, rate_constants_mut, S_mut, Nl, ...
            rho_mut, idx_basic_species, v, ind_one, max_counter, 0);
    nlpc_mut.elapse_time = toc(time_init);
    clear tim_init
    disp('Done!')
    
    x_eq_mut = nlpc_mut.x;
    
    % 3.2. Add drugs to CRN 

    [CMIM_drug, n_new_species1] = f_add_drug_Raf_from_file(CMIM, drug1);
    [CMIM_drug, n_new_species2] = f_add_drug_Raf_from_file(CMIM_drug, drug2);
    n_new_species = n_new_species1 + n_new_species2;
    
    k1_drug = 0.106 * 1e-3; k2_drug = 0.593 * 1e-4; k3_drug = k1_drug; k4_drug = 1.2296 * 1e-3;
    k5_drug = k1_drug; k6_drug = k4_drug; k7_drug = 0.1 * 1e-1; k8_drug = 0.33 * 1e-2;
    k = [k1_drug k2_drug k3_drug k4_drug k5_drug k6_drug k7_drug k8_drug];
    
    [~, idx_k1] = ismember('cd_1', CMIM_drug.rates.names); [~, idx_k2] = ismember('cd_2', CMIM_drug.rates.names);
    [~, idx_k3] = ismember('cd_3', CMIM_drug.rates.names); [~, idx_k4] = ismember('cd_4', CMIM_drug.rates.names);
    [~, idx_k5] = ismember('cd_5', CMIM_drug.rates.names); [~, idx_k6] = ismember('cd_6', CMIM_drug.rates.names);
    [~, idx_k7] = ismember('cd_7', CMIM_drug.rates.names); [~, idx_k8] = ismember('cd_8', CMIM_drug.rates.names);
    idx_k = [idx_k1 idx_k2 idx_k3 idx_k4 idx_k5 idx_k6 idx_k7 idx_k8];

    %   3.3. Drugs initial values
    [~, idx_d1] = ismember(drug1, CMIM_drug.species.names);
    [~, idx_d2] = ismember(drug2, CMIM_drug.species.names);
    idx_basic_species_drug = [idx_basic_species; idx_d1; idx_d2];
    init_drug_all1 = linspace(0, 100, 11); init_drug_all1 = init_drug_all1(2:end);
    init_drug_all2 = linspace(0, 2000, 11); init_drug_all2 = init_drug_all2(2:end);
    init_drug_all2 = [init_drug_all2, 240];
    
    for id1 = 1:numel(init_drug_all1)
        for id2 = 1:numel(init_drug_all2)
            init_drug1 = init_drug_all1(id1);
            init_drug2 = init_drug_all2(id2);
            fprintf('*******  Analysing perc = %2.2f, init_drug1 = %2.2f, init_drug2 = %2.2f \n ********', ...
                    perc, init_drug1, init_drug2)

            x_0_mut_drug = padarray(x_eq_mut,[n_new_species 0],0,'post');
            x_0_mut_drug(idx_d1) = init_drug1;
            x_0_mut_drug(idx_d2) = init_drug2;

            rate_constants_mut_drug = padarray(rate_constants_mut,[n_new_species 0],0,'post');
            rate_constants_mut_drug(idx_k) = k;

            ind_one = CMIM_drug.matrix.ind_one;
            Nl = CMIM_drug.matrix.Nl;
            rho_mut_drug = Nl * x_0_mut_drug;
            S_mut_drug = CMIM_drug.matrix.S;
            v_drug = CMIM_drug.matrix.v;
            
            disp('*****     Drug    *****')
            disp('Solving through NLPC... ')
            time_init = tic;
            nlpc_mut_drug = f_NLPC_restart(x_0_mut_drug, rate_constants_mut_drug, S_mut_drug, Nl, ...
                    rho_mut_drug, idx_basic_species_drug, v_drug, ind_one, max_counter, 0);
            nlpc_drug_mut.elapse_time = toc(time_init);
            clear tim_init
            disp('Done!')
            
            x_eq_mut_drug = nlpc_mut_drug.x;

            %% Step 5. Simulate combined effect of mutation and drug
            x_0_combo = padarray(x_eq_phys,[n_new_species 0],0,'post');
            x_0_combo(idx_d1) = init_drug1; x_0_combo(idx_d2) = init_drug2;
            rate_constants_combo = rate_constants_phys; rate_constants_combo(idx_k) = k;

            [x_0_combo, rate_constants_combo] = ...
                f_define_mutated_condition(mut_prot, x_0_combo, ...
                rate_constants_combo, CMIM_drug, perc);

            disp('*****     Combined effect of drug and mutation    *****')
            disp('Solving through NLPC... ')
            time_init = tic;
            rho_combo = Nl*x_0_combo;
            nlpc_combo = f_NLPC_restart(x_0_combo, rate_constants_combo, S_mut_drug, Nl, ...
                    rho_combo, idx_basic_species_drug, v_drug, ind_one, max_counter, 0);
            nlpc_combo.elapse_time = toc(time_init);
            clear aux_mut tim_init
            disp('Done!')
            
            x_eq_combo = nlpc_combo.x;

            %% Step 6. Save
            ris_drug.mut_prot = mut_prot;
            ris_drug.k1 = k1_drug;
            ris_drug.k2 = k2_drug;
            ris_drug.k3 = k3_drug;
            ris_drug.k4 = k4_drug;
            ris_drug.k5 = k5_drug;
            ris_drug.k6 = k6_drug;
            ris_drug.k7 = k7_drug;
            ris_drug.k8 = k8_drug;
            ris_drug.init_drug1 = init_drug1;
            ris_drug.init_drug2 = init_drug2;
            ris_drug.max_counter = max_counter;

            ris_drug.rate_constants_mut = rate_constants_mut;
            ris_drug.rate_constants_mut_drug = rate_constants_mut_drug;
            ris_drug.rate_constants_combo = rate_constants_combo;

            ris_drug.x_eq_mut = x_eq_mut;
            ris_drug.x_eq_mut_drug = x_eq_mut_drug;
            ris_drug.x_eq_combo = x_eq_combo;

            aux_save = sprintf('nlpc_%s_on_mut_%s_%2.2f_%2.2f.mat', drug, ...
                mut_prot, init_drug1, init_drug2);    

            save(fullfile(folder_results, 'drugs', aux_save), 'ris_drug')
            
            clearvars -except max_counter target_folder id1 id2 ...
                idx_basic_species_drug folder_results perc mut_prot perc_all ...
                init_drug_all1 init_drug_all2 k k1_drug k2_drug k3_drug k4_drug ...
                k5_drug k6_drug k7_drug k8_drug idx_k idx_k1 idx_k2 idx_k3 idx_k4 ...
                idx_k5 idx_k6 idx_k7 idx_k8 idx_d1 idx_d2 CMIM_drug ris_phys ...
                x_eq_phys rate_constants_phys x_eq_mut rate_constants_mut ...
                cons_law_target drug n_new_species

        end
    end
    clear x_eq_mut rate_constants_mut
end
