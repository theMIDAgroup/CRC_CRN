clc
clear
close all

% This code computes the equilibrium points of CR-CRNs affected by a GoF of k-Ras
% after adding the drug Trametinib (with different dosages) to the network

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
drug = 'TMT';

perc_all = 0; 

%% Step 2. Load and store data
load(file_mim_clean, 'CMIM');
load(file_ris_phys, 'nlpc_phys') 

x_eq_phys = nlpc_phys(1).x;
rate_constants_phys = CMIM.rates.std_values;

%%

for ip = 1:numel(perc_all)

    perc = perc_all(ip);

%% Step 3. Simulate mutation and drug one after the other
    % 3.1. Mutation
    x_0_mut = x_eq_phys; rate_constants_mut = rate_constants_phys;
    
    [x_0_mut, rate_constants_mut] = ...
        f_define_mutated_condition(mut_prot, x_0_mut, rate_constants_mut,...
        CMIM, perc); 
    
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

    % 3.2. Add drug to CRN 
    [CMIM_drug, n_new_species] = f_add_drug_Raf_from_file(CMIM, drug); %target_prot

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
    
    %   3.3. Drug initial value
    [~, idx_d] = ismember('TMT', CMIM_drug.species.names);
    idx_basic_species_drug = [idx_basic_species; idx_d];
    init_drug_all1 = linspace(0, 2000, 11); init_drug_all1 = init_drug_all1(2:end); 
    init_drug_all2 = [1400, 1080, 700, 300];
    init_drug_all = [init_drug_all1 init_drug_all2];

    for id = 1:numel(init_drug_all)

        init_drug = init_drug_all(id);
        fprintf('*******  Analysing perc = %2.2f, init_drug = %2.2f \n ********', ...
                perc, init_drug)
        % 3.3. Drug
        x_0_mut_drug = padarray(x_eq_mut,[n_new_species 0],0,'post');
        x_0_mut_drug(idx_d) = init_drug;

        rate_constants_mut_drug = padarray(rate_constants_mut,[3 0],0,'post');
        rate_constants_mut_drug(idx_k) = k;

        ind_one = CMIM_drug.matrix.ind_one;
        Nl = CMIM_drug.matrix.Nl;
        rho_mut_drug = Nl * x_0_mut_drug;
        S_mut_drug = CMIM_drug.matrix.S;
        v_drug = CMIM_drug.matrix.v;

        disp('*****     Drug    *****')
        disp('Solving with NLPC... ')
        time_init = tic;

        nlpc_mut_drug = f_NLPC_restart(x_0_mut_drug, rate_constants_mut_drug, S_mut_drug, Nl, ...
                rho_mut_drug, idx_basic_species_drug, v_drug, ind_one, max_counter, 0);
        nlpc_mut_drug.elapse_time = toc(time_init);
        clear tim_init
        disp('Done!')

        x_eq_mut_drug = nlpc_mut_drug.x;

        %% Step 5. Simulate combined effect of mutation and drug
        x_0_combo = padarray(x_eq_phys,[n_new_species 0],0,'post'); x_0_combo(idx_d) = init_drug;
        rate_constants_combo = rate_constants_phys;
        rate_constants_combo(idx_k) = k;

        [x_0_combo, rate_constants_combo] = f_define_mutated_condition(mut_prot, ...
            x_0_combo, rate_constants_combo, CMIM_drug, perc);

        disp('*****     Combined effect of drug and mutation    *****')
        disp('Solving with NLPC... ')
        time_init = tic;
        rho_combo = Nl*x_0_combo;
        nlpc_combo = f_NLPC_restart(x_0_combo, rate_constants_combo, S_mut_drug, Nl, ...
                rho_combo, idx_basic_species_drug, v_drug, ind_one, max_counter, 0);
        nlpc_combo.elapse_time = toc(time_init);
        clear tim_init
        disp('Done!')

        x_eq_combo = nlpc_combo.x;

        %% Step 6. Save
        ris_drug.mut_prot = mut_prot;
        ris_drug.k1 = k1_drug; ris_drug.k2 = k2_drug; ris_drug.k3 = k3_drug;
        ris_drug.k4 = k4_drug; ris_drug.k5 = k5_drug;
        ris_drug.k6 = k6_drug; ris_drug.init_drug = init_drug;
        ris_drug.max_counter = max_counter;

        ris_drug.rate_constants_mut = rate_constants_mut;
        ris_drug.rate_constants_mut_drug = rate_constants_mut_drug;
        ris_drug.rate_constants_combo = rate_constants_combo;

        ris_drug.x_eq_mut = x_eq_mut;
        ris_drug.x_eq_mut_drug = x_eq_mut_drug;
        ris_drug.x_eq_combo = x_eq_combo;


        aux_save = sprintf('nlpc_%s_on_mut_%s_%2.2f.mat', drug, mut_prot,...
            init_drug_all(id));    

        save(fullfile(folder_results, 'drugs', aux_save), 'ris_drug')

        clearvars -except max_counter x_eq_mut_drug target_folder id idx_basic_species_drug...
            folder_results  perc mut_prot perc_all init_drug_all k k1_drug k2_drug k3_drug k4_drug k5_drug k6_drug ...
            idx_k idx_k1 idx_k2 idx_k3 idx_k4 idx_k5 idx_k6 idx_d CMIM_drug ris_phys x_eq_phys ...
            rate_constants_phys x_eq_mut perc_clt rate_constants_mut target_prot drug n_new_species

    end

    clear x_eq_mut rate_constants_mut
end
