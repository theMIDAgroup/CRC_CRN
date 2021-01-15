clc
clear
close all

% NOTE: Run main_script_single_mutation.m first to compute the values of the 
% protein concentrations in the physiological network as well in the network
% affected by a mutation resulting in the gain of function of KRAS.

%% Step 1. Define general parameters
% 1.1. Data
target_folder = './data';
file_mim = fullfile(target_folder, 'CRC_CRN.mat');
% Data are available upon reasonable request to Dr. Sara Sommariva
% (sommariva at dima.unige.it)

% 1.2. Starting mutation 
mut_prot = 'Ras';
target_prot = 'Raf';

% 1.2. Other files and folders
folder_results = fullfile('.', './results');
file_ris_phys = fullfile(folder_results, 'results_physiological.mat');

max_t = 5*10^7;
perc_all = [0]; % Level of the mutation
perc_clt = [1, 0.75, 0.5, 0.25]; % Drug initial concentration (defined as
                  % function of the constant aggregate concentration within
                  % Braf conservation law=

%% Step 2. Load and store data
load(file_mim, 'CMIM');
load(file_ris_phys, 'ris_phys')

x_eq_phys = ris_phys.x_eq;
rate_constants_phys = CMIM.rates.std_values;

cons_laws = CMIM.matrix.Nl;
[~, idx_target] = ismember(target_prot, CMIM.species.names);
cons_law_target = cons_laws(cons_laws(:, idx_target)==1, :);

%% Step 3. Define drug parameters
%   3.1. Rate constants
k1_drug = 0.5;
k2_drug = k1_drug * 10^-2;
[~, idx_k1] = ismember('cd_1', CMIM.rates.names);
[~, idx_k2] = ismember('cd_2', CMIM.rates.names);

%   3.2. Drug initial value
[~, idx_d] = ismember('DRUG', CMIM.species.names);

for ip = 1:numel(perc_all)

    perc = perc_all(ip);

%% Step 4. Simulate mutation and drug one after the other
    % 4.1. Mutation
    x_0_mut = x_eq_phys; rate_constants_mut = rate_constants_phys;
    [x_0_mut, rate_constants_mut] = ...
        f_define_mutated_condition(mut_prot, x_0_mut, ...
        rate_constants_mut, CMIM, perc);

    disp('*****     Mutation    *****')
    disp('Solving ODE system... ')
    [time_mut, x_t_mut] = ode15s(@(t_, x_) f_odefun_MIM(...
        t_, x_, rate_constants_mut, CMIM, 'Sv'), [0 max_t], x_0_mut);
    disp('Done!')

    x_t_mut = x_t_mut';
    x_eq_mut = x_t_mut(:, end);

    tot_conc_clt = cons_law_target * x_eq_mut;
    init_drug_all = perc_clt * tot_conc_clt;

for id = 1:numel(init_drug_all)

    init_drug = init_drug_all(id);
    fprintf('*******  Analysing perc = %2.2f, init_drug = %2.2f \n ********', ...
            perc, init_drug)

    % 4.2. Drug
    x_0_mut_drug = x_eq_mut; 
    x_0_mut_drug(idx_d) = init_drug;

    rate_constants_mut_drug = rate_constants_mut;
    rate_constants_mut_drug(idx_k1) = k1_drug;
    rate_constants_mut_drug(idx_k2) = k2_drug;

    disp('*****     Drug    *****')
    disp('Solving ODE system... ')
    [time_mut_drug, x_t_mut_drug] = ode15s(@(t_, x_) f_odefun_MIM(...
        t_, x_, rate_constants_mut_drug, CMIM, 'Sv'), [0 max_t], x_0_mut_drug);
    disp('Done!')

    x_t_mut_drug = x_t_mut_drug';
    x_eq_mut_drug = x_t_mut_drug(:, end);

    %% Step 5. Simulate combined effect of mutation and drug
    x_0_combo = x_eq_phys; rate_constants_combo = rate_constants_phys;

    x_0_combo(idx_d) = init_drug;
    rate_constants_combo(idx_k1) = k1_drug;
    rate_constants_combo(idx_k2) = k2_drug;

    [x_0_combo, rate_constants_combo] = ...
        f_define_mutated_condition(mut_prot, x_0_combo, ...
        rate_constants_combo, CMIM, perc);

    disp('*****     Combined effect of drug and mutation    *****')
    disp('Solving ODE system... ')
    [time_combo, x_t_combo] = ode15s(@(t_, x_) f_odefun_MIM(...
        t_, x_, rate_constants_combo, CMIM, 'Sv'), [0 max_t], x_0_combo);
    disp('Done!')

    x_t_combo = x_t_combo';
    x_eq_combo = x_t_combo(:, end);

    %% Step 6. Save
    ris_drug.mut_prot = mut_prot;
    ris_drug.k1 = k1_drug;
    ris_drug.k2 = k2_drug;
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

    aux_save = sprintf('results_drug_on_mut_%s_alpha_%2.2f.mat', ...
        mut_prot, perc_clt(id));
    save(fullfile(folder_results, aux_save), 'ris_drug')
    
    clearvars -except target_folder folder_results max_t perc mut_prot ...
        perc_all init_drug_all k1_drug k2_drug idx_k1 idx_k2 idx_d ...
        CMIM ris_phys x_eq_phys rate_constants_phys ...
        time_mut x_t_mut x_eq_mut perc_clt rate_constants_mut ...
        cons_law_target
    
end

    clear time_mut x_t_mut x_eq_mut rate_constants_mut
end














