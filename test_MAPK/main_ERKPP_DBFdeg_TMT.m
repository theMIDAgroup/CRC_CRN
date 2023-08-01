clc
clear
close all

% This code computes the dynamic of CR-CRNs affected by a GoF of k-Ras
% adding first the drug Dabrafenib (with different dosages) to the network
% and then (after reaching the equilibrium) the drug TMT

addpath(fullfile('../funcs'))

%% Step 1. Define general parameters
% 1.1. Data
target_folder = '../data';
file_mim = fullfile(target_folder, 'CRC_CRN_nodrug_complete.mat');

% 1.2. Folders and files
folder_results = './results_paper';
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');

% 1.3. Mutation + Drugs
mut_prot = 'Ras';
perc = 0;
drug1 = 'DBF'; drug2 = 'TMT';
drug = strcat(drug1, '_', drug2);

max_counter = 500;
proj = 1;
max_t = 5*10^7;

%% Step 2. Load and store data
load(file_mim, 'CMIM');
load(file_ris_phys, 'nlpc_phys')
idx_basic_species = find(CMIM.species.std_initial_values > 0);
x_eq_phys = nlpc_phys(1).x;
rate_constants_phys = CMIM.rates.std_values;

ind_pl = find(strcmp(CMIM.species.names, 'ERKPP'));
ind_ERK = find(strcmp(CMIM.species.names, 'ERK'));
N = CMIM.matrix.Nl;
ind_cl_ERK = find(N(:,ind_ERK) == 1);
ERK_tot = N(ind_cl_ERK,:) * x_eq_phys; 

%% Step 3. Add drug DBF and define drug parameters
% 3.1. Add DBF to CRN
[CMIM_drug, n_new_species1] = f_add_drug_Raf_from_file_deg(CMIM, drug1);
n_new_species = n_new_species1;
k1_drug = 0.106 * 1e-3; k2_drug = 0.593 * 1e-4; k3_drug = 0.579 * 1e-5;
k = [k1_drug k2_drug k3_drug];
[~, idx_k1] = ismember('cd_1', CMIM_drug.rates.names); 
[~, idx_k2] = ismember('cd_2', CMIM_drug.rates.names);
[~, idx_k3] = ismember('cd_3', CMIM_drug.rates.names);

idx_k = [ idx_k1 idx_k2 idx_k3 ];

[~, idx_d1] = ismember(drug1, CMIM_drug.species.names);
init_drug_all1 = 50;

% 3.1. Add TMT to CRN
%[CMIM_drug_1_2, n_new_species2] = f_add_drug_Raf_from_file_after(CMIM_drug, drug2, 1);
[CMIM_drug_1_2, n_new_species2] = f_add_drug_Raf_from_file(CMIM_drug, drug2);
n_new_species_1_2 = n_new_species1 + n_new_species2;

% 3.2. Rate constants
k4_drug = k1_drug; k5_drug = 1.2296 * 1e-3;
k6_drug = k1_drug; k7_drug = k4_drug; 
k8_drug = 0.1 * 1e-1; k9_drug = 0.33 * 1e-2;
k_1_2=[ k4_drug k5_drug k6_drug k7_drug k8_drug k9_drug];

[~, idx_k4] = ismember('cd_4', CMIM_drug_1_2.rates.names); [~, idx_k5] = ismember('cd_5', CMIM_drug_1_2.rates.names); 
[~, idx_k6] = ismember('cd_6', CMIM_drug_1_2.rates.names); [~, idx_k7] = ismember('cd_7', CMIM_drug_1_2.rates.names); 
[~, idx_k8] = ismember('cd_8', CMIM_drug_1_2.rates.names); [~, idx_k9] = ismember('cd_9', CMIM_drug_1_2.rates.names);

idx_k_1_2 = [idx_k4 idx_k5 idx_k6 idx_k7 idx_k8 idx_k9];

%   3.3. Drug initial value
[~, idx_d2] = ismember(drug2, CMIM_drug_1_2.species.names);
idx_basic_species_drug_1_2 = [idx_basic_species; idx_d2];
init_drug_all2 = 240;

%% Step 4. Simulate mutation and drug one after the other
% 4.1. Equilibrium when mutation
x_0_mut = x_eq_phys; rate_constants_mut = rate_constants_phys;
[x_0_mut, rate_constants_mut] = f_define_mutated_condition(mut_prot, x_0_mut,...
    rate_constants_mut, CMIM, perc);

disp('*****     Mutation    *****')
disp('Solving through NLPC... ')
ris_mut = f_NLPC_restart(x_0_mut, rate_constants_mut, CMIM.matrix.S, ...
    CMIM.matrix.Nl,CMIM.matrix.Nl*x_0_mut, idx_basic_species,...
    CMIM.matrix.v, CMIM.matrix.ind_one, max_counter, proj);
disp('Done!')

x_eq_mut=ris_mut.x;

% 4.2 Equilibrium when mutation + DBF
init_drug1 = init_drug_all1;
x_0_mut_drug = [x_eq_mut;zeros(n_new_species,1)];
x_0_mut_drug(idx_d1) = init_drug1; 

rate_constants_mut_drug = [rate_constants_mut; zeros(numel(k),1)];
rate_constants_mut_drug(idx_k) = k;

disp('*****     Drug  DBF  *****')
disp('Solving ODE system... ')
[time_mut_drug, x_t_mut_drug] = ode15s(@(t_, x_) f_odefun_MIM(...
    t_, x_, rate_constants_mut_drug, CMIM_drug, 'Sv'), [0 max_t], x_0_mut_drug);
disp('Done!')

x_t_mut_drug = x_t_mut_drug';

t_min = find(abs((x_t_mut_drug(ind_pl,:)/ERK_tot)-(x_eq_phys(ind_pl)/ERK_tot))<0.01, 1, 'first');


% 4.3 Equilibrium when drug TMT is added
init_drug2 = init_drug_all2;

%fprintf('*******  Analysing perc = %2.2f, init_drug1 = %2.2f, init_drug2 = %2.2f \n ********', ...
%                perc, init_drug1, init_drug2)
rate_constants_mut_drug_1_2 = [rate_constants_mut_drug; zeros(numel(k_1_2),1)];
rate_constants_mut_drug_1_2(idx_k_1_2) = k_1_2;

x_0_mut_drug_1_2 = [x_t_mut_drug(:, t_min);zeros(n_new_species2,1)];
x_0_mut_drug_1_2(idx_d2) = init_drug2; 

disp('*****     Drug TMT    *****')
disp('Solving NLPC... ')
ris_drug_1_2 = f_NLPC_restart(x_0_mut_drug_1_2, rate_constants_mut_drug_1_2, CMIM_drug_1_2.matrix.S, ...
    CMIM_drug_1_2.matrix.Nl,CMIM_drug_1_2.matrix.Nl*x_0_mut_drug_1_2, idx_basic_species_drug_1_2,...
    CMIM_drug_1_2.matrix.v, CMIM_drug_1_2.matrix.ind_one, max_counter, proj);
disp('Done!')

x_eq_mut_drug_1_2 = ris_drug_1_2.x;

disp('Solving ODE system t_0=0... ')
[time_mut_drug_1_2_0, x_t_mut_drug_1_2_0] = ode15s(@(t_, x_) f_odefun_MIM(...
    t_, x_, rate_constants_mut_drug_1_2, CMIM_drug_1_2, 'Sv'), [0 max_t], x_0_mut_drug_1_2);
disp('Done!')

x_t_mut_drug_1_2_0 = x_t_mut_drug_1_2_0';

disp('Solving ODE system t_0=t*... ')
[time_mut_drug_1_2_t_star, x_t_mut_drug_1_2_tstar] = ode15s(@(t_, x_) f_odefun_MIM(...
    t_, x_, rate_constants_mut_drug_1_2, CMIM_drug_1_2, 'Sv'), [time_mut_drug(t_min) max_t], x_0_mut_drug_1_2);
disp('Done!')

x_t_mut_drug_1_2_tstar = x_t_mut_drug_1_2_tstar';

% Simulate mutation and drugs combined (degradation DBF + TMT)
x_0_mut_drug_1_2_combo = [x_eq_mut;zeros(n_new_species_1_2,1)];
x_0_mut_drug_1_2_combo(idx_d1) = init_drug1; 
x_0_mut_drug_1_2_combo(idx_d2) = init_drug2;


disp('Solving ODE system DBF deg + TMT... ')
[time_mut_drug_1_2_combo, x_t_mut_drug_1_2_combo] = ode15s(@(t_, x_) f_odefun_MIM(...
    t_, x_, rate_constants_mut_drug_1_2, CMIM_drug_1_2, 'Sv'), time_mut_drug , x_0_mut_drug_1_2_combo);
disp('Done!')

x_t_mut_drug_1_2_combo = x_t_mut_drug_1_2_combo';

%% Step 5. Save
ris_drug_test.drug1 = drug1;
ris_drug_test.drug2 = drug2;

ris_drug_test.init_drug_1 = init_drug1;
ris_drug_test.init_drug_2 = init_drug2;
ris_drug_test.max_t = max_t;
ris_drug_test.max_counter = max_counter;
ris_drug_test.proj = proj;
ris_drug_test.rate_constants_mut = rate_constants_mut;
ris_drug_test.rate_constants_mut_drug = rate_constants_mut_drug;
ris_drug_test.rate_constants_mut_drug_1_2 = rate_constants_mut_drug_1_2;
   
ris_drug_test.t_min = t_min;
ris_drug_test.time_mut_drug = time_mut_drug;
ris_drug_test.t0_0.time_mut_drug_1_2 = time_mut_drug_1_2_0;
ris_drug_test.t0_tstar.time_mut_drug_1_2 = time_mut_drug_1_2_t_star;
ris_drug_test.time_mut_drug_1_2_combo = time_mut_drug_1_2_combo;

ris_drug_test.x_t_mut_drug = x_t_mut_drug;
ris_drug_test.t0_0.x_t_mut_drug_1_2 = x_t_mut_drug_1_2_0;    
ris_drug_test.t0_tstar.x_t_mut_drug_1_2 = x_t_mut_drug_1_2_tstar;
ris_drug_test.x_e_mut_drug_1_2 = x_eq_mut_drug_1_2;
ris_drug_test.x_0_mut_drug_1_2.time5 = x_0_mut_drug_1_2;
ris_drug_test.x_t_mut_drug_1_2_combo = x_t_mut_drug_1_2_combo;

ris_drug_test.indexERKPP = ind_pl;
ris_drug_test.ERK_tot = ERK_tot;

aux_save = sprintf('test_review_ERKPP_TMT_DBF_deg.mat');
    
save(fullfile(folder_results, 'drugs', aux_save), 'ris_drug_test')
