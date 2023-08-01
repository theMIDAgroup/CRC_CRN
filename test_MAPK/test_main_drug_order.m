clc
clear
close all

% This code computes the dynamic of CR-CRNs affected by a GoF of k-Ras
% adding first the drug Dabrafenib (with different dosages) to the network
% and secondly (after reaching the equilibrium) the drug TMT

addpath(fullfile('../funcs'))

%% Step 1. Define general parameters
% 1.1. Mutation + Drugs
mut_prot = 'Ras';
perc = 0;
drug1 = 'DBF'; drug2 = 'TMT';
drug = strcat(drug1, '_', drug2);

% 1.2. Data
target_folder = '../data';
file_mim = fullfile(target_folder, 'CRC_CRN_nodrug_complete.mat');

% 1.3. Folders and files
folder_results = './results_paper';
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');

folder_results_drug = './results_paper/drugs';
folder_ris_drug = fullfile(folder_results_drug, sprintf('nlpc_DBF_TMT_on_mut_%s_40.00_240.00', mut_prot));

max_counter = 500;
proj = 1;
max_t = 5*10^7;
n_dyn = 5;

%% Step 2. Load and store data
load(file_mim, 'CMIM'); 
load(file_ris_phys, 'nlpc_phys')
idx_basic_species = find(CMIM.species.std_initial_values>0);
x_eq_phys = nlpc_phys(1).x;
rate_constants_phys = CMIM.rates.std_values;

%% Step 3. Add drug DBF and define drug parameters
% 3.1. Add DBF to CRN
[CMIM_drug, n_new_species1] = f_add_drug_Raf_from_file(CMIM, drug1);
n_new_species = n_new_species1;

k1_drug = 0.106 * 1e-3; k2_drug = 0.593 * 1e-4;
k = [k1_drug k2_drug];
[~, idx_k1] = ismember('cd_1', CMIM_drug.rates.names); [~, idx_k2] = ismember('cd_2', CMIM_drug.rates.names);
idx_k = [ idx_k1 idx_k2 ];

[~, idx_d1] = ismember(drug1, CMIM_drug.species.names);
idx_basic_species_drug = [idx_basic_species; idx_d1];
init_drug_all1 = 40;

% 3.2. Add TMT to CRN
[CMIM_drug_1_2, n_new_species2] = f_add_drug_Raf_from_file(CMIM_drug, drug2);
n_new_species_1_2 = n_new_species1 + n_new_species2;

k3_drug = k1_drug; k4_drug = 1.2296 * 1e-3;
k5_drug = k1_drug; k6_drug = k4_drug; k7_drug = 0.1 * 1e-1; k8_drug = 0.33 * 1e-2;
k_1_2=[ k3_drug k4_drug k5_drug k6_drug k7_drug k8_drug];

[~, idx_k3] = ismember('cd_3', CMIM_drug_1_2.rates.names); [~, idx_k4] = ismember('cd_4', CMIM_drug_1_2.rates.names);
[~, idx_k5] = ismember('cd_5', CMIM_drug_1_2.rates.names); [~, idx_k6] = ismember('cd_6', CMIM_drug_1_2.rates.names);
[~, idx_k7] = ismember('cd_7', CMIM_drug_1_2.rates.names); [~, idx_k8] = ismember('cd_8', CMIM_drug_1_2.rates.names);
idx_k_1_2 = [ idx_k3 idx_k4 idx_k5 idx_k6 idx_k7 idx_k8];

[~, idx_d2] = ismember(drug2, CMIM_drug_1_2.species.names);
idx_basic_species_drug_1_2 = [idx_basic_species_drug; idx_d2];
init_drug_all2 = 240;

%% Step 4. Simulate mutation and drug one after the other
% 4.1. Equilibrium - Mutation
x_0_mut = x_eq_phys; rate_constants_mut = rate_constants_phys;
[x_0_mut, rate_constants_mut] = f_define_mutated_condition(mut_prot, x_0_mut,...
    rate_constants_mut, CMIM, perc);

disp('*****     Mutation    *****')
disp('Solving through NLPC... ')
ris_mut=f_NLPC_restart(x_0_mut, rate_constants_mut, CMIM.matrix.S, ...
    CMIM.matrix.Nl,CMIM.matrix.Nl*x_0_mut, idx_basic_species,...
    CMIM.matrix.v, CMIM.matrix.ind_one, max_counter, proj);
disp('Done!')

x_eq_mut=ris_mut.x;

% 4.2. Equilibrium - Mutation + DBF
init_drug1 = init_drug_all1;
x_0_mut_drug = [x_eq_mut;zeros(n_new_species,1)];
x_0_mut_drug(idx_d1) = init_drug1; 

rate_constants_mut_drug = [rate_constants_mut; zeros(numel(k),1)];
rate_constants_mut_drug(idx_k) = k;

disp('*****     Drug  DBF  *****')

disp('Solving through NLPC... ')
ris_drug_1=f_NLPC_restart(x_0_mut_drug, rate_constants_mut_drug, CMIM_drug.matrix.S, ...
    CMIM_drug.matrix.Nl,CMIM_drug.matrix.Nl*x_0_mut_drug, idx_basic_species_drug,...
    CMIM_drug.matrix.v, CMIM_drug.matrix.ind_one, max_counter, proj);
disp('Done!')
x_eq_mut_drug(:,1)=ris_drug_1.x;

disp('Solving ODE system... ')
[time_mut_drug, x_t_mut_drug] = ode15s(@(t_, x_) f_odefun_MIM(...
    t_, x_, rate_constants_mut_drug, CMIM_drug, 'Sv'), [0 max_t], x_0_mut_drug);
disp('Done!')

x_t_mut_drug=x_t_mut_drug';

% 4.3 Adding drug TMT
init_drug2 = init_drug_all2;

%fprintf('*******  Analysing perc = %2.2f, init_drug1 = %2.2f, init_drug2 = %2.2f \n ********', ...
%                perc, init_drug1, init_drug2)

rate_constants_mut_drug_1_2 = [rate_constants_mut_drug; zeros(numel(k_1_2),1)];
rate_constants_mut_drug_1_2(idx_k_1_2) = k_1_2;

times=sort(randi([1, numel(time_mut_drug)],[5,1]));
x_eq_drug_1_2=zeros(numel(x_eq_mut_drug)+n_new_species2,n_dyn+1);

for i = 1:n_dyn+1
if i <= n_dyn
    x_0_mut_drug_1_2 = [x_t_mut_drug(:, times(i)); zeros(n_new_species2, 1)];
else
    x_0_mut_drug_1_2 = [x_eq_mut_drug; zeros(n_new_species2, 1)];
end
x_0_mut_drug_1_2(idx_d2) = init_drug2; 

disp('*****     Drug TMT    *****')
disp('Solving through NLPC... ')
ris_drug_1_2=f_NLPC_restart(x_0_mut_drug_1_2, rate_constants_mut_drug_1_2, CMIM_drug_1_2.matrix.S, ...
    CMIM_drug_1_2.matrix.Nl, CMIM_drug_1_2.matrix.Nl*x_0_mut_drug_1_2, idx_basic_species_drug_1_2,...
    CMIM_drug_1_2.matrix.v, CMIM_drug_1_2.matrix.ind_one, max_counter, proj);
disp('Done!')

x_eq_drug_1_2(:,i)=ris_drug_1_2.x;

end

%% Step 6. Save
ris_drug_order.drug1 = drug1;
ris_drug_order.drug2 = drug2;
    
ris_drug_order.init_drug_1 = init_drug1;
ris_drug_order.init_drug_2 = init_drug2;
ris_drug_order.max_t = max_t;
ris_drug_order.max_counter = max_counter;
ris_drug_order.proj = proj;
ris_drug_order.rate_constants_mut = rate_constants_mut;
ris_drug_order.rate_constants_mut_drug = rate_constants_mut_drug;
ris_drug_order.rate_constants_mut_drug_1_2 = rate_constants_mut_drug_1_2;
    
ris_drug_order.x_0_mut_drug_1_2.time1 = x_t_mut_drug(:,1);
ris_drug_order.x_0_mut_drug_1_2.time2 = x_t_mut_drug(:,2);
ris_drug_order.x_0_mut_drug_1_2.time3 = x_t_mut_drug(:,3);
ris_drug_order.x_0_mut_drug_1_2.time4 = x_t_mut_drug(:,4);
ris_drug_order.x_0_mut_drug_1_2.time5 = x_t_mut_drug(:,5);

ris_drug_order.x_eq_mut = x_eq_mut;
ris_drug_order.x_eq_mut_drug = x_eq_mut_drug;
ris_drug_order.x_eq_mut_drug_1_2.nlpc = x_eq_drug_1_2(:,6);
ris_drug_order.x_eq_mut_drug_1_2.time1 = x_eq_drug_1_2(:,1);
ris_drug_order.x_eq_mut_drug_1_2.time2 = x_eq_drug_1_2(:,2);
ris_drug_order.x_eq_mut_drug_1_2.time3 = x_eq_drug_1_2(:,3);
ris_drug_order.x_eq_mut_drug_1_2.time4 = x_eq_drug_1_2(:,4);
ris_drug_order.x_eq_mut_drug_1_2.time5 = x_eq_drug_1_2(:,5);
  
ris_drug_order.time_mut_drug = time_mut_drug;
ris_drug_order.x_t_mut_drug = x_t_mut_drug;

ris_drug_order.chosen_times = times;
    
aux_save = sprintf('test_review_order_%.2f_%.2f.mat', init_drug1, init_drug2);
   
save(fullfile(folder_results, 'drugs', aux_save), 'ris_drug_order')




