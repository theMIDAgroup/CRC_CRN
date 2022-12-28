clear
close all
clc

% Here we use the NLPC method for computing equilibrium points on the 
% physiological and mutated stoichiometric surfaces, using the novel 
% non linear projector. 
% Starting points: x0_all selected with 'main_extract_x0.m'.
% Function used for the simulation: 'f_NLPC_restart' with proj = 0.
% PHYSIOLOGICAL CASE: results in results/'nlpc_phys_testproj.mat'
% MUTATED CASES: results in results/'nlpc_mut_protein.mat', where
% 'protein' represents the mutation considered in each situation.

addpath(fullfile('..', 'funcs'))

do_phys = 1;
do_mutation = 1; perc = 0;

%% Step1. Define path and load
target = fullfile('..', 'data');

path_mim = fullfile(target, 'CRC_CRN_nodrug.mat');
folder_st_points = './results/starting_points';
folder_results = './results/';
file_x0_phys = fullfile(folder_st_points, 'x0_phys.mat');
aux_file_x0_mut = 'x0_%s.mat';

load(path_mim, 'new_CMIM'); CRN = new_CMIM;

%% Step 2. Define general parameters of the network
v = CRN.matrix.v;
Nl = CRN.matrix.Nl;
idx_basic_species = find(CRN.species.std_initial_values>0);
n_species = size(CRN.matrix.S, 1);
ind_one = n_species + 1;

max_counter = 250;
toll_cond_init_point = 10^17; % This should be the same inside f_NLPC_restart

jacobian_v = f_compute_analytic_jacobian_v(v, n_species, ind_one);

lof_mutations = {'APC', 'AKT', 'SMAD4',  'PTEN'};
gof_mutations = {'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = [lof_mutations, gof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);

%% Step 3. Run NLPC on the network for physiological cells
if do_phys
    
% 3.1. Define network parameters specific to physiological cells
rates_phys = CRN.rates.std_values;
S_phys = CRN.matrix.S;
rho_phys = Nl*CRN.species.std_initial_values;

% 3.2. Load initial points
load(file_x0_phys, 'x0_all');
n_runs = size(x0_all, 2);

% 3.3. Run the algorithm
for ir = 1:n_runs
    fprintf('Physiological run = %d \n', ir)
    time_init = tic;
    aux_phys = f_NLPC_restart(x0_all(:, ir), rates_phys, S_phys, Nl, ...
        rho_phys, idx_basic_species, v, ind_one, max_counter, 0);
    aux_phys.elapse_time = toc(time_init);
    nlpc_phys(ir) = aux_phys;

    clear aux_phys tim_init
end

% 3.4. Save
save(fullfile(folder_results, 'nlpc_phys.mat'), 'nlpc_phys');
clear nlpc_phys_aml x0_all rates_phys S_phys rho_phys n_runs

end
%% Step 4. Run NLPC on the network of mutated cells
if do_mutation 
    
    x0_phys = CRN.species.std_initial_values;
    rates_phys = CRN.rates.std_values;
    S_mut = CRN.matrix.S; % The code is setting such that GOF mutations
                          % are defined zeroing proper rates
    
for im = 1:n_mutations
    
    protein = all_mutations{im};
    
    % 4.1. Load parameters of mutated network and initial points
    file_x0_mut = fullfile(folder_st_points, sprintf(aux_file_x0_mut, protein));
    load(file_x0_mut, 'x0_all', 'par');
    rho_mut = par.rho_mut;
    rates_mut = par.rates_mut;
    n_runs = size(x0_all, 2);
    
    % 4.2. Run algorithm
    for ir = 1:n_runs
        fprintf('Mutation %s run = %d \n', protein, ir)
        time_init = tic;
        aux_mut = f_NLPC_restart(x0_all(:, ir), rates_mut, S_mut, Nl, ...
            rho_mut, idx_basic_species, v, ind_one, max_counter, 0);
        aux_mut.elapse_time = toc(time_init);
        nlpc_mut(ir) = aux_mut;
        clear aux_mut tim_init
    end

    % 4.4. Save
    save(fullfile(folder_results, ...
        sprintf('nlpc_mut_%s.mat', protein)), 'nlpc_mut', 'x0_all')

    clear protein nlpc_mut x0_all rates_mut rho_mut par n_runs

end

end