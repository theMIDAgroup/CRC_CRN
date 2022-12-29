clc
clear
close all

% Here we use the dynamic approach for computing equilibria using the 
% points x0_all as initial conditions. For comparison purposes, we have 
% used the same points for as initial point of NLPC method 
% in 'main_NLPC.m' and 'main_NLPC_orthogonal.m' even though in that case
% some of those vectors didn't lead the algorithm to convergence.
% Those starting points have been selected with 'main_extract_x0.m'.

addpath(fullfile('..', 'funcs'))

do_phys = 1;
do_mutation = 1;

%% Step1. Define path and load
target = fullfile('..', 'data');
folder_results = './results';
folder_st_points = './results/starting_points';

path_mim = fullfile(target, 'CRC_CRN_nodrug.mat'); % Network
load(path_mim, 'new_CMIM'); CRN = new_CMIM;

file_x0_phys = fullfile(folder_st_points, 'x0_phys.mat');
aux_file_x0_mut = 'x0_%s.mat';

%% Step 2. Define general parameters of the network

x0_phys = CRN.species.std_initial_values;
rates_phys = CRN.rates.std_values;
max_t = 2.5*10^7;
perc = 0;

lof_mutations = {'APC', 'AKT', 'SMAD4',  'PTEN'};
gof_mutations = { 'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = [gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);


%% Step 3. Solve the dynamical system for the physiological network
if do_phys

% 3.1. Load initial condition
load(file_x0_phys, 'x0_all');
n_runs = size(x0_all, 2);

% 3.2. Solve the system for all the defined initial conditions
for ir = 1:n_runs
    fprintf('Physiological run = %d \n', ir)
    time_init = tic;
    [~, aux_sol] = ode15s(@(t_, x_) f_odefun_MIM(...
        t_, x_, rates_phys, CRN, 'Sv'), [0 max_t], x0_all(:, ir));
    elapse_time = toc(time_init);
    aux.x0 = x0_all(:, ir);
    aux_sol = aux_sol';
    aux.x = aux_sol(:, end);
    aux.elapse_time = elapse_time;
    dyn_phys(ir) = aux;
    
    clear aux aux_sol elapse_time time_init
    
end

% 3.3. Save results
save(fullfile(folder_results, 'dyn_phys.mat'), 'dyn_phys', 'x0_all')

clear x0_all dyn_phys n_runs
 
end

%% Step 4. Solve the dynamical system for the mutated network
if do_mutation
for im = 1:n_mutations

    protein = all_mutations{im};

% 4.1. Define the mutated network 
    file_x0_mut = fullfile(folder_st_points, sprintf(aux_file_x0_mut, protein));
    load(file_x0_mut, 'x0_all', 'par');
    rho_mut = par.rho_mut;
    rates_mut = par.rates_mut;
    n_runs = size(x0_all, 2);
    
% 4.3. Solve the system for all the defined initial conditions
    for ir = 1:n_runs
        fprintf('Mutation of %s run = %d \n', protein, ir)
        time_init = tic;
        [~, aux_sol] = ode15s(@(t_, x_) f_odefun_MIM(...
            t_, x_, rates_mut, CRN, 'Sv'), [0 max_t], x0_all(:, ir));
        elapse_time = toc(time_init);
        aux.x0 = x0_all(:, ir);
        aux_sol = aux_sol';
        aux.x = aux_sol(:, end);
        aux.elapse_time = elapse_time;
        dyn_mut(ir) = aux;

        clear aux aux_sol elapse_time time_init

    end
    
% 4.4. Save results
    save(fullfile(folder_results, ...
        sprintf('dyn_mut_%s.mat', protein)), 'dyn_mut')
    
    clear protein rates_mut x0_all dyn_mut par
    
end
end