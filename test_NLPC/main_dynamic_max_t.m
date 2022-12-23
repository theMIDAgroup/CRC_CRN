
clc
clear
close all


% Here we use the dynamic approach for computing equilibrium points
% considering different time windows, to see if it's better to study the
% dynamical system at higher times to have more precise results.
% max_t considered:
% - 1*10^7
% - 2*10^7
% - 3*10^7
% - 4*10^7
% - 5*10^7
% PHYSIOLOGICAL CASE: results stored in 'results/dyn_phys_maxt.mat', where
% 'maxt' represents the value considered for max_t each time.
% There is no mutated case.


addpath(fullfile('..', 'funcs'))

do_phys = 1;
do_mutation = 0;

%% Step1. Define path and load
target = fullfile('..', 'data');
folder_results = './results';

path_mim = fullfile(target, 'CRC_CRN_nodrug.mat'); % Network
load(path_mim, 'new_CMIM'); CRN = new_CMIM;

file_x0_phys = fullfile(folder_results, 'x0_phys.mat');
aux_file_x0_mut = 'x0_%s.mat';

%% Step 2. Define general parameters of the network
x0_phys = CRN.species.std_initial_values;
rates_phys = CRN.rates.std_values;
max_t = [1 2 3 4 5]*10^7;
perc = 0;

lof_mutations = {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT', 'ARF'};
gof_mutations = { 'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = [gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);

%% Step 3. Solve the dynamical system for the physiological network
if do_phys
    
    for i=1:length(max_t)
        % 3.1. Load initial condition
        load(file_x0_phys, 'x0_all');
        n_runs = size(x0_all, 2);

        % 3.2. Solve the system for all the defined initial conditions
        for ir = 1:n_runs
            fprintf('Physiological run = %d \n', ir)
            time_init = tic;
            [~, aux_sol] = ode15s(@(t_, x_) f_odefun_MIM(...
                t_, x_, rates_phys, CRN, 'Sv'), [0 max_t(i)], x0_all(:, ir));
            elapse_time = toc(time_init);
            aux.x0 = x0_all(:, ir);
            aux_sol = aux_sol';
            aux.x = aux_sol(:, end);
            aux.elapse_time = elapse_time;
            dyn_phys_t(ir) = aux;

            clear aux aux_sol elapse_time time_init

        end

        % 3.3. Save results
        save(sullfile(folder_results, sprintf('dyn_phys_%s.mat', max_t(i))), 'dyn_phys_t', 'x0_all')

        clear x0_all dyn_phys n_runs
        end
end

