clear all
close all
clc

% Here we use the dynamic approach for computing equilibrium points using as
% initial conditions the points the PNG method managed to reach convergence with.
% So we are using the ones saved in './results/png_phys_testproj.mat'.
% The reason why we do this is showing that the results we obtain with the dynamic
% approach doesn't depend on the specific point we're starting from, but only
% on the stoichiometric surface we're working on. 

do_phys = 1;
do_mutation = 0;

%% Path, folders & load - physiological case

addpath(fullfile('..', 'funcs'))

target = fullfile('..', 'data');
folder_results = './results';

path_mim = fullfile(target, 'CRC_CRN_nodrug.mat'); % Network
load(path_mim, 'new_CMIM'); CRN = new_CMIM;

load(fullfile(folder_results, 'png_phys_testproj.mat'), 'png_phys');
aux_file_x0_mut = 'x0_%s.mat';


%% Define general parameters of the network

x0_phys = CRN.species.std_initial_values;
rates_phys = CRN.rates.std_values;
max_t = 2.5*10^7;
n_species = size(CRN.matrix.S, 1);

perc = 0;

lof_mutations = {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT', 'ARF'};
gof_mutations = { 'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = [gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);

%% Starting points for the algorithm - physiological case

x0_vector = zeros(n_species, length(png_phys));

for i=1:length(png_phys)
    x0_vector(:, i) = png_phys(i).x0;
end

n_runs = size(x0_vector, 2);


%% Solve the dynamical system for the physiological network

if do_phys
    

% Solve the system for all the defined initial conditions
for ir = 1:n_runs
    fprintf('Physiological run = %d \n', ir)
    time_init = tic;
    [~, aux_sol] = ode15s(@(t_, x_) f_odefun_MIM(...
        t_, x_, rates_phys, CRN, 'Sv'), [0 max_t], x0_vector(:, ir));
    elapse_time = toc(time_init);
    aux.x0 = x0_vector(:, ir);
    aux_sol = aux_sol';
    aux.x = aux_sol(:, end);
    aux.elapse_time = elapse_time;
    dyn_phys(ir) = aux;
    
    clear aux aux_sol elapse_time time_init
    
end


%% Save results
save('dyn_x0_png_phys.mat', 'dyn_phys', 'x0_vector')

clear x0_all dyn_phys n_runs

end


%% Step 4. Solve the dynamical system for the mutated network
if do_mutation
for im = 1:n_mutations

    protein = all_mutations{im};

% 4.1. Define the mutated network 
    file_x0_mut = fullfile(folder_results, sprintf(aux_file_x0_mut, protein));
    load(file_x0_mut, 'x0_all', 'par');
    
    x0_vector = zeros(n_species, length(png_phys));

    for i=1:length(png_phys)
        load(fullfile(folder_results, sprintf('png_mut_%s.mat', protein)));
        x0_vector(:, i) = png_mut.x0;
    end
    
    rho_mut = par.rho_mut;
    rates_mut = par.rates_mut;
    n_runs = size(x0_vector, 2);
    
% 4.3. Solve the system for all the defined initial conditions
    for ir = 1:n_runs
        fprintf('Mutation of %s run = %d \n', protein, ir)
        time_init = tic;
        [~, aux_sol] = ode15s(@(t_, x_) f_odefun_MIM(...
            t_, x_, rates_mut, CRN, 'Sv'), [0 max_t], x0_vector(:, ir));
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
        sprintf('dyn_x0_png_mut_%s.mat', protein)), 'dyn_mut')
    
    clear protein rates_mut x0_all dyn_mut par
    
end
end