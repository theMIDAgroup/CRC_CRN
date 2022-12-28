% Here we work on different stoichiometric surfaces: the one corresponding to
% the physiological status and others related to specific mutations.
% On each one of these surfaces, we select 50 points.
% Considered mutations are the following:
% Loss of Function: {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT', 'ARF'}
% Gain of Function: {'Raf', 'Ras', 'PI3K', 'BetaCatenin'}
% Loss of Function (type2): {'TP53'}.
% Computed points are saved in 'results/starting_points/'

clc
clear
close 

addpath(fullfile('..', 'funcs'))

n_runs = 50;

%% Step 1. Define general parameters and folders
% 1.a. Define parameters of the analysis
lof_mutations = {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT', 'ARF'};
gof_mutations = {'Raf', 'Ras', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = [gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);
perc = 0;

% 1.b. Define CRN path and load
target = fullfile('..', 'data');
folder_results = './results/starting_points';
if ~exist(folder_results, 'dir')
   mkdir(folder_results)
end
path_mim = fullfile(target, 'CRC_CRN_nodrug.mat');
load(path_mim, 'new_CMIM'); CRN = new_CMIM;

%% Step 2. Define general parameters of the network
v = CRN.matrix.v;
Nl = CRN.matrix.Nl;
idx_basic_species = find(CRN.species.std_initial_values>0);
n_species = size(CRN.matrix.S, 1);
ind_one = n_species + 1;

jacobian_v = f_compute_analytic_jacobian_v(v, n_species, ind_one);

%% Step 3. Define initial points for the physiological state
% 3.1. Define network parameters specific to physiological cells
rates_phys = CRN.rates.std_values;
S_phys = CRN.matrix.S;
rho_phys = Nl*CRN.species.std_initial_values;

% 3.2. Draw initial points
x0_all = zeros(n_species, n_runs);
toll_cond_init_point = 10^17;
for ir = 1:n_runs
    aux_cond = Inf; 
    while aux_cond > toll_cond_init_point
        aux = f_draw_from_ssurf(Nl, rho_phys, idx_basic_species, [-3, 3]);
        aux_cond = cond(f_evaluate_jacobian(rates_phys, aux, ...
                    S_phys, idx_basic_species, jacobian_v, Nl));
    end
    x0_all(:, ir) = aux;
end

save(fullfile(folder_results, 'x0_phys.mat'), 'x0_all')

clear x0_all

%% Step 4. Define initial points for the mutated states
x0_phys = CRN.species.std_initial_values;
rates_phys = CRN.rates.std_values;
S_mut = CRN.matrix.S;

for im = 1:n_mutations

    protein = all_mutations{im};

    % 4.1. Define network parameters specific to mutated cells
    [x0_mut, rates_mut] = f_define_mutated_condition(protein, ...
                                    x0_phys, rates_phys, CRN, perc);
    rho_mut = Nl*x0_mut;
    
    % 4.2. Draw initial points
    x0_all = zeros(n_species, n_runs);
    for ir = 1:n_runs
        aux_cond = Inf; 
        while aux_cond > toll_cond_init_point
            aux = f_draw_from_ssurf(Nl, rho_mut, idx_basic_species, [-3, 3]);
            aux_cond = cond(f_evaluate_jacobian(rates_mut, aux, ...
                        S_mut, idx_basic_species, jacobian_v, Nl));
        end
        x0_all(:, ir) = aux;
    end

    par.rho_mut = rho_mut;
    par.rates_mut = rates_mut;
    
    % Save
    save(fullfile(folder_results, sprintf('x0_%s.mat', protein)), ...
        'x0_all', 'par')
    clear x0_mut rho_mut x0_all par
end
