clear
close all
clc

addpath(fullfile('..', 'funcs'))

do_phys = 1;
do_mutation = 1; perc = 0;

%% Step1. Define path and load
target = fullfile('..', 'data');   
path_mim = fullfile(target, 'CRC_CRN_nodrug.mat');
load(path_mim, 'new_CMIM'); CRN = new_CMIM; 

%% Step 2. Define general parameters of the network
v = CRN.matrix.v;
Nl = CRN.matrix.Nl;
idx_basic_species = find(CRN.species.std_initial_values>0);
n_species = size(CRN.matrix.S, 1);
ind_one = n_species + 1;

n_runs = 50;
max_counter = 500;
toll_cond_init_point = 10^17; % This should be the same inside f_PNG_restart

jacobian_v = f_compute_analytic_jacobian_v(v, n_species, ind_one);

lof_mutations = {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT', 'ARF'};
gof_mutations = {'Raf', 'Ras', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = [gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);

%% Step 3. Run PNG on the network for physiological cells
if do_phys
    
% 3.1. Define network parameters specific to physiological cells
rates_phys = CRN.rates.std_values;
S_phys = CRN.matrix.S;
rho_phys = Nl*CRN.species.std_initial_values;

% 3.2. Define initial points
x0_all = zeros(n_species, n_runs);
for ir = 1:n_runs
    aux_cond = Inf; 
    while aux_cond > toll_cond_init_point
        aux = f_draw_from_ssurf(Nl, rho_phys, idx_basic_species, [-3, 3]);
        aux_cond = cond(f_evaluate_jacobian(rates_phys, aux, ...
                    S_phys, idx_basic_species, jacobian_v, Nl));
    end
    x0_all(:, ir) = aux;
end

% 3.3. Run algorithm
for ir = 1:n_runs
    fprintf('Physiological run = %d \n', ir)
    time_init = tic;
    aux_phys = f_PNG_restart(x0_all(:, ir), rates_phys, S_phys, Nl, ...
        rho_phys, idx_basic_species, v, ind_one, max_counter);
    aux_phys.elapse_time = toc(time_init);
    png_phys(ir) = aux_phys;
    clear aux_phys tim_init
end

% 3.4. Save
save('png_phys_test1.mat', 'png_phys', 'x0_all')

clear png_phys x0_all rates_phys S_phys rho_phys

end
%% Step 4. Run PNG on the network of mutated cells
if do_mutation 
    
    x0_phys = CRN.species.std_initial_values;
    rates_phys = CRN.rates.std_values;
    S_mut = CRN.matrix.S; % The code is setting such that GOF mutations
                          % are defined zeroing proper rates
    
for im = 1:n_mutations
    
    protein = all_mutations{im};
    
    % 4.1. Define network parameters specific to physiological cells
    [x0_mut, rates_mut] = f_define_mutated_condition(protein, ...
                                    x0_phys, rates_phys, CRN, perc);
    rho_mut = Nl*x0_mut;
    
    % 4.2. Define initial points
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
    
    % 4.3. Run algorithm
    for ir = 1:n_runs
        fprintf('Mutation %s run = %d \n', protein, ir)
        time_init = tic;
        aux_mut = f_PNG_restart(x0_all(:, ir), rates_mut, S_mut, Nl, ...
            rho_mut, idx_basic_species, v, ind_one, max_counter);
        aux_mut.elapse_time = toc(time_init);
        png_mut(ir) = aux_mut;
        clear aux_mut tim_init
    end

    % 4.4. Save
    save(sprintf('png_mut_%s_test1.mat', protein), 'png_mut', 'x0_all')

    clear protein png_mut x0_all rates_mut rho_mut

end

end

%% Step 5. Check results (physiological)
% path_results_dyn = fullfile('..', 'results');
% 
% file_phys = fullfile(path_results_dyn, 'results_physiological.mat');
% load(file_phys, 'ris_phys')
% ir = 1; 
% delta = (png_phys(ir).x  - ris_phys.x_eq(1:end-2)) ./ ris_phys.x_eq(1:end-2);
% figure
% plot(delta)

%% Step 6. Check results (mutation)
protein = 'BetaCatenin';
ir = 1; % Run of PNG to be consider

path_results_dyn = fullfile('..', 'results');

file_dyn = fullfile('..', 'results', ...
    sprintf('results_mutation_%s_perc_0.0.mat', protein));
file_png = fullfile(sprintf('png_mut_%s_test1.mat', protein));
file_dyn_phys = fullfile('..', 'results', ...
    'results_physiological.mat');
file_png_phys = fullfile('png_phys_test1.mat');

load(file_dyn, 'ris_mutated')
load(file_png, 'png_mut')
load(file_dyn_phys, 'ris_phys')
load(file_png_phys, 'png_phys')

x_eq_png = png_mut(ir).x;
x_eq_dyn_1 = ris_mutated.x_mut_eq(1:n_species);
x_eq_dyn_2 = ris_mutated.x_xemut_eq(1:n_species);
x_eq_phys_dyn = ris_phys.x_eq(1:n_species);
x_eq_phys = png_phys(1).x;

delta_png = (x_eq_png - x_eq_phys) ./ x_eq_phys;
delta_dyn_1 = (x_eq_dyn_1 - x_eq_phys_dyn) ./ x_eq_phys_dyn;
delta_dyn_2 = (x_eq_dyn_2 - x_eq_phys) ./ x_eq_phys;

figure
subplot(2, 1, 1)
plot(delta_png, 'k', 'linewidth', 2)
hold on
plot(delta_dyn_1, 'r--', 'linewidth', 1.5)
symlog('y')

subplot(2, 1, 2)
plot(delta_png, 'k', 'linewidth', 2)
hold on
plot(delta_dyn_2, 'r--', 'linewidth', 1.5)
symlog('y')

%% Save results

% save(fullfile(folder_results, '20211015_test_newtonGD_restart.mat'), 'struct')
