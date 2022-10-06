clc
clear
close all


% Here we create boxplots for comparing the PNG method (associated to non-projector
% operator) with the dynamic approach.
% They are compared in terms of:
% - precision, i.e. norm of F in the equilibrium points;
% - time elapsed to have convergence.


addpath('../funcs')

lof_mutations = {'APC', 'AKT', 'SMAD4',  'PTEN'};
gof_mutations = {'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = ['phys', gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);

%% Step 1. Define general parameters
% 1.1. Folder and files
folder_data = '../data';
folder_results = './results';
folder_results_dyn = './results/dinamica';
folder_figures = './results/figures';

file_CRN = fullfile(folder_data, 'CRC_CRN_nodrug.mat');

aux_dyn_phys = 'dyn_%s.mat';
aux_png_phys = 'png_%s_new.mat';
aux_dyn_mut = 'dyn_mut_%s.mat';
aux_png_mut = 'png_mut_new_%s.mat';

% 1.2. Define general parameters of the network
load(file_CRN); CRN = new_CMIM; clear new_CMIM
rates_phys = CRN.rates.std_values;
x0_phys = CRN.species.std_initial_values;
idx_basic_species = find(x0_phys>0);
Nl = CRN.matrix.Nl;
S_phys = CRN.matrix.S;
v_phys = CRN.matrix.v;
ind_one = size(CRN.matrix.S, 1) + 1;

perc = 0; n_runs = 50;

%% Step 2. Make boxplots
% 2.1. Initialize
bp_elapsed_time_png = zeros(n_runs, n_mutations);
bp_elapsed_time_dyn = zeros(n_runs, n_mutations);
bp_fxeq_png = zeros(n_runs, n_mutations);
bp_fxeq_dyn = zeros(n_runs, n_mutations);


% 
for im = 1:n_mutations
   mutation = all_mutations{im};
   cond_name{im} = mutation;

   %% Physiological network
   if strcmp(mutation, 'phys')
       
       rho_phys = Nl*x0_phys; 
       
       % Load
       load(fullfile(folder_results_dyn, sprintf(aux_dyn_phys, mutation)), 'dyn_phys')
       load(fullfile(folder_results, sprintf(aux_png_phys, mutation)), 'png_phys')
       
       % Store results
       bp_elapsed_time_dyn(:, im) = [dyn_phys.elapse_time]; 
       bp_elapsed_time_png(:, im) = [png_phys.elapse_time];
       
       f_dyn = zeros(n_runs, 1);
       f_png = zeros(n_runs, 1);
       for ir = 1:n_runs
            f_dyn(ir) = norm(f_evaluate_mim(rates_phys, dyn_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one));
            f_png(ir) = norm(f_evaluate_mim(rates_phys, png_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one));                        
       end
       bp_fxeq_dyn(:, im) = f_dyn;
       bp_fxeq_png(:, im) = f_png;
       
       clear dyn_phys png_phys
       
   else
   
   %% Mutated network
   
    % Define parameters of the newtork
    [x0_mut, rates_mut] = f_define_mutated_condition(mutation, ...
                                x0_phys, rates_phys, CRN, perc);
    rho_mut = Nl*x0_mut;
    % Load
    load(fullfile(folder_results_dyn, sprintf(aux_dyn_mut, mutation)), 'dyn_mut');
    load(fullfile(folder_results, sprintf(aux_png_mut, mutation)), 'png_mut');
    
    bp_elapsed_time_dyn(:, im) = [dyn_mut.elapse_time];
    bp_elapsed_time_png(:, im) = [png_mut.elapse_time];
    for ir = 1:n_runs
        f_dyn(ir) = norm(f_evaluate_mim(rates_mut, dyn_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one));
        f_png(ir) = norm(f_evaluate_mim(rates_mut, png_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one));                        
    end
    bp_fxeq_dyn(:, im) = f_dyn;
    bp_fxeq_png(:, im) = f_png;
   
   end
   
end

%% Make boxplots
f_bp_time = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
aux_ = [bp_elapsed_time_png; bp_elapsed_time_dyn];
group_idx = [ones(1, n_runs), 2*ones(1, n_runs)];
group_names = {'PNG', 'Dynamics'};
h = daboxplot(aux_,'groups', group_idx, 'outsymbol','ko',...
    'fill',0,'legend',group_names(1:2), 'xtlabels', cond_name, 'mean', 0,...
    'color', [1 0 0; 0 0.447 0.7410]);
grid on
ylabel('Elapsed time [sec]')
xtickangle(30)
set(gca, 'Fontsize', 20, 'Yscale', 'linear')%, 'log')
ylim([-5, 601])
saveas(f_bp_time, fullfile(folder_figures, 'bp_time.png'))

f_bp_fxe = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
aux_ = [bp_fxeq_png; bp_fxeq_dyn];
group_idx = [ones(1, n_runs), 2*ones(1, n_runs)];
group_names = {'PNG', 'Dynamics'};
h = daboxplot(aux_,'groups', group_idx, 'outsymbol','ko',...
    'fill',0,'legend',group_names(1:2), 'xtlabels', cond_name, 'mean', 0,...
    'color', [1 0 0; 0 0.447 0.7410]);
grid on
ylabel('||f(x_e)||')
xtickangle(30)
set(gca, 'Fontsize', 20, 'YScale', 'log')
set(h.lg, 'Location', 'East')
ylim([0, 1])
saveas(f_bp_fxe, fullfile(folder_figures, 'bp_fxe.png'))

%%
mean_time_png = mean(bp_elapsed_time_png, 1);
[fast_value, fast_idx] = min(mean_time_png);
[slow_value, slow_idx] = max(mean_time_png);

mean_time_dyn = mean(bp_elapsed_time_dyn, 1);
[fast_value_dyn, fast_idx_dyn] = min(mean_time_dyn);
[slow_value_dyn, slow_idx_dyn] = max(mean_time_dyn);




