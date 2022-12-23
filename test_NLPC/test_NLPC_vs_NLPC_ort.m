clc
clear
close all


% Here we create boxplots for comparing NLPC method with orthogonal
% projector with the NLPC method with non-projector operator.
% They are compared in terms of:
% - number of starting point the algorithm has to be applied with to have
%   convergence;
% - time elapsed to have convergence.


addpath('../funcs')

lof_mutations = {'APC',  'AKT', 'SMAD4', 'PTEN'};
gof_mutations = {'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = ['phys', gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);

%% Step 1. Define general parameters
% 1.1. Folder and files
folder_data = '../data';
folder_results1 = './results/31_10_orthproj_40alpha_10-2';
folder_results2 = './results/31_10_classicmethod3mod_40alpha_10-2';
folder_figures = './results/figures/31_10';

file_CRN = fullfile(folder_data, 'CRC_CRN_nodrug.mat');

aux_nlpc_phys = 'nlpc_%s.mat';
aux_nlpc_ort_phys = 'nlpc_ort_%s.mat';
aux_nlpc_mut = 'nlpc_mut_%s.mat';
aux_nlpc_ort_mut = 'nlpc_ort_mut_%s.mat';

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
bp_elapsed_time_nlpc = zeros(n_runs, n_mutations);
bp_elapsed_time_nlpc_ort = zeros(n_runs, n_mutations);
bp_attempts_nlpc = zeros(n_runs, n_mutations);
bp_attempts_nlpc_ort = zeros(n_runs, n_mutations);

mean_attempts_nlpc = zeros(1, n_mutations);
mean_attempts_nlpc_ort = zeros(1, n_mutations);

ratio_attempts = zeros(n_runs);
mean_ratio_attempts = zeros(1, n_mutations);

ratio_attempts_all = zeros(1, n_runs * n_mutations);

% 
for im = 1:n_mutations
   mutation = all_mutations{im};
   cond_name{im} = mutation;

   %% Physiological network
   if strcmp(mutation, 'phys')
       
       rho_phys = Nl*x0_phys; 
       
       % Load
       load(fullfile(folder_results2, sprintf(aux_nlpc_phys, mutation)), 'nlpc_phys')
       disp(fullfile(folder_results2, sprintf(aux_nlpc_phys, mutation)))
       load(fullfile(folder_results1, sprintf(aux_nlpc_ort_phys, mutation)), 'nlpc_ort_phys')
       
       % Store results
       bp_elapsed_time_nlpc(:, im) = [nlpc_phys.elapse_time]; 
       bp_elapsed_time_nlpc_ort(:, im) = [nlpc_ort_phys.elapse_time];
       
       bp_attempts_nlpc(:, im) = [nlpc_phys.num_trials];
       bp_attempts_nlpc_ort(:, im) = [nlpc_ort_phys.num_trials];
       
       clear nlpc_phys nlpc_ort_phys
       
   else
   
   %% Mutated network
   
    % Define parameters of the network
    [x0_mut, rates_mut] = f_define_mutated_condition(mutation, ...
                                x0_phys, rates_phys, CRN, perc);
    rho_mut = Nl*x0_mut;
    % Load
    load(fullfile(folder_results2, sprintf(aux_nlpc_mut, mutation)), 'nlpc_mut');
    load(fullfile(folder_results1, sprintf(aux_nlpc_ort_mut, mutation)), 'nlpc_ort_mut');
    
    bp_elapsed_time_nlpc(:, im) = [nlpc_mut.elapse_time];
    bp_elapsed_time_nlpc_ort(:, im) = [nlpc_ort_mut.elapse_time];
    
    bp_attempts_nlpc(:, im) = [nlpc_mut.num_trials];
    bp_attempts_nlpc_ort(:, im) = [nlpc_ort_mut.num_trials];
   
   end
   
%    ratio_attempts = bp_attempts_nlpc_ort(:,im)./bp_attempts_nlpc(:,im);
%    mean_ratio_attempts(im) = mean(ratio_attempts);
%    
%    mean_attempts_nlpc(im) = mean(bp_attempts_nlpc(:,im));
%    mean_attempts_nlpc_ort(im) = mean(bp_attempts_nlpc_ort(:,im));
%    
    t = (im-1)*n_runs;
    ratio_attempts_all((t+1):(t+n_runs)) = bp_attempts_nlpc_ort(:,im)./bp_attempts_nlpc(:,im);
    
end


%% Make boxplots 
for im = 1:numel(cond_name)
    mutation = cond_name{im};
    if strcmp(mutation, 'Ras')
        cond_name(im) = {'k-Ras'};
    end
    if strcmp(mutation, 'BetaCatenin')
        cond_name(im) = {'Betacatenin'};
    end
    if strcmp(mutation, 'TP53')
        cond_name(im) = {'p53'};
    end
end

bp_proj_time = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
aux_ = [bp_elapsed_time_nlpc; bp_elapsed_time_nlpc_ort];
group_idx = [ones(1, n_runs), 2*ones(1, n_runs)];
group_names = {'Operator $\mathcal{P}$', 'Projector \textit{P}'};
h = daboxplot(aux_, 'groups', group_idx, 'outsymbol', 'ko',...
    'fill', 0,'legend', group_names(1:2), 'xtlabels', cond_name, 'mean', 0,...
    'color', [1 0 0; 0 0.5 0]);
grid on
ylabel('Elapsed time [sec]', 'Interpreter', 'Latex')
xtickangle(30)
set(gca, 'Fontsize', 20, 'TickLabelInterpreter','latex')
set(h.lg, 'Location', 'North', 'Interpreter', 'Latex')
ylim([-5, 2000])
saveas(bp_proj_time, fullfile(folder_figures, 'bp_proj_time2.png'))

f_bp_attempts = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
aux_ = [bp_attempts_nlpc; bp_attempts_nlpc_ort];
group_idx = [ones(1, n_runs), 2*ones(1, n_runs)];
group_names = {'Operator $\mathcal{P}$', 'Projector \textit{P}'};
h = daboxplot(aux_,'groups', group_idx, 'outsymbol','ko',...
    'fill',0,'legend',group_names(1:2), 'xtlabels', cond_name, 'mean', 0,...
    'color', [1 0 0; 0 0.5 0]);
grid on
ylabel('N. Restarts', 'Interpreter', 'Latex')
xtickangle(30)
set(gca, 'Fontsize', 20, 'TickLabelInterpreter','latex')
set(h.lg, 'Location', 'North', 'Interpreter', 'Latex')
ylim([0, 60])
saveas(f_bp_attempts, fullfile(folder_figures, 'bp_attempts2.png'))


%% Calcolo rapporto fra numero di restarts con l'ortogonale e con il "non-proiettore"

ratio = bp_attempts_nlpc_ort ./ bp_attempts_nlpc;
mean_ratio = mean(ratio(:));

disp("Mean of the ratio of number of attempts (NLPC/NLPC ort): " + mean_ratio);


%% 

mean_restarts_nlpc = mean(bp_attempts_nlpc);
mean_restarts_nlpc_ort = mean(bp_attempts_nlpc_ort);
ratio = mean_restarts_nlpc_ort ./ mean_restarts_nlpc;

mean_time_nlpc = mean(bp_elapsed_time_nlpc, 1);
[fast_value, fast_idx] = min(mean_time_nlpc);
[slow_value, slow_idx] = max(mean_time_nlpc);

mean_time_nlpc = mean(bp_elapsed_time_nlpc_ort, 1);
[fast_value_nlpc_ort, fast_idx_nlpc_ort] = min(mean_time_nlpc_ort);
[slow_value_nlpc_ort, slow_idx_nlpc_ort] = max(mean_time_nlpc_ort);




