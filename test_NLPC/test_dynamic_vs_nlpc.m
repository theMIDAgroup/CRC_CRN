clc
clear
close all


% Here we create boxplots for comparing the NLPC method (associated to non-projector
% operator) with the dynamic approach.
% They are compared in terms of:
% - precision, i.e. norm of F in the equilibrium points;
% - time elapsed to have convergence.


addpath('../funcs')

lof_mutations = {'APC', 'AKT', 'SMAD4', 'PTEN'};
gof_mutations = {'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = ['phys', gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);

n_runs = 50;

%% Step 1. Define general parameters
% 1.1. Folder and files
folder_data = '../data';
folder_results = './results/18_11';
folder_results_dyn = './results/dinamica';
folder_figures = './results/figures/18_11';

file_CRN = fullfile(folder_data, 'CRC_CRN_nodrug.mat');

aux_dyn_phys = 'dyn_%s.mat';
aux_nlpc_phys = '/nlpc_%s.mat';
aux_dyn_mut = 'dyn_mut_%s.mat';
aux_nlpc_mut = '/nlpc_mut_%s.mat';

% 1.2. Define general parameters of the network
load(file_CRN); CRN = new_CMIM; clear new_CMIM
rates_phys = CRN.rates.std_values;
x0_phys = CRN.species.std_initial_values;
idx_basic_species = find(x0_phys>0);
Nl = CRN.matrix.Nl;
S_phys = CRN.matrix.S;
v_phys = CRN.matrix.v;
ind_one = size(CRN.matrix.S, 1) + 1;


%% Step 2. Make boxplots
% 2.1. Initialize
bp_elapsed_time_nlpc = zeros(n_runs, n_mutations);
bp_elapsed_time_dyn = zeros(n_runs, n_mutations);
bp_fxeq_nlpc = zeros(n_runs, n_mutations);
bp_fxeq_dyn = zeros(n_runs, n_mutations);

diff = zeros(1, n_runs);
mim_dyn = 0;
mim_nlpc = 0;

for im = 1:n_mutations
   mutation = all_mutations{im};
   cond_name{im} = mutation;

   %% Physiological network
   if strcmp(mutation, 'phys')
       
       rho_phys = Nl*x0_phys; 
       
       % Load
       load(fullfile(folder_results_dyn, sprintf(aux_dyn_phys, mutation)), 'dyn_phys')
       load(fullfile(folder_results, sprintf(aux_nlpc_phys, mutation)), 'nlpc_phys')
       
       % Store results
       bp_elapsed_time_dyn(:, im) = [dyn_phys.elapse_time]; 
       bp_elapsed_time_nlpc(:, im) = [nlpc_phys.elapse_time];
       
       f_dyn = zeros(n_runs, 1);
       f_nlpc = zeros(n_runs, 1);
       for ir = 1:n_runs
            mim_dyn = f_evaluate_mim(rates_phys, dyn_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one);
            f_dyn(ir) = norm(mim_dyn);
            mim_nlpc = f_evaluate_mim(rates_phys, nlpc_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one);
            diff(ir) = norm(mim_dyn - mim_nlpc);
            f_nlpc(ir) = norm(mim_nlpc);                        
       end
       bp_fxeq_dyn(:, im) = f_dyn;
       bp_fxeq_nlpc(:, im) = f_nlpc;
       
       figure; semilogy(diff, '*');
       title ("Comparison between equilibria - " + mutation);
       
       saveas(im, fullfile(folder_figures, sprintf('diff_eq_%s_dyn_nlpc.png', mutation)))
    
       mim_dyn = 0;
       mim_nlpc = 0;
       
       clear dyn_phys nlpc_phys
       
   else
   
   %% Mutated network
   
    % Define parameters of the newtork
    perc = 0; 
    [x0_mut, rates_mut] = f_define_mutated_condition(mutation, ...
                                x0_phys, rates_phys, CRN, perc);
    rho_mut = Nl*x0_mut;
    % Load
    load(fullfile(folder_results_dyn, sprintf(aux_dyn_mut, mutation)), 'dyn_mut');
    load(fullfile(folder_results, sprintf(aux_nlpc_mut, mutation)), 'nlpc_mut');
    
    bp_elapsed_time_dyn(:, im) = [dyn_mut.elapse_time];
    bp_elapsed_time_nlpc(:, im) = [nlpc_mut.elapse_time];
    for ir = 1:n_runs
        mim_dyn = f_evaluate_mim(rates_mut, dyn_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one);
        f_dyn(ir) = norm(mim_dyn);
        mim_nlpc = f_evaluate_mim(rates_mut, nlpc_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one);
        f_nlpc(ir) = norm(mim_nlpc);  
        diff(ir) = norm(mim_dyn - mim_nlpc);
    end
    bp_fxeq_dyn(:, im) = f_dyn;
    bp_fxeq_nlpc(:, im) = f_nlpc;
    
    im = figure; semilogy(diff, '*');
    title ("Comparison between equilibria - " + mutation);
    
    saveas(im, fullfile(folder_figures, sprintf('diff_eq_%s_dyn_NLPC.png', mutation)))
    
    mim_dyn = 0;
    mim_nlpc = 0;
   end
   
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
   

f_bp_time = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
aux_ = [bp_elapsed_time_nlpc; bp_elapsed_time_dyn];
group_idx = [ones(1, n_runs), 2*ones(1, n_runs)];
group_names = {'NLPC', 'Dynamics'};
h = daboxplot(aux_,'groups', group_idx, 'outsymbol','ko',...
    'fill',0,'legend',group_names(1:2), 'xtlabels', cond_name, 'mean', 0,...
    'color', [1 0 0; 0 0.447 0.7410]);
grid on
ylabel('Elapsed time [sec]', 'Interpreter', 'Latex')
xtickangle(30)
set(gca, 'Fontsize', 20, 'Yscale', 'linear', 'TickLabelInterpreter','latex')
set(h.lg, 'Interpreter', 'Latex')
ylim([-5, 601])
saveas(f_bp_time, fullfile(folder_figures, 'bp_time_grad2.png'))

f_bp_fxe = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
aux_ = [bp_fxeq_nlpc; bp_fxeq_dyn];
group_idx = [ones(1, n_runs), 2*ones(1, n_runs)];
group_names = {'NLPC', 'Dynamics'};
h = daboxplot(aux_,'groups', group_idx, 'outsymbol','ko',...
    'fill',0,'legend',group_names(1:2), 'xtlabels', cond_name, 'mean', 0,...
    'color', [1 0 0; 0 0.447 0.7410]);
grid on
ylabel('$|| \textbf{f}(\textbf{x}_{eq}) ||$', 'Interpreter', 'Latex')
xtickangle(30)
set(gca, 'Fontsize', 20, 'YScale', 'log', 'TickLabelInterpreter','latex')
set(h.lg, 'Location', 'East', 'Interpreter', 'Latex')
ylim([0, 1])
saveas(f_bp_fxe, fullfile(folder_figures, 'bp_fxe_grad2.png'))





