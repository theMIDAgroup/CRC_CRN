clc
clear
close all


% Here we compare the NLPC-BB method (associated to the new non-linear 
% projector) with the NLPC one and the dynamic approach.
% The three are compared in terms of elapsed time to have convergence;
% then the two versions of NLPC are compared in terms of accuracy, i.e.
% norm of F in the equilibrium points, with respect to time elapsed to reach
% equilibrium.


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
folder_results = './results_paper';
folder_figures = './figures';
if ~exist(folder_figures, 'dir')
   mkdir(folder_figures)
end

file_CRN = fullfile(folder_data, 'CRC_CRN_nodrug.mat');

aux_dyn_phys = 'dyn_%s.mat';
aux_nlpc_phys = 'nlpc_%s.mat';
aux_nlpc_bb_phys = 'nlpc_bb_%s.mat';
aux_dyn_mut = 'dyn_mut_%s.mat';
aux_nlpc_mut = 'nlpc_mut_%s.mat';
aux_nlpc_bb_mut = 'nlpc_bb_mut_%s.mat';

% 1.2. Define general parameters of the network
load(file_CRN); CRN = new_CMIM; clear new_CMIM
rates_phys = CRN.rates.std_values;
x0_phys = CRN.species.std_initial_values;
idx_basic_species = find(x0_phys>0);
Nl = CRN.matrix.Nl;
S_phys = CRN.matrix.S;
v_phys = CRN.matrix.v;
ind_one = size(CRN.matrix.S, 1) + 1;


%% Step 2. Make boxplots - NLPC vs NLPC-BB vs Dynamics
% 2.1. Initialize
bp_elapsed_time_nlpc = zeros(n_runs, n_mutations);
bp_elapsed_time_nlpc_bb = zeros(n_runs, n_mutations);
bp_elapsed_time_dyn = zeros(n_runs, n_mutations);
bp_fxeq_nlpc = zeros(n_runs, n_mutations);
bp_fxeq_nlpc_bb = zeros(n_runs, n_mutations);
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
       load(fullfile(folder_results, sprintf(aux_dyn_phys, mutation)), 'dyn_phys')
       load(fullfile(folder_results, sprintf(aux_nlpc_bb_phys, mutation)), 'nlpc_phys'); nlpc_bb_phys = nlpc_phys; clear nlpc_phys
       load(fullfile(folder_results, sprintf(aux_nlpc_phys, mutation)), 'nlpc_phys')

       % Store results
       bp_elapsed_time_dyn(:, im) = [dyn_phys.elapse_time]; 
       bp_elapsed_time_nlpc(:, im) = [nlpc_phys.elapse_time];
       bp_elapsed_time_nlpc_bb(:, im) = [nlpc_bb_phys.elapse_time];
       
       f_dyn = zeros(n_runs, 1);
       f_nlpc = zeros(n_runs, 1);
       f_nlpc_bb = zeros(n_runs, 1);
       for ir = 1:n_runs
            mim_dyn = f_evaluate_mim(rates_phys, dyn_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one);
            f_dyn(ir) = norm(mim_dyn);
            mim_nlpc = f_evaluate_mim(rates_phys, nlpc_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one);
            diff(ir) = norm(mim_dyn - mim_nlpc);
            f_nlpc(ir) = norm(mim_nlpc);  
            mim_nlpc_bb = f_evaluate_mim(rates_phys, nlpc_bb_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one);
            f_nlpc_bb(ir) = norm(mim_nlpc_bb);
       end
       bp_fxeq_dyn(:, im) = f_dyn;
       bp_fxeq_nlpc(:, im) = f_nlpc;
       
 
       mim_dyn = 0;
       mim_nlpc = 0;
       mim_nlpc_bb = 0;
       
       clear dyn_phys nlpc_phys
       
   else
   
   %% Mutated network
   
    % Define parameters of the newtork
    perc = 0; 
    [x0_mut, rates_mut] = f_define_mutated_condition(mutation, ...
                                x0_phys, rates_phys, CRN, perc);
    rho_mut = Nl*x0_mut;
    % Load
    load(fullfile(folder_results, sprintf(aux_dyn_mut, mutation)), 'dyn_mut');
    load(fullfile(folder_results, sprintf(aux_nlpc_bb_mut, mutation)), 'nlpc_mut'); nlpc_bb_mut = nlpc_mut; clear nlpc_mut
    load(fullfile(folder_results, sprintf(aux_nlpc_mut, mutation)), 'nlpc_mut'); 

    bp_elapsed_time_dyn(:, im) = [dyn_mut.elapse_time];
    bp_elapsed_time_nlpc(:, im) = [nlpc_mut.elapse_time];
    bp_elapsed_time_nlpc_bb(:, im) = [nlpc_bb_mut.elapse_time];

    for ir = 1:n_runs
        mim_dyn = f_evaluate_mim(rates_mut, dyn_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one);
        f_dyn(ir) = norm(mim_dyn);
        mim_nlpc = f_evaluate_mim(rates_mut, nlpc_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one);
        f_nlpc(ir) = norm(mim_nlpc);  
        diff(ir) = norm(mim_dyn - mim_nlpc);
        mim_nlpc_bb = f_evaluate_mim(rates_mut, nlpc_bb_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one);
        f_nlpc_bb(ir) = norm(mim_nlpc_bb);
    end
    bp_fxeq_dyn(:, im) = f_dyn;
    bp_fxeq_nlpc(:, im) = f_nlpc;
    bp_fxeq_nlpc_bb(:, im) = f_nlpc_bb;
    
      
    mim_dyn = 0;
    mim_nlpc = 0;
    mim_nlpc_bb = 0;
   end
   
end

%% Step 3. Make figures
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
   
%% F1. Boxplots elapsed time
f_bp_time = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
aux_ = [bp_elapsed_time_nlpc; bp_elapsed_time_nlpc_bb; bp_elapsed_time_dyn];
group_idx = [ones(1, n_runs), 2*ones(1, n_runs), 3*ones(1, n_runs)];
group_names = {'NLPC', 'NLPC with BB', 'Dynamics'};
h = daboxplot(aux_,'groups', group_idx, 'outsymbol','ko',...
    'fill',0,'legend',group_names(1:3), 'xtlabels', cond_name, 'mean', 0,...
    'color', [1 0 0; 0.6350 0.0780 0.1840; 0 0.447 0.7410]);
grid on
ylabel('Elapsed time [sec]', 'Interpreter', 'Latex')
xtickangle(30)
set(gca, 'Fontsize', 20, 'Yscale', 'linear', 'TickLabelInterpreter','latex')
set(h.lg, 'Interpreter', 'Latex')
ylim([-5, 601])
% saveas(f_bp_time, fullfile(folder_figures, 'bp_time_grad_bb.png'))

%% F2. Scatter plot Accuracy vs Elapsed time - NLPC vs NLPC-BB
a = [2, 3, 9];
el_time_nlpc_a = bp_elapsed_time_nlpc(:,a);
el_time_nlpc_a = el_time_nlpc_a(:);
fxeq_nlpc_a = bp_fxeq_nlpc(:,a);
fxeq_nlpc_a = fxeq_nlpc_a(:);
a_leg_nlpc = '';
for i = 1:length(a)
    if i ~= 1
        a_leg_nlpc = a_leg_nlpc + ', ';
    end
    mutation = all_mutations{a(i)};
    if strcmp(mutation, 'Ras')
        all_mutations(a(i)) = {'k-Ras'};
    end
    if strcmp(mutation, 'BetaCatenin')
        all_mutations(a(i)) = {'Betacatenin'};
    end
    if strcmp(mutation, 'TP53')
        all_mutations(a(i)) = {'p53'};
    end
    a_leg_nlpc = a_leg_nlpc + string(all_mutations(a(i)));
end

b = 4;
el_time_nlpc_b = bp_elapsed_time_nlpc(:,b);
el_time_nlpc_b = el_time_nlpc_b(:);
fxeq_nlpc_b = bp_fxeq_nlpc(:,b);
fxeq_nlpc_b = fxeq_nlpc_b(:);
b_leg_nlpc = '';
for i = 1:length(b)
    if i ~= 1
        b_leg_nlpc = b_leg_nlpc + ', ';
    end
    mutation = all_mutations{b(i)};
    if strcmp(mutation, 'Ras')
        all_mutations(b(i)) = {'k-Ras'};
    end
    if strcmp(mutation, 'BetaCatenin')
        all_mutations(b(i)) = {'Betacatenin'};
    end
    if strcmp(mutation, 'TP53')
        all_mutations(b(i)) = {'p53'};
    end
    b_leg_nlpc = b_leg_nlpc + string(all_mutations(b(i)));
end

c = 1:10;
c = c(~ismember(c,a));
c = c(~ismember(c,b));
el_time_nlpc_c = bp_elapsed_time_nlpc(:,c);
el_time_nlpc_c = el_time_nlpc_c(:);
fxeq_nlpc_c = bp_fxeq_nlpc(:, c);
fxeq_nlpc_c = fxeq_nlpc_c(:);
c_leg_nlpc = '';
for i = 1:length(c)
    if i ~= 1
        c_leg_nlpc = c_leg_nlpc + ', ';
    end
    mutation = all_mutations{c(i)};
    if strcmp(mutation, 'Ras')
        all_mutations(c(i)) = {'k-Ras'};
    end
    if strcmp(mutation, 'BetaCatenin')
        all_mutations(c(i)) = {'Betacatenin'};
    end
    if strcmp(mutation, 'TP53')
        all_mutations(c(i)) = {'p53'};
    end
    c_leg_nlpc = c_leg_nlpc + string(all_mutations(c(i)));
end



el_time_nlpc_bb_a = bp_elapsed_time_nlpc_bb(:,a);
el_time_nlpc_bb_a = el_time_nlpc_bb_a(:);
fxeq_nlpc_bb_a = bp_fxeq_nlpc_bb(:,a);
fxeq_nlpc_bb_a = fxeq_nlpc_bb_a(:);

el_time_nlpc_bb_b = bp_elapsed_time_nlpc_bb(:,b);
el_time_nlpc_bb_b = el_time_nlpc_bb_b(:);
fxeq_nlpc_bb_b = bp_fxeq_nlpc_bb(:,b);
fxeq_nlpc_bb_b = fxeq_nlpc_bb_b(:);

el_time_nlpc_bb_c = bp_elapsed_time_nlpc_bb(:,c);
el_time_nlpc_bb_c = el_time_nlpc_bb_c(:);
fxeq_nlpc_bb_c = bp_fxeq_nlpc_bb(:, c);
fxeq_nlpc_bb_c = fxeq_nlpc_bb_c(:);


el_time_dyn_a = bp_elapsed_time_dyn(:,a);
el_time_dyn_a = el_time_dyn_a(:);
fxeq_dyn_a = bp_fxeq_dyn(:,a);
fxeq_dyn_a = fxeq_dyn_a(:);

el_time_dyn_b = bp_elapsed_time_dyn(:,b);
el_time_dyn_b = el_time_dyn_b(:);
fxeq_dyn_b = bp_fxeq_dyn(:,b);
fxeq_dyn_b = fxeq_dyn_b(:);

el_time_dyn_c = bp_elapsed_time_dyn(:,c);
el_time_dyn_c = el_time_dyn_c(:);
fxeq_dyn_c = bp_fxeq_dyn(:, c);
fxeq_dyn_c = fxeq_dyn_c(:);

%% 
f_prec_nlpc_vs_dyn = figure('units','normalized','outerposition',[0 0 1 1]);

col(1, :) = [0.8500, 0.3250, 0.0980];
col(2, :) = [0.4940, 0.1840, 0.5560];
col(3, :) = [0.9290, 0.6940, 0.1250];


subplot(2,2,1);
line1 = semilogy(el_time_nlpc_a, fxeq_nlpc_a,'diamond', 'Linewidth', 1, 'Markersize', 10, ...
    'DisplayName', a_leg_nlpc, 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', col(1,:)); 
hold on
line2 = semilogy(el_time_nlpc_b, fxeq_nlpc_b, '.', 'MarkerSize', 30, ...
   'DisplayName', b_leg_nlpc, 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', col(2,:));  
line3 = semilogy(el_time_nlpc_c, fxeq_nlpc_c, '+', 'Linewidth', 3, 'Markersize', 10, ...
    'DisplayName', c_leg_nlpc, 'MarkerFaceColor', col(3,:), 'MarkerEdgeColor', col(3,:)); 
grid on
hold off
title('\textbf{NLPC}', 'Interpreter','latex')
xlabel ('Elapsed time [sec]', 'Interpreter', 'Latex');
ylabel ('Accuracy - $||\textbf{f}(\textbf{x}_{nlpc})||$', 'Interpreter', 'Latex');
set(gca, 'Fontsize', 20, 'TickLabelInterpreter','latex')

subplot(2,2,2);
line1 = semilogy(el_time_dyn_a, fxeq_dyn_a,'diamond', 'Linewidth', 1, 'Markersize', 10, ...
      'DisplayName', a_leg_nlpc, 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', col(1,:)); 
hold on
line2 = semilogy(el_time_dyn_b, fxeq_dyn_b, '.', 'MarkerSize', 30, ...
     'DisplayName', b_leg_nlpc, 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', col(2,:));
line3 = semilogy(el_time_dyn_c, fxeq_dyn_c,'+', 'Linewidth', 3, 'Markersize', 10, ...
    'DisplayName', c_leg_nlpc, 'MarkerFaceColor', col(3,:), 'MarkerEdgeColor', col(3,:));
grid on
hold off
title('\textbf{Dynamics}', 'Interpreter','latex')
xlabel ('Elapsed time [sec]', 'Interpreter', 'Latex')
ylabel ('Accuracy - $||\textbf{f}(\textbf{x}_{dyn})||$', 'Interpreter', 'Latex')
set(gca, 'Fontsize', 20, 'TickLabelInterpreter','latex')

ax = subplot(2,5,8,'Visible','off');
axPos = ax.Position;
delete(ax)

% Construct a Legend with the data from the sub-plots
hL = legend([line1,line2,line3]);
% Move the legend to the position of the extra axes
hL.Position(2:3) = axPos(1:2);
hL.Position(2) = 0.36;
hL.Interpreter = 'Latex';

saveas(f_prec_nlpc_vs_dyn, fullfile(folder_figures, 'scatter_nlpc_dyn.png'))


f_prec_nlpc_vs_nlpcbb = figure('units','normalized','outerposition',[0 0 1 1]);

col(1, :) = [0.8500, 0.3250, 0.0980];
col(2, :) = [0.4940, 0.1840, 0.5560];
col(3, :) = [0.9290, 0.6940, 0.1250];


subplot(2,2,1);
line1 = semilogy(el_time_nlpc_a, fxeq_nlpc_a,'diamond', 'Linewidth', 1, 'Markersize', 10, ...
    'DisplayName', a_leg_nlpc, 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', col(1,:)); 
hold on
line2 = semilogy(el_time_nlpc_b, fxeq_nlpc_b, '.', 'MarkerSize', 30, ...
   'DisplayName', b_leg_nlpc, 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', col(2,:));  
line3 = semilogy(el_time_nlpc_c, fxeq_nlpc_c, '+', 'Linewidth', 3, 'Markersize', 10, ...
    'DisplayName', c_leg_nlpc, 'MarkerFaceColor', col(3,:), 'MarkerEdgeColor', col(3,:)); 
grid on
hold off
title('\textbf{NLPC}', 'Interpreter','latex')
xlabel ('Elapsed time [sec]', 'Interpreter', 'Latex');
ylabel ('Accuracy - $||\textbf{f}(\textbf{x}_{nlpc})||$', 'Interpreter', 'Latex');
set(gca, 'Fontsize', 20, 'TickLabelInterpreter','latex')

subplot(2,2,2);
line1 = semilogy(el_time_nlpc_bb_a, fxeq_nlpc_bb_a,'diamond', 'Linewidth', 1, 'Markersize', 10, ...
      'DisplayName', a_leg_nlpc, 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', col(1,:)); 
hold on
line2 = semilogy(el_time_nlpc_bb_b, fxeq_nlpc_bb_b, '.', 'MarkerSize', 30, ...
     'DisplayName', b_leg_nlpc, 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', col(2,:));
line3 = semilogy(el_time_nlpc_bb_c, fxeq_nlpc_bb_c,'+', 'Linewidth', 3, 'Markersize', 10, ...
    'DisplayName', c_leg_nlpc, 'MarkerFaceColor', col(3,:), 'MarkerEdgeColor', col(3,:));
grid on
hold off
title('\textbf{NLPC with BB}', 'Interpreter','latex')
xlabel ('Elapsed time [sec]', 'Interpreter', 'Latex')
ylabel ('Accuracy - $||\textbf{f}(\textbf{x}_{dyn})||$', 'Interpreter', 'Latex')
set(gca, 'Fontsize', 20, 'TickLabelInterpreter','latex')

ax = subplot(2,5,8,'Visible','off');
axPos = ax.Position;
delete(ax)

% Construct a Legend with the data from the sub-plots
hL = legend([line1,line2,line3]);
% Move the legend to the position of the extra axes
hL.Position(2:3) = axPos(1:2);
hL.Position(2) = 0.36;
hL.Interpreter = 'Latex';

saveas(f_prec_nlpc_vs_nlpcbb, fullfile(folder_figures, 'scatter_nlpc_nlpcbb.png'))

%% F3 (additional) Box plots for Accuracy
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






