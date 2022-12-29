clc
clear
close all

% Here we compare NLPC method with the classical orthogonal
% projector and NLPC method with the novel non-linear projector.
% They are compared in terms of:
% - number of starting point the algorithm has to be applied with to have
%   convergence;
% - time elapsed to have convergence.

% Then we create a table where the two projectors are compared in terms of:
% - conditioning number of the jacobian function J_F at each iteration of
% the algorithm
% - percentage of null components of each reconstruction of the solution at
% each iteration of the algorithm


addpath('../funcs')

lof_mutations = {'APC',  'AKT', 'SMAD4', 'PTEN'};
gof_mutations = {'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = ['phys', gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);

%% Step 1. Define general parameters
% 1.1. Folder and files
folder_data = '../data';
folder_results = './results_paper';
folder_figures = './figures';
if ~exist(folder_figures, 'dir')
   mkdir(folder_figures)
end

file_CRN = fullfile(folder_data, 'CRC_CRN_nodrug.mat');

aux_nlpc_phys = 'nlpc_phys.mat';
aux_nlpc_ort_phys = 'nlpc_ort_phys.mat';
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

perc = 0; n_runs = 50; n_species = 419;

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
       load(fullfile(folder_results, sprintf(aux_nlpc_phys, mutation)), 'nlpc_phys')
       disp(fullfile(folder_results, sprintf(aux_nlpc_phys, mutation)))
       load(fullfile(folder_results, sprintf(aux_nlpc_ort_phys, mutation)), 'nlpc_ort_phys')
       
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
    load(fullfile(folder_results, sprintf(aux_nlpc_mut, mutation)), 'nlpc_mut');
    load(fullfile(folder_results, sprintf(aux_nlpc_ort_mut, mutation)), 'nlpc_ort_mut');
    
    bp_elapsed_time_nlpc(:, im) = [nlpc_mut.elapse_time];
    bp_elapsed_time_nlpc_ort(:, im) = [nlpc_ort_mut.elapse_time];
    
    bp_attempts_nlpc(:, im) = [nlpc_mut.num_trials];
    bp_attempts_nlpc_ort(:, im) = [nlpc_ort_mut.num_trials];
   
   end
   
    t = (im-1)*n_runs;
    ratio_attempts_all((t+1):(t+n_runs)) = bp_attempts_nlpc_ort(:,im)./bp_attempts_nlpc(:,im);
    
end

mean_ratio = mean(mean(ratio_attempts_all()));

disp("Mean of the ratio of number of attempts (PNG/PNG ort): " + mean_ratio);

% Change mutation labels
for im = 1:numel(cond_name)
    mutation = cond_name{im};
    if strcmp(mutation, 'Ras')
        cond_name(im)= {'k-Ras'};
    end
    if strcmp(mutation, 'BetaCatenin')
        cond_name(im) = {'Betacatenin'};
    end
    if strcmp(mutation, 'TP53')
        cond_name(im) = {'p53'};
    end
end

%% F1 (additional) Boxplots for elapsed time
bp_proj_time = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
aux_ = [bp_elapsed_time_nlpc; bp_elapsed_time_nlpc_ort];
group_idx = [ones(1, n_runs), 2*ones(1, n_runs)];
group_names = {'Novel proj. $\mathcal{P}$', 'Orth. proj. \textit{P}'};
h = daboxplot(aux_, 'groups', group_idx, 'outsymbol', 'ko',...
    'fill', 0,'legend', group_names(1:2), 'xtlabels', cond_name, 'mean', 0,...
    'color', [1 0 0; 0 0.5 0]);
grid on
ylabel('Elapsed time [sec]', 'Interpreter', 'Latex')
xtickangle(30)
set(gca, 'Fontsize', 20, 'TickLabelInterpreter','latex')
set(h.lg, 'Location', 'North', 'Interpreter', 'Latex')
ylim([-5, 2000])
saveas(bp_proj_time, fullfile(folder_figures, 'bp_elapsed_time_ort.png'))

%% F2 Boxplots for the number of restarts
f_bp_attempts = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
aux_ = [bp_attempts_nlpc; bp_attempts_nlpc_ort];
group_idx = [ones(1, n_runs), 2*ones(1, n_runs)];
group_names = {'Novel proj. $\mathcal{P}$', 'Orth. proj. \textit{P}'};
h = daboxplot(aux_,'groups', group_idx, 'outsymbol','ko',...
    'fill',0,'legend',group_names(1:2), 'xtlabels', cond_name, 'mean', 0,...
    'color', [1 0 0; 0 0.5 0]);
grid on
ylabel('N. Restarts', 'Interpreter', 'Latex')
xtickangle(30)
set(gca, 'Fontsize', 20, 'TickLabelInterpreter','latex')
set(h.lg, 'Location', 'North', 'Interpreter', 'Latex')
ylim([0, 68])
saveas(f_bp_attempts, fullfile(folder_figures, 'bp_restarts_ort.png'))


%% Table proj. comparison - physiological status


%% Mean max zeros 

max_zeri_ort = zeros(1, n_runs);
max_zeri_class = zeros(1, n_runs);

% Load
load(fullfile(folder_results, sprintf(aux_nlpc_phys, mutation)), 'nlpc_phys')
load(fullfile(folder_results, sprintf(aux_nlpc_ort_phys, mutation)), 'nlpc_ort_phys')
       
for i=1:n_runs
    zeri_x0 = sum(nlpc_ort_phys(i).x0 == 0);
    zeri_ort(i).n = [zeri_x0 nlpc_ort_phys(i).zeri(1).n];
    zeri_class(i).n = [zeri_x0 nlpc_phys(i).zeri(1).n];

    max_zeri_ort(i) = max(zeri_ort(i).n);
    max_zeri_class(i) = max(zeri_class(i).n);
end   
   
mean_max_zeri_ort = mean(max_zeri_ort);
mean_max_zeri_class = mean(max_zeri_class);
   
std_max_zeri_ort = std(max_zeri_ort);
std_max_zeri_class = std(max_zeri_class);

mean_max_zeri_ort_perc = mean_max_zeri_ort*(100/n_species);
std_max_zeri_ort_perc = std_max_zeri_ort*(100/n_species);
mean_max_zeri_class_perc = mean_max_zeri_class*(100/n_species);
std_max_zeri_class_perc = std_max_zeri_class*(100/n_species);


%% Conditioning number

max_cond_ort = zeros(1, n_runs);
max_cond_class = zeros(1, n_runs);
max_log_cond_ort = zeros(1, n_runs);
 
% Load
load(fullfile(folder_results, sprintf(aux_nlpc_phys, mutation)), 'nlpc_phys')
load(fullfile(folder_results, sprintf(aux_nlpc_ort_phys, mutation)), 'nlpc_ort_phys')
      
for i=1:n_runs
    cond_ort(i).n = [nlpc_ort_phys(i).cond_number(1).n];
    cond_class(i).n = [nlpc_phys(i).cond_number(1).n];
    log_cond_ort(i).n = log10(cond_ort(i).n);
    log_cond_class(i).n = log10(cond_class(i).n);
     
    max_cond_ort(i) = max(cond_ort(i).n);
    max_cond_class(i) = max(cond_class(i).n);
    max_log_cond_ort(i) = max(log_cond_ort(i).n);
    max_log_cond_class(i) = max(log_cond_class(i).n);
end

mean_max_cond_ort = mean(max_cond_ort);
mean_max_cond_class = mean(max_cond_class);   

mean_max_log_cond_ort = mean(max_log_cond_ort);
mean_max_log_cond_class = mean(max_log_cond_class);   

std_max_cond_ort = std(max_cond_ort);
std_max_cond_class = std(max_cond_class);   

std_max_log_cond_ort = std(max_log_cond_ort);
std_max_log_cond_class = std(max_log_cond_class);   


%% Create Latex table

string_table1 = "\begin{table}[ht]" + newline + " \centering " + newline + ...
    "	\begin{tabular}{c|*{2}{c|}}	\cline{2-3} " ...
    + "	& \textbf{Novel proj.} $\mathcal{P}$ & \textbf{Orth. proj. $P$} \\ " +...
    newline +	"\hline	" + newline + "\multicolumn{1}{|c|}{\textbf{Num. null components (\%)}} 	&";
string_table2 = " \end{tabular}	" + newline + "  \end{table}";
string_results = sprintf('%0.2f', mean_max_zeri_class_perc(1)) +" " + ' $\pm$ ' ...
    + " " + sprintf('%0.2f', std_max_zeri_class_perc(1)) ...
    + "& " + sprintf('%0.2f', mean_max_zeri_ort_perc(1)) +" " + ' $\pm$ ' + " " ...
    + sprintf('%0.2f', std_max_zeri_ort_perc(1)) + " \\" ...
    + "\hline" + newline ...
    + "\multicolumn{1}{|c|}{\textbf{Cond Number $\mathbf{J}_{\mathbf{f}}$ (log. scale)}} & " ...
    + sprintf('%d', round(mean_max_log_cond_class(1))) +" " ...
    + ' $\pm$ ' + " " + sprintf('%d', round(std_max_log_cond_class(1))) ...
    + "& " + sprintf('%d', round(mean_max_log_cond_ort(1))) +" " + ' $\pm$ ' + " " ...
    + sprintf('%d', round(std_max_log_cond_ort(1)));


string_results = string_results + "\\  " + newline + " \hline  " + newline + "";

string_table = string_table1 + string_results + string_table2;

disp(string_table)

