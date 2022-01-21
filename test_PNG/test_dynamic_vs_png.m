clc
clear
close all

addpath('../funcs')

%% Step 1. Define general parameters
% 1.1. Folder and files
folder_data = '../data';
folder_res = './results';

file_CRN = fullfile(folder_data, 'CRC_CRN_nodrug.mat');

file_dyn = fullfile(folder_res, 'dyn_phys_test1.mat');
file_png = fullfile(folder_res, 'png_phys_test1.mat');

% 1.2. Define general parameters of the network
load(file_CRN); CRN = new_CMIM; clear new_CMIM

rates_phys = CRN.rates.std_values;
x0_phys = CRN.species.std_initial_values;
idx_basic_species = find(x0_phys>0);
Nl = CRN.matrix.Nl;
S_phys = CRN.matrix.S;
v_phys = CRN.matrix.v;
ind_one = size(CRN.matrix.S, 1) + 1;

% 1.3. Mutations to be considered
perc = 0;
lof_mutations = {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT', 'ARF'};
gof_mutations = {'Raf', 'Ras', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = [gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);

aux_dyn_file = 'dyn_mut_%s_test1.mat';
aux_png_file = 'png_mut_%s_test1.mat';

file_table = fullfile(folder_res, 'table_dyn_vs_png.txt');

%% Step 1. physiological cell
% 1.1. Load
load(file_dyn)
load(file_png)

% 1.2. Elapse time
elapsed_time_dyn = [dyn_phys.elapse_time];
elapsed_time_png = [png_phys.elapse_time];

%       Average and std over runs 
mean_elapsed_dyn = mean(elapsed_time_dyn);
std_elapsed_dyn = std(elapsed_time_dyn);
mean_elapsed_png = mean(elapsed_time_png);
std_elapsed_png = std(elapsed_time_png);

% 1.3. Value of f in the equilibrium point
n_runs = numel(dyn_phys);
f_dyn = zeros(n_runs, 1);
f_png = zeros(n_runs, 1);

rho_phys = Nl*dyn_phys(1).x0;

for ir = 1:n_runs
    f_dyn(ir) = norm(f_evaluate_mim(rates_phys, dyn_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one));
    f_png(ir) = norm(f_evaluate_mim(rates_phys, png_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one));                        
end

%       Average and std over runs
mean_f_dyn = mean(f_dyn);
std_f_dyn = std(f_dyn);
mean_f_png = mean(f_png);
std_f_png = std(f_png);

% 1.4. Plot
figure
subplot(2, 1, 1)
plot(1:numel(elapsed_time_dyn), elapsed_time_dyn, 'b*', 'linewidth', 2, ...
    'Displayname', 'Dynamic')
hold on
plot(1:numel(elapsed_time_png), elapsed_time_png, 'r*', 'linewidth', 2, ...
    'Displayname', 'PNG')
legend('show')
xlabel('Runs')
ylabel('Elsapsed time [sec]')
title('Physiological')
subplot(2, 1, 2)
plot(1:numel(elapsed_time_dyn), f_dyn, 'b*', 'linewidth', 2, ...
    'Displayname', 'Dynamic')
hold on
plot(1:numel(elapsed_time_png), f_png, 'r*', 'linewidth', 2, ...
    'Displayname', 'PNG')
legend('show')
xlabel('Runs')
ylabel('f(x_{eq})')

%       Table
fileID = fopen(file_table,'w');
fprintf(fileID, 'Physiological &  %2.2e $\\pm$ %2.2e & %2.2e $\\pm$ %2.2e & %2.2f $\\pm$ %2.2f & %2.2f $\\pm$ %2.2f \\\\ \n', ...
    mean_f_png, std_f_png, mean_f_dyn, std_f_dyn, ...
    mean_elapsed_png, std_elapsed_png, mean_elapsed_dyn, std_elapsed_dyn);

%% Step 2. Mutated cell
for im = 1:n_mutations
    protein = all_mutations{im}; 
    
% 2.1. Load
    file_mut_dyn = fullfile(folder_res, sprintf(aux_dyn_file, protein));
    file_mut_png = fullfile(folder_res, sprintf(aux_png_file, protein));

    load(file_mut_dyn, 'dyn_mut')
    load(file_mut_png, 'png_mut')

% 2.2. Define parameters of the newtork
    [x0_mut, rates_mut] = f_define_mutated_condition(protein, ...
                                x0_phys, rates_phys, CRN, perc);

% 2.3. Elapse time
    elapsed_time_dyn = [dyn_mut.elapse_time];
    elapsed_time_png = [png_mut.elapse_time];
    
%       Average and std over runs 
    mean_elapsed_dyn = mean(elapsed_time_dyn);
    std_elapsed_dyn = std(elapsed_time_dyn);
    mean_elapsed_png = mean(elapsed_time_png);
    std_elapsed_png = std(elapsed_time_png);    

% 2.4. Value of f in the equilibrium point
    n_runs = numel(dyn_mut);
    f_dyn = zeros(n_runs, 1);
    f_png = zeros(n_runs, 1);

    rho_mut = Nl*x0_mut;

    for ir = 1:n_runs
        f_dyn(ir) = norm(f_evaluate_mim(rates_mut, dyn_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one));
        f_png(ir) = norm(f_evaluate_mim(rates_mut, png_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one));                        
    end
    
    %       Average and std over runs
    mean_f_dyn = mean(f_dyn);
    std_f_dyn = std(f_dyn);
    mean_f_png = mean(f_png);
    std_f_png = std(f_png);

% 2.5. Plot
    f_mut = figure;
    subplot(2, 1, 1)
    plot(1:numel(elapsed_time_dyn), elapsed_time_dyn, 'b*', 'linewidth', 2, ...
        'Displayname', 'Dynamic')
    hold on
    plot(1:numel(elapsed_time_png), elapsed_time_png, 'r*', 'linewidth', 2, ...
        'Displayname', 'PNG')
    legend('show')
    xlabel('Runs')
    ylabel('Elsapsed time [sec]')
    title(sprintf('Mutation %s', protein))
    subplot(2, 1, 2)
    plot(1:numel(elapsed_time_dyn), f_dyn, 'b*', 'linewidth', 2, ...
        'Displayname', 'Dynamic')
    hold on
    plot(1:numel(elapsed_time_png), f_png, 'r*', 'linewidth', 2, ...
        'Displayname', 'PNG')
    legend('show')
    xlabel('Runs')
    ylabel('f(x_{eq})')

    %       Table
    fprintf(fileID, '%s & %2.2e $\\pm$ %2.2e & %2.2e $\\pm$ %2.2e & %2.2f $\\pm$ %2.2f & %2.2f $\\pm$ %2.2f \\\\ \n', ...
    protein, mean_f_png, std_f_png, mean_f_dyn, std_f_dyn, ...
    mean_elapsed_png, std_elapsed_png, mean_elapsed_dyn, std_elapsed_dyn);
    
    clear elapsed_time_dyn elapsed_time_png dyn_mut png_mut f_dyn f_png ...
        mean_elapsed_dyn mean_elapsed_png std_elapsed_dyn std_elapsed_dyn ...
        mean_f_dyn mean_f_png std_f_dyn std_f_png
    
    pause
    
end
 fclose(fileID);


% 
% %       Average and std over runs
% mean_f_dyn = mean(f_dyn);
% std_f_dyn = std(f_dyn);
% mean_f_png = mean(f_png);
% std_f_png = std(f_png);
% 
% % 1.4. Plot
% figure
% subplot(2, 1, 1)
% plot(1:numel(elapsed_time_dyn), elapsed_time_dyn, 'b*', 'linewidth', 2, ...
%     'Displayname', 'Dynamic')
% hold on
% plot(1:numel(elapsed_time_png), elapsed_time_png, 'r*', 'linewidth', 2, ...
%     'Displayname', 'PNG')
% legend('show')
% xlabel('Runs')
% ylabel('Elsapsed time [sec]')
% title('Physiological')
% subplot(2, 1, 2)
% plot(1:numel(elapsed_time_dyn), f_dyn, 'b*', 'linewidth', 2, ...
%     'Displayname', 'Dynamic')
% hold on
% plot(1:numel(elapsed_time_png), f_png, 'r*', 'linewidth', 2, ...
%     'Displayname', 'PNG')
% legend('show')
% xlabel('Runs')
% ylabel('f(x_{eq})')
% 
% %       Table
% fprintf('Physiological & %2.2f $\\pm$ %2.2f & %2.2f $\\pm$ %2.2f &  %2.2e $\\pm$ %2.2e & %2.2e $\\pm$ %2.2e &  \n ', ...
%     mean_elapsed_png, std_elapsed_png, mean_elapsed_dyn, std_elapsed_dyn, ...
%     mean_f_png, std_f_png, mean_f_dyn, std_f_dyn)

