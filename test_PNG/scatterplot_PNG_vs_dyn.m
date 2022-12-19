clc
clear
close all


% Scatterplot precision - time elapsed of dynamic approach vs PNG
% algorithm. Results obtained fixing 11 different networks (physiological
% and mutated ones)

addpath('../funcs')

lof_mutations = {'APC', 'AKT', 'SMAD4', 'PTEN'};
gof_mutations = {'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = ['phys', gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);

%% Step 1. Define general parameters
% 1.1. Folder and files
folder_data = '../data';
folder_results = './results/18_11';
folder_results_dyn = './results/dinamica';
folder_figures = './results/figures/31_10';

file_CRN = fullfile(folder_data, 'CRC_CRN_nodrug.mat');

aux_dyn_phys = 'dyn_%s.mat';
aux_png_phys = '/png_%s.mat';
aux_dyn_mut = 'dyn_mut_%s.mat';
aux_png_mut = '/png_mut_%s.mat';

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

perc = 0; n_runs = 50;

elapsed_time_png = zeros(n_runs, n_mutations); 
elapsed_time_dyn = zeros(n_runs, n_mutations);
fxeq_png = zeros(n_runs, n_mutations);
fxeq_dyn = zeros(n_runs, n_mutations);

diff = zeros(1, n_runs);
mim_dyn = 0;
mim_png = 0;

for im = 1:n_mutations
   mutation = all_mutations{im};

   %% Physiological network
   if strcmp(mutation, 'phys')
       
       rho_phys = Nl*x0_phys; 
       % Load
       load(fullfile(folder_results_dyn, sprintf(aux_dyn_phys, mutation)), 'dyn_phys')
       load(fullfile(folder_results, sprintf(aux_png_phys, mutation)), 'png_phys')
       % Store results
       elapsed_time_dyn(:, im) = [dyn_phys.elapse_time]; 
       elapsed_time_png(:, im) = [png_phys.elapse_time];
       
       f_dyn = zeros(n_runs, 1);
       f_png = zeros(n_runs, 1);
       for ir = 1:n_runs
            mim_dyn = f_evaluate_mim(rates_phys, dyn_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one);
            f_dyn(ir) = norm(mim_dyn);
            mim_png = f_evaluate_mim(rates_phys, png_phys(ir).x, idx_basic_species, ...
                                Nl, rho_phys, S_phys, v_phys, ind_one);
            diff(ir) = norm(mim_dyn - mim_png);
            f_png(ir) = norm(mim_png);                        
       end
       fxeq_dyn(:, im) = f_dyn;
       fxeq_png(:, im) = f_png;
       
       mim_dyn = 0;
       mim_png = 0;
       
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
    %Store results
    elapsed_time_dyn(:, im) = [dyn_mut.elapse_time];
    elapsed_time_png(:, im) = [png_mut.elapse_time];
    for ir = 1:n_runs
        mim_dyn = f_evaluate_mim(rates_mut, dyn_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one);
        f_dyn(ir) = norm(mim_dyn);
        mim_png = f_evaluate_mim(rates_mut, png_mut(ir).x, idx_basic_species, ...
                                    Nl, rho_mut, S_phys, v_phys, ind_one);
        f_png(ir) = norm(mim_png);  
        diff(ir) = norm(mim_dyn - mim_png);
    end
    fxeq_dyn(:, im) = f_dyn;
    fxeq_png(:, im) = f_png;
    
    mim_dyn = 0;
    mim_png = 0;
   end
   
end


%% Graphs

% Graph divided by times (all networks)

a = [2, 3, 9];
el_time_png_a = elapsed_time_png(:,a);
el_time_png_a = el_time_png_a(:);
fxeq_png_a = fxeq_png(:,a);
fxeq_png_a = fxeq_png_a(:);
a_leg_png = '';
for i = 1:length(a)
    if i ~= 1
        a_leg_png = a_leg_png + ', ';
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
    a_leg_png = a_leg_png + string(all_mutations(a(i)));
end

b = 4;
el_time_png_b = elapsed_time_png(:,b);
el_time_png_b = el_time_png_b(:);
fxeq_png_b = fxeq_png(:,b);
fxeq_png_b = fxeq_png_b(:);
b_leg_png = '';
for i = 1:length(b)
    if i ~= 1
        b_leg_png = b_leg_png + ', ';
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
    b_leg_png = b_leg_png + string(all_mutations(b(i)));
end

c = 1:10;
c = c(~ismember(c,a));
c = c(~ismember(c,b));
el_time_png_c = elapsed_time_png(:,c);
el_time_png_c = el_time_png_c(:);
fxeq_png_c = fxeq_png(:, c);
fxeq_png_c = fxeq_png_c(:);
c_leg_png = '';
for i = 1:length(c)
    if i ~= 1
        c_leg_png = c_leg_png + ', ';
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
    c_leg_png = c_leg_png + string(all_mutations(c(i)));
end



el_time_dyn_a = elapsed_time_dyn(:,a);
el_time_dyn_a = el_time_dyn_a(:);
fxeq_dyn_a = fxeq_dyn(:,a);
fxeq_dyn_a = fxeq_dyn_a(:);

el_time_dyn_b = elapsed_time_dyn(:,b);
el_time_dyn_b = el_time_dyn_b(:);
fxeq_dyn_b = fxeq_dyn(:,b);
fxeq_dyn_b = fxeq_dyn_b(:);

el_time_dyn_c = elapsed_time_dyn(:,c);
el_time_dyn_c = el_time_dyn_c(:);
fxeq_dyn_c = fxeq_dyn(:, c);
fxeq_dyn_c = fxeq_dyn_c(:);

%% 
f_prec = figure('units','normalized','outerposition',[0 0 1 1]);

col(1, :) = [0.8500, 0.3250, 0.0980];
col(2, :) = [0.4940, 0.1840, 0.5560];
col(3, :) = [0.9290, 0.6940, 0.1250];


subplot(2,2,1);
line1 = semilogy(el_time_png_a, fxeq_png_a,'diamond', 'Linewidth', 1, 'Markersize', 10, ...
    'DisplayName', a_leg_png, 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', col(1,:)); 
hold on
line2 = semilogy(el_time_png_b, fxeq_png_b, '.', 'MarkerSize', 30, ...
   'DisplayName', b_leg_png, 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', col(2,:));  
line3 = semilogy(el_time_png_c, fxeq_png_c, '+', 'Linewidth', 3, 'Markersize', 10, ...
    'DisplayName', c_leg_png, 'MarkerFaceColor', col(3,:), 'MarkerEdgeColor', col(3,:)); 
grid on
hold off
title('\textbf{PNG}', 'Interpreter','latex')
xlabel ('Elapsed time [sec]', 'Interpreter', 'Latex');
ylabel ('Accuracy - $||\textbf{f}(\textbf{x}_{png})||$', 'Interpreter', 'Latex');
set(gca, 'Fontsize', 20, 'TickLabelInterpreter','latex')


subplot(2,2,2);
line1 = semilogy(el_time_dyn_a, fxeq_dyn_a,'diamond', 'Linewidth', 1, 'Markersize', 10, ...
      'DisplayName', a_leg_png, 'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', col(1,:)); 
hold on
line2 = semilogy(el_time_dyn_b, fxeq_dyn_b, '.', 'MarkerSize', 30, ...
     'DisplayName', b_leg_png, 'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', col(2,:));
line3 = semilogy(el_time_dyn_c, fxeq_dyn_c,'+', 'Linewidth', 3, 'Markersize', 10, ...
    'DisplayName', c_leg_png, 'MarkerFaceColor', col(3,:), 'MarkerEdgeColor', col(3,:));
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

%% Save
saveas(f_prec, fullfile(folder_figures, 'scatter_png_dyn.png'))

