clc
clear
close all

% This code allows to reproduce the figures on the effect of the action of 
% optimal drug on the concentration on a CRC-CRN affected by a gain of function
% mutation of KRAS. In particular, it shows the dynamic of ppERK under the
% effect of Dabrafenib until time t* and under of the combined effect with
% Trametinib after t*.
% t* is determined as the time where the concentration of ERKPP under
% mutation + DBF is similar to the physiological case
%
% 'In-silico modelling of the mitogen-activated kinase (MAPK) pathway in 
% colorectal cancer: mutations and targeted therapy.'
% Sommariva et al. 2023, preprint at 
% https://www.biorxiv.org/content/10.1101/2023.04.18.537359v1.full.pdf
% The results of the simulations required to run this code are provided 
% in the folder './results_paper'

set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex'); 

%% Step 1. Define general parameters
% 1.1. Mutation + Drugs
mut_prot = 'Ras';
drug1 = 'DBF';
drug2 = 'TMT';
drug = strcat(drug1, '_', drug2);

% 1.2. Files and folders
folder_results = './results_paper';
folder_figures = './figures_paper';
folder_drug = './results_paper/drugs';
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');
file_ris_drug = fullfile(folder_drug, 'test_review_ERKPP_TMT_DBF_deg.mat');

%% Step 2. Load data
load(file_ris_phys, 'nlpc_phys')
load(file_ris_drug, 'ris_drug_test')

%% Step 3. Store data
x_eq_phys = nlpc_phys(1).x;

init_drug_1 = ris_drug_test.init_drug_1;
init_drug_2 = ris_drug_test.init_drug_2;

ind_ppERK = ris_drug_test.indexERKPP;
ERK_tot = ris_drug_test.ERK_tot;

x_eq_phys_ERK = x_eq_phys(ind_ppERK);

t_star = ris_drug_test.t_min;  
time_DBF_deg = ris_drug_test.time_mut_drug;
time_DBF = time_DBF_deg(1:t_star-1);
time_TMT_add = ris_drug_test.t0_0.time_mut_drug_1_2;
time_combo_1_2_0 = ris_drug_test.time_mut_drug_1_2_combo;

x_t_DBF_deg = ris_drug_test.x_t_mut_drug;
x_t_DBF = x_t_DBF_deg(:,1:t_star-1);
x_t_TMT = ris_drug_test.t0_0.x_t_mut_drug_1_2;
x_t_TMT_tstar = ris_drug_test.t0_tstar.x_t_mut_drug_1_2;

x_t_combo_ERK = [x_t_DBF(ind_ppERK,:), x_t_TMT(ind_ppERK,:)]';

x_t_combo_1_2_0 = ris_drug_test.x_t_mut_drug_1_2_combo;

%% Step 4. Figure. Comparison: dynamic of ppERK under the effect of DBF (with degradation) or DBF+TMT (TMT added at time t*)

col = zeros(3,3);
col(1, :) = [0.8500, 0.3250, 0.0980]; 
col(2, :) = [0,128,0] / 256;
col(3, :) = [222, 49, 99] / 256;

set_p1_name = {'$\frac{\mathrm{p-p-ERK}}{\mathrm{ERK_{tot}}}$'};

aux_time = [time_DBF; time_TMT_add + time_DBF_deg(t_star)];
aux_time_hours_deg = time_DBF_deg / 60;
aux_time_hours = aux_time / 60;

xlim_left = aux_time_hours(1);
xlim_right = 3*1e5;

fig_erk = figure('units','normalized','outerposition', [0.1 0.1 0.7 0.6]);

subplot(1, 2, 1)
semilogx(aux_time_hours, x_t_combo_ERK / ERK_tot, ...
        'Linewidth', 3, 'Markersize', 10, ...
    'MarkerFaceColor', col(1,:), 'MarkerEdgeColor', col(1,:),...
    'Displayname', strcat('$c_D$ = ', num2str(init_drug_1), 'nM and $c_T$ = ', num2str(init_drug_2), 'nM'))
hold on
semilogx(aux_time_hours_deg, x_t_DBF_deg(ind_ppERK,:) / ERK_tot, ...
          'Linewidth', 3, 'Linestyle', '--',...
      'Markersize', 10, ...
    'MarkerFaceColor', col(2,:), 'MarkerEdgeColor', col(2,:),...
    'Displayname', strcat('$c_D$ = ', num2str(init_drug_1), ' nM'))
hold on
semilogx(aux_time_hours, x_eq_phys_ERK * ones(numel(aux_time_hours),1)/ERK_tot, ...
        'k--', 'Linewidth', 3,...
        'Displayname', 'phys')

lgd = legend('show');
my_symlog('y')
xlabel('Time (min)', 'Fontsize', 20, 'Interpreter', 'Latex')
ylabel(set_p1_name, 'Fontsize', 20, 'Interpreter', 'Latex')
grid on

set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 20)
xlim([xlim_left xlim_right])   
xticks([0.05, 1, aux_time_hours_deg(t_star), 4800])
xticklabels({'0.05', '1', 't*', '4800'})
yticks([0.8, 1.6])
yticklabels({'0.25', '0.5'})
hL = subplot(1, 2, 2);
poshL = get(hL, 'position');
axis(hL, 'off')
set(lgd,'position', [poshL(1), poshL(2) * 0.1 + 0.5, poshL(3), poshL(4) * 0.2], 'Interpreter', 'Latex')
set(gcf, 'PaperPositionMode', 'auto');
hold off

%% Step 5. Figure. Comparison: dynamic of ppERK under the effect of DBF (with degradation) or DBF+TMT (TMT added at time 0)

aux_time_TMT_0 = time_combo_1_2_0;
aux_time_hours_TMT_0 = aux_time_TMT_0/60;

xlim_left = aux_time_hours(1);
xlim_right = 3*1e5;

fig_erk_at_0 = figure('units','normalized','outerposition',[0 0 0.65 0.35]);

subplot(1, 2, 1)
semilogx(aux_time_hours_deg, x_t_DBF_deg(ind_ppERK,:)/ERK_tot, ...
          'Linewidth', 4.5 , ...
          'Displayname', 'no TMT')
hold on
semilogx(aux_time_hours_TMT_0, x_t_combo_1_2_0(ind_ppERK, :)/ERK_tot, ...
        'Linewidth', 4, 'Color', col(2,:),...
    'Displayname', sprintf('TMT at $t^*= %1.1f$ min', 0))
semilogx(aux_time_hours, x_t_combo_ERK(:)/ERK_tot, ...
          'Linewidth', 2.5,  'Color', col(3,:),...
    'Displayname',sprintf('TMT at $t^*= %1.1f$ min', aux_time_hours(t_star)))    
semilogx(aux_time_hours, x_eq_phys_ERK*ones(numel(aux_time_hours),1)/ERK_tot, ...
        'k--', 'Linewidth', 3,...
        'Displayname', 'phys')

disp(x_t_DBF_deg(ind_ppERK,end) /  ERK_tot)
disp(x_t_combo_ERK(end) / ERK_tot)
disp(x_t_combo_1_2_0(ind_ppERK, end)/ERK_tot)
disp(x_eq_phys_ERK / ERK_tot)
    
lgd = legend('show');
my_symlog('y')
xlabel('Time (min)', 'Fontsize', 20, 'Interpreter', 'Latex')
ylabel(set_p1_name, 'Fontsize', 20, 'Interpreter', 'Latex')
grid on

set(gca, 'TickLabelInterpreter','latex', 'Fontsize', 20)
xlim([xlim_left xlim_right])   
xticks([0.05, 1, 4800])
xticklabels({'0.05', '1',  '4800'})
yticks([0.75 1.65])
yticklabels({'0.25', '0.5'})
hL = subplot(1, 2, 2);
poshL = get(hL, 'position');
axis(hL, 'off')
set(lgd,'position', [poshL(1), poshL(2)*0.1+0.5, poshL(3), poshL(4)*0.2], 'Interpreter', 'Latex')
set(gcf, 'PaperPositionMode', 'auto');

hold off

%% Step 7. Save
saveas(fig_erk, fullfile(folder_figures, 'test_ERKPP_TMT_DBF_.jpg'))
saveas(fig_erk_at_0, fullfile(folder_figures, 'test_ERKPP_add_TMT_DBF_at_0.jpg'))

