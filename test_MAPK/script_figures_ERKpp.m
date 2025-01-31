clc
clear
close all

% This code allows to reproduce the figures on the effect of the action of 
% different drugs on the concentration on a CRC-CRN affected by a gain of function
% mutation of KRAS. In particular, it shows the dynamic of ppERK under the
% effect of
% 1. Dabrafenib
% 2. Trametinib
% 3. Their combination
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
% 1.1. Starting mutation 
mut_prot = 'Ras';
drug1 = 'DBF';
drug2 = 'TMT';
drug = strcat(drug1, '_', drug2);
init_drug_all1 = [50, 37.5, 25, 12.5];
init_drug_all2 = [200, 150, 100, 50];
n_d = numel(init_drug_all1);

% 1.2. Files and folders
folder_results = './results_paper'; %togliere old
folder_figures = './figures_paper';
folder_data = '../data';
file_mim = fullfile(folder_data, 'CRC_CRN_nodrug_complete.mat');
file_species_names = fullfile(folder_data, 'CRC_CRN_species_names.mat');
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');
file_ris_mut = fullfile(folder_results, 'mutations', 'nlpc_mut_KRAS_perc_0.0.mat');

%% Step 2. Load and store data model.
load(file_species_names, 'species_names');
load(file_ris_phys, 'nlpc_phys')
load(file_mim, 'CMIM');
load(file_ris_mut, 'nlpc_mut')

x_eq_phys = nlpc_phys(1).x;
x_eq_mut_null = nlpc_mut(1).x;
delta_mut_null = (x_eq_mut_null - x_eq_phys) ./ x_eq_phys;

n_species = numel(species_names);

%% Step 3. Maximum ERK concentration
ind_ERK = find(strcmp(CMIM.species.names, 'ERK'));
N = CMIM.matrix.Nl;
ind_cl_ERK = find(N(:,ind_ERK)==1);
ERK_tot = N(ind_cl_ERK,:) * x_eq_phys;

%% Step 4. Load and store drug data
% DBF
DBF_dynamics = cell(1, numel(init_drug_all1));
DBF_time_dynamics = cell(1,numel(init_drug_all1));
DBF_delta_mut_drug = zeros(n_species, numel(init_drug_all1));
DBF_delta_combo = zeros(n_species, numel(init_drug_all1));

for id = 1:numel(init_drug_all1)    
    file_ris_drug = fullfile(folder_results, 'drugs',  ...
        sprintf('dyn_%s_on_mut_%s_%2.2f.mat', drug1, mut_prot, init_drug_all1(id)));
    load(file_ris_drug, 'ris_drug')
    x_eq_mut_drug = ris_drug.x_eq_mut_drug(1:n_species);
    x_eq_combo = ris_drug.x_eq_combo(1:n_species);
    
    DBF_delta_mut_drug(:, id) = (x_eq_mut_drug - x_eq_phys) ./ x_eq_phys;  
    DBF_delta_combo(:, id) = (x_eq_combo - x_eq_phys) ./ x_eq_phys;
    
    DBF_dynamics{id} = ris_drug.x_t_mut_drug;
    DBF_time_dynamics{id} = ris_drug.time_mut_drug;
end
x_eq_mut = ris_drug.x_eq_mut(1:n_species);
DBF_delta_mut = (x_eq_mut - x_eq_phys) ./ x_eq_phys;

% TMT
TMT_dynamics = cell(1, numel(init_drug_all2));
TMT_time_dynamics = cell(1, numel(init_drug_all2));
TMT_delta_mut_drug = zeros(n_species, numel(init_drug_all2));
TMT_delta_combo = zeros(n_species, numel(init_drug_all2));

for id = 1:numel(init_drug_all2)    
    file_ris_drug = fullfile(folder_results, 'drugs', ...
        sprintf('dyn_%s_on_mut_%s_%2.2f.mat', drug2, ...
        mut_prot, init_drug_all2(id)));
    load(file_ris_drug, 'ris_drug')
    x_eq_mut_drug = ris_drug.x_eq_mut_drug(1:n_species);
    x_eq_combo = ris_drug.x_eq_combo(1:n_species);
    
    TMT_delta_mut_drug(:, id) = (x_eq_mut_drug - x_eq_phys) ./ x_eq_phys;  
    TMT_delta_combo(:, id) = (x_eq_combo - x_eq_phys) ./ x_eq_phys;
    
    TMT_dynamics{id} = ris_drug.x_t_mut_drug;
    TMT_time_dynamics{id} = ris_drug.time_mut_drug;
end

x_eq_mut = ris_drug.x_eq_mut(1:n_species);
TMT_delta_mut = (x_eq_mut - x_eq_phys) ./ x_eq_phys;

% Combination
DBF_TMT_dynamics = cell(numel(init_drug_all1), numel(init_drug_all2));
DBF_TMT_time_dynamics = cell(numel(init_drug_all1), 1);
DBF_TMT_delta_mut_drug = zeros(n_species, numel(init_drug_all1), numel(init_drug_all2));
DBF_TMT_delta_combo = zeros(n_species, numel(init_drug_all1), numel(init_drug_all2));

for i=1:4
    file_ris_drug = fullfile(folder_results, 'drugs',...
        sprintf('dyn_%s_%s_on_mut_%s_%2.2f_%2.2f.mat', ...
        drug1,  drug2, mut_prot, init_drug_all1(i), init_drug_all2(i)));
    load(file_ris_drug, 'ris_drug')
    DBF_TMT_dynamics{i} = ris_drug.x_t_mut_drug;
    DBF_TMT_time_dynamics{i} = ris_drug.time_mut_drug;
end

x_eq_mut = ris_drug.x_eq_mut(1:n_species);
DBF_TMT_delta_mut = (x_eq_mut - x_eq_phys) ./ x_eq_phys;

%% P1. Dynamics
set_p1 = {'ERKPP'};
set_p1_name = {'$\frac{\mathrm{p-p-ERK}}{\mathrm{ERK_{tot}}}$'};

[~, idx_p1] = ismember(set_p1, species_names);

is = 1;
xlim_left = 10^-2;
xlim_right = 100;

f_dyn_p1 = figure('units','normalized','outerposition',[0 0 0.6 0.9]);

% DBF
subplot(3, 2, 1);
aux_is = idx_p1(is);
for id = 1:n_d
    aux_time = DBF_time_dynamics{id};
    aux_time_hours = aux_time/60;
    semilogx(aux_time_hours(1:end), DBF_dynamics{id}(aux_is, 1:end)/ERK_tot, ...
        'Linewidth', 3, ...
        'Displayname', sprintf('$c_D = %1.1f$ nM', init_drug_all1(id)))
    hold on
end

semilogx(aux_time_hours(1:end), x_eq_phys(aux_is)*ones(size(aux_time))/ERK_tot, ...
   'k--', 'Linewidth', 3, 'Displayname', 'phys')
grid on
ylabel(set_p1_name{is}, 'Fontsize', 20, 'Interpreter', 'Latex')
lgd = legend('show');
set(gca, 'TickLabelInterpreter','latex', 'Fontsize', 20)
set(lgd, 'Fontsize', 20, 'Interpreter', 'Latex')
xlim([xlim_left, xlim_right])
xticks([0.05, 1, 13, 60])
yticks([0.25, 0.50])
ttl = title('(A) DBF');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';

hL = subplot(3, 2, 2);
poshL = get(hL, 'position');
axis(hL, 'off')
set(lgd,'position', [poshL(1), poshL(2), poshL(3)*0.6, poshL(4)], 'Interpreter', 'Latex')

% TMT
subplot(3, 2, 3);
aux_is = idx_p1(is);
for id = 1:n_d
    aux_time = TMT_time_dynamics{id};
    aux_time_hours = aux_time/60;
    semilogx(aux_time_hours(1:end), TMT_dynamics{id}(aux_is, 1:end)/ERK_tot, ...
        'Linewidth', 3, ...
        'Displayname', sprintf('$c_T = %1.0f$ nM', init_drug_all2(id)))
    hold on
end

semilogx(aux_time_hours(1:end), x_eq_phys(aux_is)*ones(size(aux_time))/ERK_tot, ...
   'k--', 'Linewidth', 3, 'Displayname', 'phys')
grid on
ylabel(set_p1_name{is}, 'Fontsize', 20, 'Interpreter', 'Latex')
lgd = legend('show');
set(gca, 'TickLabelInterpreter','latex', 'Fontsize', 20)
set(lgd, 'Fontsize', 20, 'Interpreter', 'Latex')
xlim([xlim_left, xlim_right])
xticks([0.05, 1, 13, 60])
yticks([0.25, 0.50])
ttl = title('(B) TMT');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';

hL = subplot(3, 2, 4);
poshL = get(hL, 'position');
axis(hL, 'off')
set(lgd,'position', [poshL(1), poshL(2), poshL(3)*0.6, poshL(4)], 'Interpreter', 'Latex')

% Combination
subplot(3, 2, 5);
aux_is = idx_p1(is);
for i=1:4
    if i==4
        dyn = DBF_TMT_dynamics{i}(aux_is,:);
        aux_time = DBF_TMT_time_dynamics{i};
        aux_time_hours = aux_time/60;
        semilogx(aux_time_hours(1:length(dyn)), dyn/ERK_tot, ...
            'Linewidth', 3, ...
            'Displayname', sprintf('$c_D = %1.1f$ nM, $c_T = %1.0f$ nM', init_drug_all1(i), init_drug_all2(i)))
    else
        aux_time = DBF_TMT_time_dynamics{i};
        aux_time_hours = aux_time/60;
        semilogx(aux_time_hours(1:end), DBF_TMT_dynamics{i}(aux_is, :)/ERK_tot, ...
            'Linewidth', 3, ...
            'Displayname', sprintf('$c_D = %1.1f$ nM, $c_T = %1.0f$ nM', init_drug_all1(i), init_drug_all2(i)))
    end 
           hold on
end
last_ind = 105;
semilogx(aux_time_hours(1:last_ind), x_eq_phys(aux_is)*ones(1,last_ind)/ERK_tot, ...
   'k--', 'Linewidth', 3, 'Displayname', 'phys')
grid on
ylabel(set_p1_name{is}, 'Fontsize', 20, 'Interpreter', 'Latex')
lgd = legend('show');
ttl = title('(C) Combination');
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0;
ttl.HorizontalAlignment = 'left';
set(gca, 'TickLabelInterpreter','latex', 'Fontsize', 20)
set(lgd, 'Fontsize', 20)
xlim([xlim_left, xlim_right])
xlabel('Time (min)', 'Fontsize', 20, 'Interpreter', 'Latex')
xticks([0.05, 1, 13, 60])
yticks([0.25, 0.50])

hL = subplot(3, 2, 6);
poshL = get(hL, 'position');
axis(hL, 'off')
set(lgd,'position', [poshL(1), poshL(2), poshL(3)*0.6, poshL(4)], 'Interpreter', 'Latex')
set(gcf, 'PaperPositionMode', 'auto');

%% Save
% saveas(f_dyn_p1, fullfile(folder_figures, 'dyn_ppERK.jpg'))

