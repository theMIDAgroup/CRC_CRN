clc
clear
close all

%% 
% This code allows to reproduce the figures on the effect of the action of 
% the combined drugs Dabrafenib and Trametinib on the CRC-CRN affected by
% a gain of function mutation of KRAS one after the other
%
% 'In-silico modelling of the mitogen-activated kinase (MAPK) pathway in
% colorectal cancer: mutations and targeted therapy.'
% Sommariva et al. 2023, preprint at 
% https://www.biorxiv.org/content/10.1101/2023.04.18.537359v1.full.pdf
%
% The results of the simulations required to run this code are provided 
% in the folder './results_paper/drugs'
%%

set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex'); 

% addpath('./data')
%% Step 1. Define general parameters
% 1.1. Starting mutation 
mut_prot = 'Ras';
drug1 = 'DBF'; drug2 = 'TMT'; 

% 1.2. Files and folders
target_folder = '../data';
folder_results = './results_paper';
folder_figures = './figures_paper';
file_species_names = fullfile(target_folder, 'CRC_CRN_species_names.mat');
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');
file_ris_drug = fullfile(folder_results, 'drugs', sprintf('test_review_order_40.00_240.00.mat'));
file_ris_paper= fullfile(folder_results, 'drugs', sprintf('nlpc_DBF_TMT_on_mut_%s_40.00_240.00.mat', mut_prot));

%% Step 2. Load and store data model
% 2.1. Load and store
load(file_species_names, 'species_names');
load(file_ris_phys, 'nlpc_phys')
load(file_ris_drug, 'ris_drug_order')
load(file_ris_paper, 'ris_drug')

n_species = numel(species_names);
x_eq_phys = nlpc_phys(1).x;
ndyn = length(ris_drug_order.chosen_times);

% 2.2. Drug dosage
init_DBF = ris_drug_order.init_drug_1;
init_TMT=ris_drug_order.init_drug_2;
chosen_time=ris_drug_order.chosen_times;

%% Step 3. Load and store drug effects data + compute modified geometric mean G
x_eq_mut = ris_drug.x_eq_mut;
x_eq_paper = ris_drug.x_eq_combo(1:n_species);

delta_mut = (x_eq_mut -  x_eq_phys) ./ x_eq_phys;
[delta_mut_sort, idx_delta_mut_sort] = sort(delta_mut, 'descend');

delta = zeros(n_species, ndyn+1);
aux_delta = (x_eq_paper - x_eq_phys) ./ x_eq_phys;
delta(:, ndyn+1) = aux_delta(idx_delta_mut_sort);

x_eq_test = zeros(numel(x_eq_phys), ndyn);

field = cell(1, ndyn);
for i = 1:ndyn
    if i < ndyn
        field{i} = strcat('time', num2str(i));
    else
        field{i} = 'nlpc';
    end
    aux = ris_drug_order.x_eq_mut_drug_1_2.(field{i});
    x_eq_test(:,i) = aux(1:n_species);
    delta(:,i) = (x_eq_test(:,i)-x_eq_phys) ./ x_eq_phys;
    delta(:,i) = delta(idx_delta_mut_sort,i);
end

%% Step 4. Plot differences paper-test
colors = { [0.0290, 0.6940, 0.1250], [0.8500, 0.3250, 0.0980], [0, 0.4470, 0.7410], [0.9290, 0.0940, 0.1250], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560]};
markers = {'-', '-', '-', '-', '-', '-'};
linewidths = [4.5, 4, 3.5, 3, 2.5, 2];
f_error_combo = figure('units', 'normalized', 'outerposition', [0 0 1 0.55]);
aux_sp = subplot(1, 4, 1:3);    
for id = 1:ndyn+1
    if id < ndyn
    plot(1:n_species, delta(:, id), markers{id}, 'Linewidth', linewidths(id), ...
        'Displayname', sprintf('$t^*=%1.1f$ min', chosen_time(id)/60), ...
    "Color", colors{id})
    elseif id == ndyn
    plot(1:n_species, delta(:, id), markers{id}, 'Linewidth', linewidths(id), ...
        'Displayname', '$t^*=t_{eq}$', ...
    "Color", colors{id})
    else
    plot(1:n_species, delta(:, id), markers{id}, 'Linewidth', linewidths(id), ...
        'Displayname', 'Concurrent', ...
    "Color", colors{id})
    end
    hold on
end

botAxs = gca;
plot(1:n_species, delta_mut_sort, 'k', 'Linewidth', 3, ...
    'Displayname', 'GoF {\it KRAS}') 
my_symlog('y', -1)
set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 15)
xlabel('Proteins $i$', 'Fontsize', 20, 'Interpreter', 'Latex')
ylabel('$d_i$ (sorted)', 'Fontsize', 20, 'Interpreter', 'Latex')
lgd = legend('show');
set(lgd, 'Fontsize', 20, 'Location', 'East')
xlim([1, n_species])

hL = subplot(1, 4, 4);
poshL = get(hL, 'position');
axis(hL, 'off')
set(lgd,'position', [poshL(1), poshL(2)+0.2, poshL(3), poshL(4)*0.4], 'Interpreter', 'Latex')
set(gcf, 'PaperPositionMode', 'auto');

%% Step 5. Save
saveas(f_error_combo, fullfile(folder_figures, 'test_drug_order.jpg'))