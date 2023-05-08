clc
clear
close all

%% 
% This code allows to reproduce the figures on the effect of the action of 
% the combined drugs Dabrafenib and Trametinib on the CRC-CRN affected by
% a gain of function mutation of KRAS presented in the paper
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

%% Step 1. Define general parameters
% 1.1. Starting mutation 
mut_prot = 'Ras';
drug1 = 'DBF'; drug2 = 'TMT'; drug = strcat(drug1, '_', drug2);

% 1.2. Files and folders
target_folder = '../data/ci_servono';
folder_results = './results_paper';
folder_figures = './figures_paper';
file_species_names = fullfile(target_folder, 'CRC_CRN_species_names.mat');
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');
file_ris_mut = fullfile(folder_results, sprintf('nlpc_mut_%s.mat', mut_prot));

%% Step 2. Load and store data model
% 2.1. Load and store
load(file_species_names, 'species_names');
load(file_ris_phys, 'nlpc_phys')
load(file_ris_mut, 'nlpc_mut')

n_species = numel(species_names);
x_eq_phys = nlpc_phys(1).x;

% 2.2. Drug dosage
init_DBF_all1 = linspace(0, 100, 11); init_DBF_all1 = init_DBF_all1(2:end);
init_TMT_all1= linspace(0, 2000, 11); init_TMT_all1 = init_TMT_all1(2:end);
n_d1 = numel(init_DBF_all1); n_d2 = numel(init_TMT_all1);

init_DBF_all2 = [40, 30, 30, 20];
init_TMT_all2 = [240, 1000, 600, 600];
n_p1 = numel(init_DBF_all2); n_p2 = numel(init_TMT_all2);
init_drug_all2_vec = [init_DBF_all2' init_TMT_all2'];

%% Step 3. Load and store drug effects data + compute modified geometric mean G
% 3.1. Load and store
delta_combo = zeros(n_species,  numel(init_DBF_all1), numel(init_TMT_all1));
geom_delta_combo = zeros(numel(init_DBF_all1), numel(init_TMT_all1));

mu = 1e-6;

for id1 = 1:n_d1  
    for id2 = 1:n_d2
            file_ris_drug = fullfile(folder_results, 'drugs', ...
                sprintf('nlpc_%s_on_mut_%s_%2.2f_%2.2f.mat', drug, ...
                mut_prot, init_DBF_all1(id1), init_TMT_all1(id2)));
            load(file_ris_drug, 'ris_drug')
            x_eq_combo = ris_drug.x_eq_combo(1:end-5);
            delta_combo(:, id1, id2) = (x_eq_combo - x_eq_phys) ./ x_eq_phys;
            geom_delta_combo(id1, id2) = geomean(abs(delta_combo(:, id1, id2))+mu) - mu;
    end
end

%% Step 4. Plot - modified geometric mean G
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560]};
f_error = figure('units','normalized','outerposition',[0 0 1 0.6]);

alignPlot = subplot(1,3,1);
topAxs = gca;
photoAxsRatio = get(topAxs,'PlotBoxAspectRatio');
imagesc(init_TMT_all1, init_DBF_all1, geom_delta_combo);
hold on
for i=1:4
    plot(init_TMT_all2(i), init_DBF_all2(i), '*', 'MarkerSize',  18, 'LineWidth', 2.5); 
    hold on
end 
colormap(flipud(parula));
colorbar; 
axis('square');
xlabel(sprintf('Concentration of %s [nM]', drug2), 'Fontsize', 20, 'Interpreter', 'Latex');
c = colorbar; pos = get(c,'Position'); 
t = sprintf('$G(\\mathbf{d})$', log10(mu), log10(mu));
c.Label.String = t; 
c.Label.Rotation = 270; c.Label.FontSize = 20; c.Label.Interpreter = 'Latex';
c.Label.Position = [pos(1)+3.5 pos(2)-0.14];
ylabel(sprintf('Concentration of %s [nM]', drug1), 'Fontsize', 20, 'Interpreter', 'Latex');
xvec = 4:4:20; xvec = 100 * xvec;
yvec = 2:2:10; yvec = 10 * yvec;
title('(A)', 'HorizontalAlignment', 'left', 'position', [177 107], 'Fontsize', 21)
set(gca,'YDir','normal','XTick',xvec,'YTick',yvec) 

%% Step 3. Load and store drug data
dynamics = cell(size(init_drug_all2_vec, 1));
time_dynamics = cell(size(init_drug_all2_vec, 1));
delta_combo = zeros(n_species, size(init_drug_all2_vec, 1));

for id = 1:length(init_DBF_all2) 
    file_ris_drug = fullfile(folder_results, 'drugs', ...
        sprintf('nlpc_%s_on_mut_%s_%2.2f_%2.2f.mat', drug, ...
        mut_prot, init_DBF_all2(id), init_TMT_all2(id)));
    load(file_ris_drug, 'ris_drug')
    x_eq_combo = ris_drug.x_eq_combo(1:n_species);
    delta_combo(:, id) = (x_eq_combo - x_eq_phys) ./ x_eq_phys;
end
x_eq_mut = ris_drug.x_eq_mut;
delta_mut = (x_eq_mut - x_eq_phys) ./ x_eq_phys;

%% Step 4. Profile of all the concentrations
markers = {'s--', 'd--', '*--', 'o--'};
    
[delta_mut_sort, idx_delta_mut_sort] = sort(delta_mut, ...
    'descend');
aux = squeeze(delta_combo);
delta_combo_sort = aux(idx_delta_mut_sort, :);
aux_sp = subplot(1,3,2:3);
for id = 1:n_p1
    plot(1:n_species, delta_combo_sort(:, id), markers{id}', 'Linewidth', 2, ...
        'Markersize', 8,'Displayname', ...
        sprintf('$c_{D} = %d$ nM, $c_{T} = %d$ nM', init_DBF_all2(id), ...
        init_TMT_all2(id)))
    hold on
end
botAxs = gca;
plot(1:n_species, delta_mut_sort, 'k', 'Linewidth', 3, ...
    'Displayname', 'GoF KRAS') 
my_symlog('y')
title('(B)', 'HorizontalAlignment', 'left', 'position', [5 5.1], 'Fontsize', 20)
xlabel('Proteins $i$', 'Fontsize', 20, 'Interpreter', 'Latex')
ylabel('$d_i$ (sorted)', 'Fontsize', 20, 'Interpreter', 'Latex')
set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 20)
lgd = legend('show');
set(lgd, 'Fontsize', 18, 'Location', 'North')
xlim([1, n_species])
pos1 = get(aux_sp, 'Position'); new_pos1 = pos1 +[0.04 0.023 0 -0.023];
set(aux_sp, 'Position', new_pos1) ;

%% Step 6. Save
saveas(f_error, fullfile(folder_figures, 'full_performance_DBF_TMT.jpg'))
