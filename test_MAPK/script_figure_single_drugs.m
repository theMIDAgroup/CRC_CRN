clc
clear
close all

%% 
% This code allows to reproduce the figures on the effect of the action of 
% the single drugs Dabrafenib and Trametinib on the CRC-CRN affected by
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
% 1.1. Mutation and drug
mut_prot = 'Ras'; prot_name = 'KRAS'; perc = 0;
drugs = {'DBF', 'TMT'};

% 1.2. Files and folders
target_folder = '../data';
folder_results = './results_paper';
folder_figures = './figures_paper';
file_species_names = fullfile(target_folder, 'CRC_CRN_species_names.mat');
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');
file_ris_mut = fullfile(folder_results, ...
    'mutations', sprintf('nlpc_mut_%s_perc_%1.1f.mat', prot_name, perc));

%% Step 2. Load and store data model.
% 2.1. Load and store
load(file_species_names, 'species_names');
load(file_ris_phys, 'nlpc_phys')
load(file_ris_mut, 'nlpc_mut')

n_species = numel(species_names);
x_eq_phys = nlpc_phys(1).x;

for j=1:numel(drugs)
    drug = drugs{j};
    % 2.2. Drug dosage
    if strcmp(drug, 'DBF')
        n_new_species = 2;
        init_drug_all1 = linspace(0, 100, 11); init_drug_all1 = init_drug_all1(2:end);
        init_drug_all2 = [75, 55, 43.5, 25];
        conc_name = 'c_D';
    elseif strcmp(drug, 'TMT')
        n_new_species = 3;
        init_drug_all1 = linspace(0, 2000, 11); init_drug_all1 = init_drug_all1(2:end);
        init_drug_all2 = [1400, 1080, 700, 300];
        conc_name = 'c_T';
    else
    end
    init_drug_all = [init_drug_all1 init_drug_all2];
    n_d = numel(init_drug_all2);

    %% Step 3. Load and store drug effects data + compute modified geometric mean G
    % 3.1. Load and store
    delta_combo = zeros(n_species,  numel(init_drug_all));

    for id = 1:numel(init_drug_all)    
        file_ris_drug = fullfile(folder_results, 'drugs', ...
            sprintf('nlpc_%s_on_mut_%s_%2.2f.mat', drug, mut_prot, init_drug_all(id)));
        load(file_ris_drug, 'ris_drug')
        x_eq_combo = ris_drug.x_eq_combo(1:end-n_new_species);
        delta_combo(:, id) = (x_eq_combo - x_eq_phys) ./ x_eq_phys;
    end

    x_eq_mut = ris_drug.x_eq_mut;
    delta_mut = (x_eq_mut - x_eq_phys) ./ x_eq_phys;

    % 3.2. Compute modified geometrical mean G
    mu = 1e-6;
    geom_delta_combo = geomean(abs(delta_combo)+mu, 1)-mu;
    geom_delta_mut = geomean(abs(delta_mut)+mu, 1)-mu;

    %% Step 4. Plot - modified geometric mean G
    f_error = figure('units','normalized','outerposition',[0 0 1.5 0.6]);
    subplot(1,3,1);
    plot([0, init_drug_all1], [geom_delta_mut, geom_delta_combo(1:numel(init_drug_all1))], ...
        '.', 'LineWidth', 2, 'Markersize', 23, 'Color', [0.4660, 0.6740, 0.1880]);
    hold on
    colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560]};
    for id = 1:numel(init_drug_all2)
        plot(init_drug_all2(id), geom_delta_combo(numel(init_drug_all1)+id), '*', 'LineWidth', 2, 'Markersize', 10, 'Color', colors{id});
    end
    title('(A)', 'HorizontalAlignment', 'left', 'position', [1.5 0.2025])
    xlabel (sprintf('Concentration of %s [nM]', drug), 'Fontsize', 20, 'Interpreter', 'Latex')
    t = '$G(\mathbf{d})$';
    ylabel (t, 'Fontsize', 20, 'Interpreter', 'Latex')
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 20)
    hold on

    %% Step 5. Plot - profile of all the concentrations
    markers = {'s--', 'd--', '*--', 'o--'};

    [delta_mut_sort, idx_delta_mut_sort] = sort(delta_mut(:), 'descend');
    aux = squeeze(delta_combo);
    delta_combo_sort = aux(idx_delta_mut_sort, :);
    aux_sp = subplot(1,3,2:3);
    for id = 1:n_d
        init_drug = init_drug_all2(id);
        if strcmp(drug, 'DBF')
        plot(1:n_species, delta_combo_sort(:, numel(init_drug_all1)+id), markers{id}', 'Linewidth', 2, ...
            'Markersize', 8, 'Displayname', sprintf('$%s = %2.1f$ nM', conc_name, init_drug))
        elseif strcmp(drug, 'TMT')
        plot(1:n_species, delta_combo_sort(:, numel(init_drug_all1)+id), markers{id}', 'Linewidth', 2, ...
            'Markersize', 8, 'Displayname', sprintf('$%s = %2.0f$ nM', conc_name, init_drug))
        end
        hold on
    end
    plot(1:n_species, delta_mut_sort, 'k', 'Linewidth', 3, 'Displayname', 'GoF {\it KRAS}') 
    title('(B)', 'HorizontalAlignment', 'left', 'position', [1.5 5.1])
    my_symlog('y')
    xlabel('Proteins $i$', 'Fontsize', 20, 'Interpreter', 'Latex')
    ylabel('$d_i$ (sorted)', 'Fontsize', 20, 'Interpreter', 'Latex')
    set(gca, 'TickLabelInterpreter', 'latex', 'Fontsize', 20)
    lgd = legend('show');
    set(lgd, 'Fontsize', 18, 'Location', 'North')
    xlim([1, n_species])
    pos1 = get(aux_sp, 'Position');
    new_pos1 = pos1 +[0 0.023 0 -0.023];
    set(aux_sp, 'Position', new_pos1) ;

    %% Step 6. Save
     saveas(f_error, fullfile(folder_figures, sprintf('full_performance_%s.jpg', drug)))
end