clc
clear
close all

%% 
% This code allows to reproduce the figures on the effect of the action of 
% a Dabrafenib on the CRC-CRN affected by a gain of function mutation of KRAS
% presented in the paper
%
% 'Computational quantification of global effects induced by mutations and 
%   drugs in signaling networks of colorectal cancer cells.'
% Sommariva et al. 2020, preprint at 
% https://www.biorxiv.org/content/biorxiv/early/2020/12/31/2020.12.30.424842.full.pdf
%
% The results of the simulations required to run this code are provided 
% in the folder './results'
%%

set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex'); 

%% Step 1. Define general parameters
% 1.1. Starting mutation 
mut_prot = 'Ras';
k1_drug = 0.5;
perc_all = [0];  % Level of the mutation
perc_clt = [0.25, 0.5, 0.75, 1.0]; % Drug initial concentration (defined as
                  % function of the constant aggregate concentration within
                  % Braf conservation law
n_p = numel(perc_all); n_d = numel(perc_clt);

% 1.2. File and folders
folder_results = './results';
folder_figures = './figures';
file_species_names = fullfile(folder_results, 'CRC_CRN_species_names.mat');
file_ris_phys = fullfile(folder_results, 'results_physiological.mat');
file_ris_mut = fullfile(folder_results, ...
                sprintf('results_mutation_%s_perc_%1.1f.mat', mut_prot, 0));


%% Step 2. Load and store data model.
load(file_species_names, 'species_names');
load(file_ris_phys, 'ris_phys')

n_species = numel(species_names);
x_eq_phys = ris_phys.x_eq;

protein_dynamics = {'Ras', 'Raf', 'Raf*', 'ERK', 'MEK'};
[~, idx_protein_dynamics] = ismember(protein_dynamics, species_names);

load(file_ris_mut, 'ris_mutated')
x_eq_mut_null = ris_mutated.x_xemut_eq;
delta_mut_null = (x_eq_mut_null - x_eq_phys) ./ x_eq_phys; 

%% Step 3. Load and store drug data
dynamics = cell(numel(perc_all), numel(perc_clt));
time_dynamics = cell(numel(perc_all), numel(perc_clt));
delta_mut = zeros(n_species, numel(perc_all));
delta_mut_drug = zeros(n_species, numel(perc_all), numel(perc_clt));
delta_combo = zeros(n_species, numel(perc_all), numel(perc_clt));

for ip = 1:numel(perc_all)
perc = perc_all(ip);
for id = 1:numel(perc_clt)    
    file_ris_drug = fullfile(folder_results, ...
        sprintf('results_drug_on_mut_%s_alpha_%2.2f.mat', ...
        mut_prot, perc_clt(id)));
    load(file_ris_drug, 'ris_drug')

    x_eq_mut_drug = ris_drug.x_eq_mut_drug;
    x_eq_combo = ris_drug.x_eq_combo;
    
    delta_mut_drug(:, ip, id) = (x_eq_mut_drug - x_eq_phys) ./ x_eq_phys;  
    delta_combo(:, ip, id) = (x_eq_combo - x_eq_phys) ./ x_eq_phys;
    
    dynamics{ip, id} = ris_drug.x_t_mut_drug;
    time_dynamics{ip, id} = ris_drug.time_mut_drug;
end
    x_eq_mut = ris_drug.x_eq_mut;
    delta_mut(:, ip) = (x_eq_mut - x_eq_phys) ./ x_eq_phys;
end

%% P1. Dynamics
set_p1 = {'Ras'};
set_p1_name = {'k-Ras'};
set_p2 = {'Raf', 'Raf*', 'MEK', 'MEKP', 'ERK', 'ERKP'};
set_p2_name = {'Raf', 'p-Raf', 'MEK', 'p-MEK', 'ERK', 'p-ERK'};

[~, idx_p1] = ismember(set_p1, species_names);
[~, idx_p2] = ismember(set_p2, species_names);

f_dyn_p1 = figure('units','normalized','outerposition',[0 0 0.6 0.35]);
for is = 1:numel(idx_p1)
    subplot(numel(idx_p1), 2, 2*is-1)
    aux_is = idx_p1(is);
    
    for id = 1:n_d
        aux_time = time_dynamics{ip, id} /3600;
        semilogx(aux_time, dynamics{ip, id}(aux_is, :), ...
            'Linewidth', 3, ...
            'Displayname', sprintf('$\\alpha = %1.2f$',perc_clt(id)))
        hold on
    end
    semilogx(aux_time, x_eq_phys(aux_is)*ones(size(time_dynamics{ip, id})), ...
       'k--', 'Linewidth', 3, 'Displayname', 'phys')
    grid on
    ylabel(sprintf('%s [nM]', set_p1_name{is}), 'Fontsize', 20)
    lgd = legend('show');
    set(gca, 'Fontsize', 20)
    set(lgd, 'Fontsize', 25)
    xlim([aux_time(2), aux_time(end)])

    xticks([10^(-5) 10^0 10^4])
    xticklabels({'10^{-5}','10^0','10^4'})

    hL = subplot(numel(idx_p1), 2, 2*is);
    poshL = get(hL, 'position');
    poshL(3) = 0.7 * poshL(3);
    poshL(4) = 0.9 * poshL(4);
    
    set(lgd, 'position', poshL);
    axis(hL, 'off')
    
   
end
set(gcf, 'PaperPositionMode', 'auto');
saveas(f_dyn_p1, fullfile(folder_figures, 'Effect_drug_Ras.png'))


f_dyn_p2 = figure('units','normalized','outerposition',[0 0 0.6 1]);
for is = 1:numel(idx_p2)
    subplot(numel(idx_p2)/2, 2, is)
    aux_is = idx_p2(is);
    
    for id = 1:n_d
        aux_time = time_dynamics{ip, id} /3600;
        semilogx(aux_time, dynamics{ip, id}(aux_is, :), ...
            'Linewidth', 3)
        hold on
    end
    semilogx(aux_time, x_eq_phys(aux_is)*ones(size(time_dynamics{ip, id})), ...
       'k--', 'Linewidth', 3, 'Displayname', 'phys')
    grid on
    ylabel(sprintf('%s [nM]', set_p2_name{is}), 'Fontsize', 20)
    set(gca, 'Fontsize', 20)
    xlim([aux_time(2), aux_time(end)])
    if is == numel(idx_p2) || is == numel(idx_p2)  - 1
        xlabel('Time [s]', 'FontSize', 20)
    end
    xticks([10^(-5) 10^0 10^4])
    xticklabels({'10^{-5}','10^0','10^4'})

end
set(gcf, 'PaperPositionMode', 'auto');
saveas(f_dyn_p2, fullfile(folder_figures, 'Effect_drug_MAPK_cascade.png'))

%% P2. Profile of all the concentrations
markers = {'s--', 'd--', '*--', 'o--'};
n_species_no_drug = n_species - 2;
    
[delta_mut_sort, idx_delta_mut_sort] = sort(delta_mut(1:n_species_no_drug, ip), ...
    'descend');
aux = squeeze(delta_combo(:, ip, :));
delta_combo_sort = aux(idx_delta_mut_sort, :);
f_eff_pub = figure('units','normalized','outerposition',[0 0 1 0.5]);
for id = 1:n_d
plot(1:n_species_no_drug, delta_combo_sort(:, id), markers{id}',...
    'Linewidth', 2, 'Markersize', 8, ...
    'Displayname', sprintf('$\\alpha = %1.2f$', perc_clt(id)))
hold on
end
plot(1:n_species_no_drug, delta_mut_sort, 'k', 'Linewidth', 3, ...
    'Displayname', 'GoF Ras') 
my_symlog('y')
xlabel('Proteins $i$', 'Fontsize', 20, 'Interpreter', 'Latex')
ylabel('$\delta_i$ (sorted)', 'Fontsize', 20, 'Interpreter', 'Latex')
set(gca, 'TickLabelInterpreter','latex', 'Fontsize', 20)
lgd = legend('show');
set(lgd, 'Fontsize', 20, 'Location', 'EastOutside')
xlim([1, n_species_no_drug])

saveas(f_eff_pub, fullfile(folder_figures, 'Effect_drug_whole_profile.tiff'))

%% Related table for the best results (supplementary material)
file_ris_drug = fullfile(folder_results, ...
        sprintf('results_drug_on_mut_%s_alpha_%2.2f.mat', ...
        mut_prot,0.75));
load(file_ris_drug, 'ris_drug')

table_all_protein_excel = fullfile(folder_results, 'results_drug_on_mutation_of_Ras.xlsx');

idx_id = find(perc_clt == 0.75);

% - excel file
aux_tab_bad = table(species_names(idx_delta_mut_sort), ...
            round(delta_mut_sort, 2), round(delta_combo_sort(:, idx_id), 2), ...
            'VariableNames', {'Protein', 'delta_i_GoF_KRAS', 'delta_i_after_drug_action'});
writetable(aux_tab_bad, table_all_protein_excel)

%% P3. Reaction fluxes analysis 
folder_data = './data';
file_mim = fullfile(folder_data, 'CRC_CRN.mat');
file_ris_drug = fullfile(folder_results, ...
        sprintf('results_drug_on_mut_%s_alpha_%2.2f.mat', ...
        mut_prot,0.75));
    
% Protein(s) of interest
all_proteins = {'Ras', 'Ras_GTP'}; n_prot = numel(all_proteins);

% Load
%   - Model 
load(file_mim, 'CMIM');
species_names = CMIM.species.names;
species_names{CMIM.matrix.ind_one} = ' ';

%   - Dynamics
load(file_ris_drug, 'ris_drug')
time = ris_drug.time_mut_drug;
x_t = ris_drug.x_t_mut_drug; n_t = size(x_t, 2);
x_t(CMIM.matrix.ind_one, :) = ones(1, n_t);
rate_const = ris_drug.rate_constants_mut_drug;

f_summary = figure('units','normalized','outerposition',[0 0 0.5 1]);

for ip = 1:n_prot

    protein = all_proteins{ip};
    fprintf('******  Dealing with %s   ****** \n', protein)

%% Step 2. Find reactions in which the protein is involved
    [~, idx_protein] = ismember(protein, CMIM.species.names);
    b_complex_reactions = CMIM.matrix.B(CMIM.matrix.Z(idx_protein, :)~=0 , :);
    idx_reactions = [];
    for ir = 1:size(b_complex_reactions, 1)
        idx_reactions = [idx_reactions, find(b_complex_reactions(ir, :))];
    end

%% Step 3. Compute flux rates
    % 3.a. Read variables
    reactions2flux_rates = CMIM.reactions.reactions2flux_rates(idx_reactions);

    v = CMIM.matrix.v(idx_reactions, :);
    rate_const_eval = rate_const(v(:,1));
    species1_eval = x_t(v(:,2), :);
    species2_eval = x_t(v(:,3), :);
    S_prot = CMIM.matrix.S(idx_protein, idx_reactions)';

    % 3.b. Compute fluxes
    fluxes = S_prot .* rate_const_eval .* species1_eval .* species2_eval;
    fluxes_names = cell(numel(idx_reactions), 1);
    for i_f = 1:numel(idx_reactions)
       fluxes_names{i_f} = sprintf('%d %s    %s   %s', ...
           S_prot(i_f) , CMIM.rates.alias{v(i_f, 1)}, ...
           species_names{v(i_f, 2)},  species_names{v(i_f, 3)});
    end

    % 3.c. Compute flux rates
    [tmp, ~, tmp_idx] = unique(reactions2flux_rates);
    fluxes_rates = zeros(numel(tmp), n_t);
    fluxes_rates_name = cell(numel(tmp), 1);
    for i_f = 1:numel(tmp)
        fluxes_rates(i_f, :) = sum(fluxes(tmp_idx == i_f, :), 1);
        aux_fl = find(tmp_idx == i_f);
        if numel(aux_fl) == 1
            fluxes_rates_name{i_f} = fluxes_names{aux_fl};
        else
            fluxes_rates_name{i_f} = sprintf('%s + %s', ...
                fluxes_names{aux_fl(1)}, fluxes_names{aux_fl(2)});
        end
    end
    
    switch protein
        case 'Ras'
            idx_fl = [1, 2];
            aux_title = 'a) Ras';
            aux_lgd = {'$V_1$', '$V_2$', '$\dot{x}_{Ras}$'};
        otherwise 
            idx_fl = [5, 6];
            aux_title = 'b) Ras\_GTP'; 
            aux_lgd = {'$V_3$', '$V_4$', '$\dot{x}_{Ras\_GTP}$'};
    end
    figure(f_summary)
    subplot(2, 1, ip)
    for i_f = 1:numel(idx_fl)
        aux_i_f = idx_fl(i_f);
        semilogx(time, fluxes_rates(aux_i_f, :), 'Linewidth', 3, ...
            'Displayname', aux_lgd{i_f})
        hold on
    end
    tmp_flux = fluxes_rates; tmp_flux(idx_fl, :) = [];
    tmp_sum = sum(tmp_flux, 1);
    semilogx(time, tmp_sum, 'Linewidth', 3, ...
            'Displayname', 'Others')
    sum_all = sum(fluxes_rates, 1);
    semilogx(time, sum_all, 'k--', 'Linewidth', 3, ...
            'Displayname', aux_lgd{3})
    lgd = legend('show');
    set(lgd, 'Interpreter', 'Latex', 'Fontsize', 18)
    grid on
    xlim([time(2), time(end)])
    
    ylabel('Reaction fluxes')
    
    if ip == 2
        xlabel('Time [s]')
    end

    ax = gca; ax.FontSize = 20;
end

saveas(f_summary, ...
    fullfile(folder_figures, 'flux_rates_Ras_Ras_GTP.png'))