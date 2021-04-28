clc
clear
close all

%% 
% This code allows to reproduce the figures on the effect of one or
% multiple mutations on the CRC-CRN presented in the paper
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

% 1.1. Considered mutations
all_proteins = {'Ras', 'APC', 'SMAD4', 'TP53'};

lof_mutations = {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT'};
gof_mutations = {'Raf', 'Ras', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};

% 1.2. Folders and files
folder_results = fullfile('.', './results');
folder_figures = fullfile('.', './figures');

file_ris_phys = fullfile(folder_results, 'results_physiological.mat');
file_species_names = fullfile(folder_results, 'CRC_CRN_species_names.mat');
file_conservation_laws = fullfile(folder_results, 'CRC_CRN_conservation_laws.mat');
aux_file_ris_mutation = 'results_mutation_%s_perc_%1.1f.mat';

excel_mut = fullfile(folder_results, ...
           'single_gene_mutations_affected_proteins_all_mut.xlsx');

% 1.3. Parameters for the analysis
perc = 0;
thresh_eff = 0.03;

%% Step 2. Load general variables
load(file_species_names, 'species_names')
load(file_conservation_laws, 'Nl')

%   - Physiological cell
load(file_ris_phys, 'ris_phys')
x_eq = ris_phys.x_eq;
n_species = numel(x_eq);

%% Plot 1. Mutation of a single gene.
% P1.1. Load results mutated cell
f_eff_mut_log = figure('units','normalized','outerposition',[0 0 0.75 1]);
for ip = 1:numel(all_proteins)

    protein = all_proteins{ip};
    switch protein        
        case [lof_mutations, lof_mutations_type2]
            type = 'lof';
        case gof_mutations
            type = 'gof';
    end
    
    file_ris_mutation = fullfile(folder_results, ...
                sprintf(aux_file_ris_mutation, protein, perc));
    load(file_ris_mutation, 'ris_mutated')

    x_eq_mutated = ris_mutated.x_xemut_eq;
    delta = (x_eq_mutated - x_eq) ./ x_eq;

    if strcmp(protein, 'Ras')
       idx_high_delta = find(abs(delta)>5*10^3);
       disp('Species')
       disp(species_names(idx_high_delta))
       disp('Delta:')
       disp(delta(idx_high_delta))
       disp('Eq phys')
       disp(x_eq(idx_high_delta))
       disp('Eq mutated')
       disp(x_eq_mutated(idx_high_delta))
    end
    
    % Find the set of species effected by the mutation
    idx_eff_species = find(abs(delta)>thresh_eff);
    [~, aux_sort] = sort(delta(idx_eff_species));
    idx_eff_species = idx_eff_species(aux_sort);
    perc_eff = numel(idx_eff_species) *100 / n_species; 
    
    % - Find which of these species are involved in the conservation law of P
    [~, idx_protein] = ismember(protein, species_names);
    cons_law_protein = find(Nl(:, idx_protein));
    species_cons_law_protein = find(Nl(cons_law_protein, :));

    [~, tmp_idx] = ismember(species_cons_law_protein, idx_eff_species);
    is_in_cl = zeros(numel(idx_eff_species), 1);
    is_in_cl(tmp_idx) = 1;

    % - Figure for the main text
    figure(f_eff_mut_log)
    subplot(numel(all_proteins), 1, ip)
    plot(delta, 'Linewidth', 3.5, 'color', [14, 77, 146]./256)
    
    my_symlog('y')
   
    xlim([0, n_species])
    if ip == 1
        ylim([-1.15, 5.01])
    elseif ip == 3 || ip == 4
        ylim([-1.15, 1.15])
    else
        ylim([-1.15, Inf])
    end
    
    if ip < 4
        set(gca,'xticklabel',{[]}, 'FontSize', 20)
    else
        set(gca, 'FontSize', 20)
    end
        
    ylabel('$\delta_i$', 'Fontsize', 35, 'Interpreter', 'Latex')
    if ip == numel(all_proteins)
        xlabel('Proteins $i$', 'Fontsize', 35, 'Interpreter', 'Latex')
    end
    
    
    switch protein
        case 'Ras'
            aux_title = 'GoF KRAS    (a)';
        case 'APC'
            aux_title = 'LoF APC      (b)';
        case 'SMAD4'
            aux_title = 'LoF SMAD4 (c)';
        case 'TP53'
            aux_title = 'LoF TP53     (d)';
    end
    t1 = title(aux_title, 'Fontsize', 27);
    pos = get(t1, 'Position');
    pos(1) = 380;
    set(t1, 'Position', pos)
    
    % - Table for the supplementary material
    if strcmp(type, 'lof')
        idx_eff_species = idx_eff_species(~is_in_cl);
    end
    switch protein
        case 'Ras'
            aux_range = 'A2';
        case 'APC'
            aux_range = 'D2';
        case 'SMAD4'
            aux_range = 'G2';
        case 'TP53'
            aux_range = 'J2';
    end
    aux_tab_mut = table(species_names(idx_eff_species), round(delta(idx_eff_species), 2), ...
        'VariableNames', {'Protein', 'delta_i'});
    writetable(aux_tab_mut, excel_mut, 'Range', aux_range)

    clear ris_mutated x_eq_mutated delta f_eff_mut_sing fileID
    
    % - Percentage of the affected protein 
    fprintf('Percentage of proteins significantly affected by mutatio of %s = %2.2f \n', ...
        protein, perc_eff)
    
end

saveas(f_eff_mut_log, fullfile(folder_figures, 'single_gene_mutations.tiff'))


%% Plot 2. singe and multiple gene mutations affecting TP53
protein_set = {'AKTP', 'AKT','MDM2P', 'MDM2'}; % Proteins for plot
protein_set_name = {'MDM2', 'p-MDM2', 'AKT', 'p-AKT'};
protein_set2 = {'TP53'}; % Protein for test
all_mut_protein = {'PI3K', 'PTEN', 'Ras', 'Raf', 'AKT'};
all_combo = { {'Raf', 'PTEN'}, {'Raf', 'PI3K'}, {'Ras', 'PI3K'}, ....
              {'Ras', 'PTEN'},{'AKT', 'PTEN'},  {'AKT', 'PI3K'}};

n_mut = numel(all_mut_protein);
n_mut_combo = numel(all_combo);

% Find indeces of the protein of interest 
[isin_protein, idx_protein] = ismember(protein_set, species_names);
idx_protein_inv = idx_protein(end:-1:1);

[isin_protein2, idx_protein2] = ismember(protein_set2, species_names);

% Deal with results from single-gene mutations
conc_after_mut = zeros(numel(idx_protein), n_mut);
delta_mut = zeros(numel(idx_protein), n_mut);
delta_mut_p53 = zeros(numel(idx_protein2), n_mut);
for im = 1:n_mut    
    % - Load
    protein = all_mut_protein{im};
    file_ris_mutation = fullfile(folder_results, ...
                sprintf(aux_file_ris_mutation, protein, perc));
    load(file_ris_mutation, 'ris_mutated');
    x_eq_mut = ris_mutated.x_mut_eq;
    % - Prepare variable for plotting
    conc_after_mut(:, im) = x_eq_mut(idx_protein);
    temp_delta_mut = (x_eq_mut - x_eq) ./ x_eq;
    delta_mut(:, im) = temp_delta_mut(idx_protein);
    delta_mut_p53(:, im) = temp_delta_mut(idx_protein2);
end

% Deal with results from multiple-gene mutations
delta_mut_p53_combo = zeros(numel(idx_protein2), n_mut_combo);
for im = 1:n_mut_combo   
    % - Load
   all_proteins = all_combo{im};
   aux_file = sprintf('results_multi_mutations_%s_%s_pm1_%1.2f_pm2_%1.2f.mat', ...
            all_proteins{1}, all_proteins{2}, 0, 0);
   file_ris_mutation = fullfile(folder_results, aux_file);
   load(file_ris_mutation, 'ris_combo');
   x_eq_mut = ris_combo.x_eq_mut_comb;
   % - Prepare variable for plotting
   temp_delta_mut = (x_eq_mut - x_eq) ./ x_eq;
   delta_mut_p53_combo(:, im) = temp_delta_mut(idx_protein2);
end


% Plot

% - Physiological network
f_hist_phys = figure('units','normalized','outerposition',[0 0 0.8 1]);
subplot(2, n_mut, 1)
bar(1:numel(idx_protein), x_eq(idx_protein_inv))
set(gca, 'xtick', 1:numel(idx_protein), ...
        'xticklabel', protein_set_name, 'Fontsize', 20, ...
        'TickLabelInterpreter', 'none')
grid on
ylabel('Concentration [nM]', 'Fontsize', 19)
xtickangle(90)
t1 = title('(a)  Physiological');
p_title = get(t1, 'Position');
p_title(1) = 1.5; p_title(2) = p_title(2) + p_title(2) * 4/50;
set(t1, 'Position', p_title)

saveas(f_hist_phys, fullfile(folder_figures, 'hist_TP53_phys.png'))

% - Mutated network
f_hist_mut = figure('units','normalized','outerposition',[0 0 0.8 1]);
sp = {'(b)', '(c)', '(d)', '(e)', '(f)'};

fprintf('*** Alteration in p53 concentration *** \n')
   
for im = 1:n_mut
%   Set protein specific parameters
    protein = all_mut_protein{im};
    switch protein        
        case 'Raf'
            aux_title = 'GoF BRAF';
            y_max = 100;
            y_max_d = 6;
        case 'Ras'
            aux_title = 'GoF KRAS';
            y_max = 100;
            y_max_d = 6;
        case 'AKT'
            aux_title = 'LoF AKT';
            y_max = 1;
            y_max_d = 6;
%             y_max = Inf;
%             y_max_d =20000;
        case 'PI3K'
            aux_title = 'GoF PI3K';
            y_max = 100;
            y_max_d = 6;
        case 'PTEN'
            aux_title = 'LoF PTEN';
            y_max = 100;
            y_max_d = 6;
    end        
    
%   Pa. Concentrations
    figure(f_hist_mut)
    subplot(2, n_mut, im)
    bar(1:numel(idx_protein), conc_after_mut(numel(idx_protein):-1:1, im), ...
        'r', 'Displayname', 'Weighted')
    if im == 1
        ylabel('Concentration [nM]', 'Fontsize', 19)
    end
    grid on
    ylim([0, y_max])
    set(gca, 'Fontsize', 20, 'xticklabel', {[]})
    t1 = title(sprintf('%s  %s', sp{im}, aux_title));
    p_title = get(t1, 'Position');
    p_title(1) = 1.5; p_title(2) = p_title(2) + p_title(2) * 4/50;
    set(t1, 'Position', p_title)
    
%  Pb Delta
   figure(f_hist_mut)
   subplot(2, n_mut, im+n_mut)
   bar(1:numel(idx_protein), delta_mut(numel(idx_protein):-1:1, im), ...
        'r', 'Displayname', 'Weighted')
   hold on
   plot(get(gca, 'XLim'), [-1 -1], 'k--', 'Linewidth', 4)
   grid on
   if im == 1
   ylabel('\delta_i', 'Fontsize', 20)
   end
   set(gca, 'xtick', 1:numel(idx_protein), ...
        'xticklabel', protein_set_name, 'Fontsize', 20, ...
        'TickLabelInterpreter', 'none')
   ylim([-1.02, y_max_d])
   xtickangle(90)
   
   if strcmp(protein, 'AKT') % Zoom of AKT plot
       f_akt_zoom = figure('units','normalized','outerposition',[0 0 0.11 0.3]);
        ax = gca; ax.YAxis.Exponent = 4;
        bar(1:numel(idx_protein), delta_mut(numel(idx_protein):-1:1, im), ...
        'r', 'Displayname', 'Weighted')
        hold on
        plot(get(gca, 'XLim'), [-1 -1], 'k--', 'Linewidth', 4)
        grid on
        xlim([1.5, numel(idx_protein)+0.5])
        set(gca, 'xtick', 2:numel(idx_protein), 'xticklabel',{[]}, ...
            'Fontsize', 25, 'TickLabelInterpreter', 'none', 'YAxisLocation', 'right')
        ylim([-1.02, 0.5])
        saveas(f_akt_zoom, fullfile(folder_figures, 'hist_TP53_mut_AKT_zoom.png'))         
   end   
   
% Pc. Alteration in p53 concentration
   fprintf('%s --> %1.2f \n', aux_title, delta_mut_p53(im));
   
end
saveas(f_hist_mut, fullfile(folder_figures, 'hist_TP53_mut.png'))

% Pd. Multipe-gene mutations
for im = 1:n_mut_combo   
   all_proteins = all_combo{im};
   fprintf('%s + %s --> %1.2f \n', all_proteins{1}, all_proteins{2}, ...
                                delta_mut_p53_combo(im));
end