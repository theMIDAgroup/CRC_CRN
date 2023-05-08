clear all
close all
clc

%% 
% This code allows to reproduce the figures on the effect of partial, single
% and combined mutations of KRAS (GoF) and PTEN (LoF) affecting the CRC-CRN
% presented in the paper
%
% 'In-silico modelling of the mitogen-activated kinase (MAPK) pathway in
% colorectal cancer: mutations and targeted therapy.'
% Sommariva et al. 2023, preprint at 
% https://www.biorxiv.org/content/10.1101/2023.04.18.537359v1.full.pdf
%
% The results of the simulations required to run this code are provided 
% in the folder './results_paper/mutations'
%%
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex'); 

%% Path & Folders
target_folder = '../data/ci_servono';
folder_results = './results_paper';
folder_figures = './figures_paper';
file_species_names = fullfile(target_folder, 'CRC_CRN_species_names.mat');
file_ris_phys = fullfile(folder_results, 'nlpc_phys.mat');

%% Data
mut_prot = {'Ras', 'PTEN', 'Ras_PTEN', 'PTEN_Ras'};
name_mut_prot = {'KRAS', 'PTEN'};
quantity = {'single', 'single', 'multi', 'multi'};
lof_mutations = {'TBRII', 'SMAD4', 'Cadh', 'APC', 'PTEN', 'AKT'};
gof_mutations = {'Raf', 'Ras', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};

for i=1:2
    switch mut_prot{i}
        case [lof_mutations, lof_mutations_type2]
            type_mut{i} = 'LoF';
        case gof_mutations
            type_mut{i} = 'GoF';
    end
    title_fig{i} = [type_mut{i}, ' ',  name_mut_prot{i}];
end
type_mut{3} = [type_mut{1}, '_' , type_mut{2}];
type_mut{4} = [type_mut{2}, '_' , type_mut{1}];
title_fig{3} = [title_fig{1}, ' + ', title_fig{2}];%strcat(title_fig{1}, ' ', title_fig{2});
title_fig{4} = [title_fig{2}, ' + ',  title_fig{1}];

%% Load and store data
load(file_species_names, 'species_names')
load(file_ris_phys, 'nlpc_phys')

x_eq_phys = nlpc_phys(1).x;
n_species = numel(species_names);
profile_mut = zeros(n_species, 4);
for i=1:4
    switch(quantity{i})
        case 'single'
            file_mut_eq = fullfile(folder_results, 'mutations', ...
                sprintf('nlpc_mut_%s.mat', mut_prot{i}));
            load(file_mut_eq, 'nlpc_mut');
            profile_mut(:, i) = nlpc_mut(1).x;
        case 'multi'
            file_mut_eq = fullfile(folder_results, 'mutations', ...
                sprintf('nlpc_mut_%s.mat', mut_prot{i}));
            load(file_mut_eq, 'nlpc_combo');
            profile_mut(:, i) = nlpc_combo.x_eq_mut_comb(1:n_species);
    end
    clear nlpc_mut nlpc_combo
end

%% Compute delta
delta = zeros(n_species, 4);
for i=1:4
    delta(:, i) = (profile_mut(:, i) - x_eq_phys) ./ x_eq_phys;
end

%% Partial single mutations

mut_prot_p = 'Ras'; type_mut_p = {'GoF'}; name_mut_prot_p = {'KRAS'};
perc = [0, 0.3, 0.6];

profile_mut_portion = zeros(n_species, numel(perc));

for i=1:numel(perc)
    file_mut_eq_portion = fullfile(folder_results, 'mutations', ...
        sprintf('nlpc_mut_%s_perc_%1.1f.mat', name_mut_prot_p{:}, perc(i)));
    load(file_mut_eq_portion, 'nlpc_mut'); 
    profile_mut_portion(:,i) = nlpc_mut(1).x;
end

delta_portion = zeros(n_species, numel(perc));
for i=1:numel(perc)
    delta_portion(:,i) =  (profile_mut_portion(:,i) - x_eq_phys) ./ x_eq_phys;
end

%% Figure
f_mut_log = figure('units','normalized','outerposition',[0 0 7.0 1]);
n_fig = {'(B)', '(C)', '(D)'};
for i=1:3
   sp = subplot(5,1,i+2); 
   pos = get(sp, 'Position' );
   pos(2) = pos(2) -i* 0.015;  
   plot(delta(:,i), 'Linewidth', 2.5, 'color', [14, 77, 146]./256);
   if i==4
       xlabel('Proteins $i$', 'Fontsize', 30, 'Interpreter', 'Latex')
   end
   xlabel(' ')
   ylabel('$\delta_i$', 'Fontsize', 30, 'Interpreter', 'Latex')
   my_symlog('y')
   xlim([0, 420])
   t1 = title({'';[n_fig{i}, '  ', title_fig{i}]}, 'Fontsize', 24);
   t1.Units = 'Normalize'; 
   t1.Position(1) = 0; % use negative values (ie, -0.1) to move further left
   t1.HorizontalAlignment = 'left';
   set(sp, 'Position', pos ) ;
end
xlabel('Proteins $i$', 'Fontsize', 30, 'Interpreter', 'Latex')

idx_delta = 1:n_species;
perc_mut = perc*100;

figure(f_mut_log)
sp = subplot(5, 1, 1:2);
pos = get(sp, 'Position' );
pos(2) = pos(2);

for j = 1:numel(perc)
    plot(delta_portion(idx_delta,j), 'Linewidth',2, 'Displayname', ...
        sprintf('$ %3.0f \\%%$', perc_mut(j)));
    
    if j<3
        hold on
    end
end
xlabel(' '+newline+' ')
xlim([0, 420])

ylabel('$\delta_i$', 'Fontsize', 30, 'Interpreter', 'Latex')
t3 = title('(A) Partial GoF KRAS', 'Fontsize', 24);
t3.Units = 'Normalize'; 
t3.Position(1) = 0;
t3.HorizontalAlignment = 'left';
my_symlog('y')
hold off

lgd = legend('show');
set(lgd, 'Fontsize', 18, 'Location', 'Northoutside', 'orientation', 'horizontal', 'Interpreter', 'Latex')
pos_lgd = get(lgd, 'Position');
lgd.Position = [pos_lgd(1)+0.08 pos_lgd(2)+0.045 pos_lgd(3) pos_lgd(4)];
set(sp, 'Position', pos ) ;
saveas(f_mut_log, fullfile(folder_figures, ...
   sprintf('mut_%s_%s.jpg', name_mut_prot{1}, name_mut_prot{2})))

%% Compute matrix x_eq
x_eq = [x_eq_phys profile_mut(:,1:3)];

%% Histograms Ras, Ras_GTP, Raf, p-Raf, MEK, pp-MEK, ERK, pp-ERK
species = {'Ras', 'Ras_GTP', 'Raf', 'Raf*', 'MEK', 'MEKPP', 'ERK', 'ERKPP'};
name_species = {'KRAS', 'k-Ras_GTP', 'Raf', 'p-Raf', 'MEK', 'p-p-MEK', 'ERK', 'p-p-ERK'};
idx_title = {'(A)', '(B)', '(C)', '(D)', '(E)', '(F)', '(G)', '(H)'};
[~, idx_species] = ismember(species, species_names);
n_sp = numel(species);
status = {'phys', [type_mut{1}, ' ', name_mut_prot{1}], [type_mut{2}, ' ', name_mut_prot{2}],...
    [type_mut{1}, ' ', name_mut_prot{1}, ' + ', type_mut{2}, ' ', name_mut_prot{2}]};
title_2 = {'         (A): k-Ras', ' (B): k-Ras\_GTP', '           (C): Raf', '      (D): p-Raf', '          (E): MEK',...
    '   (F): p-p-MEK', '         (G): ERK', '  (H): p-p-ERK'};

f_hist_species = figure('units','normalized','outerposition',[0 0 0.7 1.5]); 
bar_color = parula(4);
for i=1:n_sp
    sp = subplot(3,4,i);
    pos = get(sp, 'Position');
    pos(1) = pos(1)+mod(i-1,4)*0.02;
    set(sp, 'Position', pos)
    for j=1:4
        bar(j, x_eq(idx_species(i),j), 'FaceColor', bar_color(j,:), 'Displayname', status{j});
        hold on;
    end
    set(gca, 'xtick', 1:4, 'xticklabel', ' ',  ...
        'Fontsize', 20, ...
        'TickLabelInterpreter', 'Latex')
    xtickangle(90)
    if ismember(i, [1, 3, 4, 5, 7])
    t2 = title(title_2{i}, 'Fontsize', 20, 'horizontalAlignment', 'Left');
    t2.HorizontalAlignment = 'center';
    else
        t2 = title(title_2{i}, 'Fontsize', 20);
    end
    if ((i==1) || (i==5))
        ylabel('Concentration [nM]', 'Fontsize', 19, 'Interpreter', 'Latex')
    end
end
grid on
lgd = legend('show');
newPosition = [0.43 0.23 0.2 0.25];
set(lgd,'Position', newPosition, 'orientation', 'horizontal');

%% Save
saveas(f_hist_species, fullfile(folder_figures, 'Hist_mut_species.jpg'));

%% Table distance phys - GoF k-Ras

x_eq_phys_gofras = x_eq(:, 1:2);
x_eq_gofras_diff = abs(x_eq_phys_gofras(:,2) - x_eq_phys_gofras(:,1));
n_max = 10;
[max_diff, idx_max_diff] = maxk(x_eq_gofras_diff, n_max);
delta_gofras = delta(idx_max_diff,1);
species_names_tex = strrep(species_names, '_', '\_');

% figure;
% plot(x_eq_gofras_diff);

table_diff = ['\begin{tabular}{|l|l|l|}' newline '\hline', ...
    '\multicolumn{1}{|l|}{\large{\textbf{$\mathbf{i}$-th protein}}} ', newline,...
    '& \large{$\mathbf{|\tilde{x}^e_i - x^e_i|}$} & \large{$\mathbf{\delta_i}$} \\ \hline', newline];

for i=1:n_max
    table_diff = [table_diff, '\multicolumn{1}{|l|} {', char(species_names_tex(idx_max_diff(i))), ...
        '} & ', sprintf('%0.2f',max_diff(i)), ' & ', sprintf('%0.2f', delta_gofras(i)), '\\ \hline', newline];
end
table_diff = [table_diff, '\end{tabular}'];

format bank;
disp(table_diff)
