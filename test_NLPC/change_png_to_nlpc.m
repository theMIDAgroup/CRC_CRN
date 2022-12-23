clc
clear
close all


lof_mutations = {'APC', 'AKT', 'SMAD4', 'PTEN'};
gof_mutations = {'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = ['phys', gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);


%% Folder and files

folder_data = '../data';
folder_results_old = './results/18_11';
folder_results = './results/22_12';
folder_results_dyn = './results/dinamica';

aux_dyn_phys = 'dyn_%s.mat';
aux_nlpc_phys = '/nlpc_%s.mat';
aux_dyn_mut = 'dyn_mut_%s.mat';
aux_nlpc_mut = '/nlpc_mut_%s.mat';


for im = 1:n_mutations
   mutation = all_mutations{im};

   %% Physiological network
   if strcmp(mutation, 'phys')
       
       load(fullfile(folder_results_old, sprintf(aux_nlpc_phys, mutation)), 'png_phys')
       nlpc_phys = png_phys;
       save(fullfile(folder_results, 'nlpc_phys.mat'), 'nlpc_phys');

   else
       
       load(fullfile(folder_results_old, sprintf(aux_nlpc_mut, mutation)), 'png_mut');
       nlpc_mut = png_mut;
       save(fullfile(folder_results, ...
            sprintf('nlpc_mut_%s.mat', mutation)), 'nlpc_mut')

   end
end