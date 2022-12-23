clear all
close all
clc


addpath('../funcs')


lof_mutations = {'APC', 'AKT', 'SMAD4',  'PTEN'};
gof_mutations = {'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = ['phys', gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);


%%

folder_data = '../data';
folder_results = './results';
folder_results_new = './results/new';
folder_figures = './results/figures';


aux_png_phys = 'png_%s.mat';
aux_png_phys_new = 'png_%s_new.mat';
aux_png_mut = 'png_mut_%s.mat';
aux_png_mut_new = 'png_mut_new_%s.mat';


%% Graphs

for im = 1:n_mutations
   mutation = all_mutations{im};
   cond_name{im} = mutation;

   %% Physiological network
   
   if strcmp(mutation, 'phys')
       
        % Load
        
        load(fullfile(folder_results_new, sprintf(aux_png_phys_new, mutation)), 'png_phys')
        png_phys_new = png_phys;
        load(fullfile(folder_results, sprintf(aux_png_phys, mutation)), 'png_phys')
        
        
        % Store
        
        num_trials_old(:, im) = [png_phys.num_trials]; 
        time_old(:, im) = [png_phys.elapse_time]; 
        num_trials(:, im) = [png_phys_new.num_trials];
        time(:, im) = [png_phys_new.elapse_time];
        
        
   else
       
   %% Mutated network
   
       % Load
       load(fullfile(folder_results_new, sprintf(aux_png_mut_new, mutation)), 'png_mut');
       png_mut_new = png_mut;
       load(fullfile(folder_results, sprintf(aux_png_mut, mutation)), 'png_mut');
       
       
       % Store
        
       num_trials_old(:, im) = [png_mut.num_trials];
       
       time_old(:, im) = [png_mut.elapse_time]; 
       num_trials(:, im) = [png_mut_new.num_trials];
       time(:, im) = [png_mut_new.elapse_time];
       
   end
   
end



%% Grafico tempo - numero iterazioni


for im = 1:n_mutations
    
    figure;
    plot(num_trials_old(:, im), time_old(:, im), 'r+')
    hold on;
    plot(num_trials(:, im), time(:, im), 'b+')
    hold off;
    
    legend('old','new')
    
end