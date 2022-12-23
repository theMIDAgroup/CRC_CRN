clear all 
close all
clc


% Here we create graphs for comparing results obtained using PNG with
% the orthogonal projector or the non-projector operator.
% We are comparing the two methods by looking at:
% - number (percentage) of null components; iteration after iteration
% - condition number of the jacobian matrix
% - number of updated components.


%% Path & Load

folder_data = '../data';
folder_results1 = './results/18_11_ort'; %31_10_orthproj_40alpha_10-2';
folder_results2 = './results/18_11'; %31_10_classicmethod3mod_40alpha_10-2';
folder_figures = './figures';

aux_png_phys = 'png_%s.mat';
aux_png_ort_phys = 'png_ort_%s.mat';
aux_png_mut = 'png_mut_%s.mat';
aux_png_ort_mut = 'png_ort_mut_%s.mat';

% Mutations
lof_mutations = {'APC',  'AKT', 'SMAD4', 'PTEN'};
gof_mutations = {'Ras', 'Raf', 'PI3K', 'BetaCatenin'};
lof_mutations_type2 = {'TP53'};
all_mutations = ['phys', gof_mutations, lof_mutations, lof_mutations_type2];
n_mutations = numel(all_mutations);


%% Define what we work with


% for i=1:n_runs
%     
%     zeri_x0 = sum(png_ort_phys(i).x0 == 0);
%     
%     upd_comp_ort(i).n = png_ort_phys(i).upd_components;
%     upd_comp_class(i).n = png_phys(i).upd_components;
%     zeri_ort(i).n = [zeri_x0 png_ort_phys(i).zeri(1).n];
%     zeri_class(i).n = [zeri_x0 png_phys(i).zeri(1).n];
%     cond_number_ort(i).n = png_ort_phys(i).cond_number(1).n;
%     cond_number_class(i).n = png_phys(i).cond_number(1).n;
% end


n_species = 419;

n_runs = 50;

%% Mean max zeros 

max_zeri_ort = zeros(1, n_runs);
max_zeri_class = zeros(1, n_runs);

mean_max_zeri_ort = zeros(1, n_mutations);
mean_max_zeri_class = zeros(1, n_mutations);
mean_max_zeri_ort_conv = zeros(1, n_mutations);
mean_max_zeri_ort_nonconv = zeros(1, n_mutations);
mean_max_zeri_class_conv = zeros(1, n_mutations);
mean_max_zeri_class_nonconv = zeros(1, n_mutations);

std_max_zeri_ort = zeros(1, n_mutations);
std_max_zeri_class = zeros(1, n_mutations);
std_max_zeri_ort_conv = zeros(1, n_mutations);
std_max_zeri_ort_nonconv = zeros(1, n_mutations);
std_max_zeri_class_conv = zeros(1, n_mutations);
std_max_zeri_class_nonconv = zeros(1, n_mutations);

conv_png_class = zeros(1, n_mutations);
conv_png_ort = zeros(1, n_mutations);

mean_max_zeri_ort_perc = zeros(1, n_mutations);
std_max_zeri_ort_perc = zeros(1, n_mutations);
mean_max_zeri_class_perc = zeros(1, n_mutations);
std_max_zeri_class_perc = zeros(1, n_mutations);

for im = 1:n_mutations
   mutation = all_mutations{im};
   cond_name{im} = mutation;
   
   %% Physiological network
   if strcmp(mutation, 'phys') 
       
       % Load
       load(fullfile(folder_results2, sprintf(aux_png_phys, mutation)), 'png_phys')
       load(fullfile(folder_results1, sprintf(aux_png_ort_phys, mutation)), 'png_ort_phys')
       
       for i=1:n_runs 
           zeri_x0 = sum(png_ort_phys(i).x0 == 0);
           zeri_ort(i).n = [zeri_x0 png_ort_phys(i).zeri(1).n];
           zeri_class(i).n = [zeri_x0 png_phys(i).zeri(1).n];

           max_zeri_ort(i) = max(zeri_ort(i).n);
           max_zeri_class(i) = max(zeri_class(i).n);
           
           conv_png_class(i) = (png_phys(i).num_trials == 1);
           conv_png_ort(i) = (png_ort_phys(i).num_trials == 1);
       end
       
   else
       
       %% Mutated Network
       
       % Load
       load(fullfile(folder_results2, sprintf(aux_png_mut, mutation)), 'png_mut');
       load(fullfile(folder_results1, sprintf(aux_png_ort_mut, mutation)), 'png_ort_mut');

   
   for i=1:n_runs 
       zeri_x0 = sum(png_ort_mut(i).x0 == 0);
       zeri_ort(i).n = [zeri_x0 png_ort_mut(i).zeri(1).n];
       zeri_class(i).n = [zeri_x0 png_mut(i).zeri(1).n];

       max_zeri_ort(i) = max(zeri_ort(i).n);
       max_zeri_class(i) = max(zeri_class(i).n);
       
       conv_png_class(i) = (png_mut(i).num_trials == 1);
       conv_png_ort(i) = (png_ort_mut(i).num_trials == 1);
   end

   end
   
   max_zeri_ort_conv = max_zeri_ort(conv_png_ort == 1);
   max_zeri_ort_nonconv = max_zeri_ort(conv_png_ort == 0);
   max_zeri_class_conv = max_zeri_class(conv_png_class == 1);
   max_zeri_class_nonconv = max_zeri_class(conv_png_class == 0);
       
   mean_max_zeri_ort(im) = mean(max_zeri_ort);
   mean_max_zeri_class(im) = mean(max_zeri_class);
   mean_max_zeri_class_conv(im) = mean(max_zeri_class_conv);
   mean_max_zeri_class_nonconv(im) = mean(max_zeri_class_nonconv);
   mean_max_zeri_ort_conv(im) = mean(max_zeri_ort_conv);
   mean_max_zeri_ort_nonconv(im) = mean(max_zeri_ort_nonconv);
   
   std_max_zeri_ort(im) = std(max_zeri_ort);
   std_max_zeri_class(im) = std(max_zeri_class);
   std_max_zeri_class_conv(im) = std(max_zeri_class_conv);
   std_max_zeri_class_nonconv(im) = std(max_zeri_class_nonconv);
   std_max_zeri_ort_conv(im) = std(max_zeri_ort_conv);
   std_max_zeri_ort_nonconv(im) = std(max_zeri_ort_nonconv);
   
end

mean_max_zeri_ort_perc_conv = mean_max_zeri_ort_conv*(100/n_species);
mean_max_zeri_ort_perc_nonconv = mean_max_zeri_ort_nonconv*(100/n_species);
std_max_zeri_ort_perc_conv = std_max_zeri_ort_conv*(100/n_species);
std_max_zeri_ort_perc_nonconv = std_max_zeri_ort_nonconv*(100/n_species);
mean_max_zeri_class_perc_conv = mean_max_zeri_class_conv*(100/n_species);
mean_max_zeri_class_perc_nonconv = mean_max_zeri_class_nonconv*(100/n_species);
std_max_zeri_class_perc_conv = std_max_zeri_class_conv*(100/n_species);
std_max_zeri_class_perc_nonconv = std_max_zeri_class_nonconv*(100/n_species);

mean_max_zeri_ort_perc = mean_max_zeri_ort*(100/n_species);
std_max_zeri_ort_perc = std_max_zeri_ort*(100/n_species);
mean_max_zeri_class_perc = mean_max_zeri_class*(100/n_species);
std_max_zeri_class_perc = std_max_zeri_class*(100/n_species);

figure;
plot(mean_max_zeri_ort_perc, '*');
hold on
plot(mean_max_zeri_class_perc, '*');
hold off;
title '% of mean of max number of zeros during PNG iterations';
legend('PNG ortogonale', 'PNG classico');
set(legend, 'Location', 'East')

% Plot divided by convergence - if reached or not


%% Flag: Newton or gradient? 

max_flag_ort = zeros(1, n_runs);
max_flag_class = zeros(1, n_runs);

conv_png_class = zeros(1, n_mutations);
conv_png_ort = zeros(1, n_mutations);

sum_flag_ort = zeros(1, n_mutations);
sum_flag_class = zeros(1, n_mutations);

sum_flag_ort_perc = zeros(1, n_mutations);
sum_flag_class_conv_perc = zeros(1, n_mutations);
sum_flag_class_perc = zeros(1, n_mutations);

grad_perc_ort = zeros(n_runs, n_mutations);
grad_perc_class = zeros(n_runs, n_mutations);


for im = 1:n_mutations
   mutation = all_mutations{im};
   cond_name{im} = mutation;
   
   %% Physiological network
   if strcmp(mutation, 'phys') 
       
       % Load
       load(fullfile(folder_results2, sprintf(aux_png_phys, mutation)), 'png_phys')
       load(fullfile(folder_results1, sprintf(aux_png_ort_phys, mutation)), 'png_ort_phys')
       
       for i=1:n_runs 
           flag_ort(i).n = png_ort_phys(i).flag(1).n;
           flag_class(i).n = png_phys(i).flag(1).n;
           
           conv_png_class(i) = (png_phys(i).num_trials == 1);
           conv_png_ort(i) = (png_ort_phys(i).num_trials == 1);

           % conto quante volte applico il gradiente
           sum_flag_ort(i) = sum(flag_ort(i).n);
           sum_flag_class(i) = sum(flag_class(i).n);
           
           % conto quante volte applico Newton
           sum_flag_ort_perc(i) = sum_flag_ort(i) * 100 / length(flag_ort(i).n);
           sum_flag_class_perc(i) = sum_flag_class(i) * 100 / length(flag_class(i).n);
           
       end
       
   else
       
       %% Mutated Network
       
       % Load
       load(fullfile(folder_results2, sprintf(aux_png_mut, mutation)), 'png_mut');
       load(fullfile(folder_results1, sprintf(aux_png_ort_mut, mutation)), 'png_ort_mut');

   
       for i=1:n_runs 
           flag_ort(i).n = png_ort_mut(i).flag(1).n;
           flag_class(i).n =  png_mut(i).flag(1).n;

           conv_png_class(i) = (png_mut(i).num_trials == 1);
           conv_png_ort(i) = (png_ort_mut(i).num_trials == 1);
           
           % conto quante volte applico il gradiente
           sum_flag_ort(i) = sum(flag_ort(i).n);
           sum_flag_class(i) = sum(flag_class(i).n);
           
           % conto quante volte applico Newton
           sum_flag_ort_perc(i) = sum_flag_ort(i) * 100 / length(flag_ort(i).n);
           sum_flag_class_perc(i) = sum_flag_class(i) * 100 / length(flag_class(i).n);
           
       end
        
       
   end
   
   grad_perc_ort(:, im) = sum_flag_ort_perc';
   grad_perc_class(:, im) = sum_flag_class_perc';
   
   %figure
   figure;
   hAxes = gca;
   subplot(1,2,1);
   xlim([1 50]);
   plot(find(conv_png_ort == 1), grad_perc_ort(conv_png_ort == 1, im), 'b*');
   hold on;
   plot(find(conv_png_ort == 0), grad_perc_ort(conv_png_ort == 0, im), 'r*');
   hold off
   legend('conv', 'non conv');
   title(string(all_mutations(im)) + ' - ort. proj.')
   subplot(1,2,2);
   plot(find(conv_png_class == 1), grad_perc_class(conv_png_class == 1, im), 'b*');
   hold on
   plot(find(conv_png_class == 0), grad_perc_class(conv_png_class == 0, im), 'r*');
   legend('conv', 'non conv');
   title(string(all_mutations(im)) + ' - op. P')

end

 

% figure;
% plot(mean_max_flag_ort_perc, '*');
% hold on
% plot(mean_max_flag_class_perc, '*');
% hold off;
% title '% of mean of max number of zeros during PNG iterations';
% legend('PNG ortogonale', 'PNG classico');
% set(legend, 'Location', 'East')



%% Mean max conditioning number

max_cond_ort = zeros(1, n_runs);
max_cond_class = zeros(1, n_runs);
max_log_cond_ort = zeros(1, n_runs);

mean_max_cond_ort = zeros(1, n_mutations);
mean_max_cond_class = zeros(1, n_mutations);
mean_max_cond_class_conv = zeros(1, n_mutations);
mean_max_cond_class_nonconv = zeros(1, n_mutations);
mean_max_cond_ort_conv = zeros(1, n_mutations);
mean_max_cond_ort_nonconv = zeros(1, n_mutations);

std_max_cond_ort = zeros(1, n_mutations);
std_max_cond_class = zeros(1, n_mutations);
std_max_cond_class_conv = zeros(1, n_mutations);
std_max_cond_class_nonconv = zeros(1, n_mutations);
std_max_cond_ort_conv = zeros(1, n_mutations);
std_max_cond_ort_nonconv = zeros(1, n_mutations);

std_max_log_cond_ort = zeros(1, n_mutations);
std_max_log_cond_class = zeros(1, n_mutations);
std_max_log_cond_class_conv = zeros(1, n_mutations);
std_max_log_cond_class_nonconv = zeros(1, n_mutations);
std_max_log_cond_ort_conv = zeros(1, n_mutations);
std_max_log_cond_ort_nonconv = zeros(1, n_mutations);

conv_png_class = zeros(1, n_mutations);
conv_png_ort = zeros(1, n_mutations);

for im = 1:n_mutations
   mutation = all_mutations{im};
   cond_name{im} = mutation;

   %% Physiological network
   if strcmp(mutation, 'phys') 
       
       % Load
       load(fullfile(folder_results2, sprintf(aux_png_phys, mutation)), 'png_phys')
       load(fullfile(folder_results1, sprintf(aux_png_ort_phys, mutation)), 'png_ort_phys')
       
       for i=1:n_runs 
           cond_ort(i).n = [png_ort_phys(i).cond_number(1).n];
           cond_class(i).n = [png_phys(i).cond_number(1).n];
           log_cond_ort(i).n = log10(cond_ort(i).n);
           log_cond_class(i).n = log10(cond_class(i).n);
         
           max_cond_ort(i) = max(cond_ort(i).n);
           max_cond_class(i) = max(cond_class(i).n);
           max_log_cond_ort(i) = max(log_cond_ort(i).n);
           max_log_cond_class(i) = max(log_cond_class(i).n);
           
           
           conv_png_class(i) = (png_mut(i).num_trials == 1);
           conv_png_ort(i) = (png_ort_mut(i).num_trials == 1);
       end
       
       
   else
       
       %% Mutated Network
       
       % Load
       load(fullfile(folder_results2, sprintf(aux_png_mut, mutation)), 'png_mut');
       load(fullfile(folder_results1, sprintf(aux_png_ort_mut, mutation)), 'png_ort_mut');

   
       for i=1:n_runs 
           cond_ort(i).n = png_ort_mut(i).cond_number(1).n;
           cond_class(i).n = png_mut(i).cond_number(1).n;
           log_cond_ort(i).n = log10(cond_ort(i).n);
           log_cond_class(i).n = log10(cond_class(i).n);

           max_cond_ort(i) = max(cond_ort(i).n);
           max_cond_class(i) = max(cond_class(i).n);
           max_log_cond_ort(i) = max(log_cond_ort(i).n);
           max_log_cond_class(i) = max(log_cond_class(i).n);
           
           conv_png_class(i) = (png_mut(i).num_trials == 1);
           conv_png_ort(i) = (png_ort_mut(i).num_trials == 1);
       end
    
   end
   
   
   
   
       figure;
       semilogy(max_cond_ort .* conv_png_ort, '*');
       hold on
       semilogy(max_cond_ort .* (1-conv_png_ort), '*');
       hold off
       legend('conv', 'non conv');
       title(mutation);
       
   
   
   if sum(conv_png_ort)>0
      max_cond_ort_conv = max_cond_ort(conv_png_ort == 1);
   else
      max_cond_ort_conv = max_cond_ort(conv_png_ort == 1);
   end
   if sum(conv_png_ort)>0
      max_cond_class_conv = max_cond_ort(conv_png_class == 1);
   else
      max_cond_class_conv = max_cond_ort(conv_png_class == 1);
   end
   max_cond_ort_nonconv = max_cond_ort(conv_png_ort == 0);
   max_cond_class_nonconv = max_cond_class(conv_png_class == 0);
   
   max_log_cond_ort_conv = max_log_cond_ort(conv_png_ort == 1);
   max_log_cond_ort_nonconv = max_log_cond_ort(conv_png_ort == 0);
   max_log_cond_class_conv = max_log_cond_class(conv_png_class == 1);
   max_log_cond_class_nonconv = max_log_cond_class(conv_png_class == 0);
   
   mean_max_cond_ort(im) = mean(max_cond_ort);
   mean_max_cond_class(im) = mean(max_cond_class);   
   mean_max_cond_ort_conv(im) = mean(max_cond_ort_conv);
   mean_max_cond_ort_nonconv(im) = mean(max_cond_ort_nonconv);
   mean_max_cond_class_conv(im) = mean(max_cond_class_conv);
   mean_max_cond_class_nonconv(im) = mean(max_cond_class_nonconv);
   
   mean_max_log_cond_ort(im) = mean(max_log_cond_ort);
   mean_max_log_cond_class(im) = mean(max_log_cond_class);   
   mean_max_log_cond_ort_conv(im) = mean(max_log_cond_ort_conv);
   mean_max_log_cond_ort_nonconv(im) = mean(max_log_cond_ort_nonconv);
   mean_max_log_cond_class_conv(im) = mean(max_log_cond_class_conv);
   mean_max_log_cond_class_nonconv(im) = mean(max_log_cond_class_nonconv);
   
   std_max_cond_ort(im) = std(max_cond_ort);
   std_max_cond_class(im) = std(max_cond_class);   
   std_max_cond_ort_conv(im) = std(max_cond_ort_conv);
   std_max_cond_ort_nonconv(im) = std(max_cond_ort_nonconv);
   std_max_cond_class_conv(im) = std(max_cond_class_conv);
   std_max_cond_class_nonconv(im) = std(max_cond_class_nonconv);
   
   std_max_log_cond_ort(im) = std(max_log_cond_ort);
   std_max_log_cond_class(im) = std(max_log_cond_class);   
   std_max_log_cond_ort_conv(im) = std(max_log_cond_ort_conv);
   std_max_log_cond_ort_nonconv(im) = std(max_log_cond_ort_nonconv);
   std_max_log_cond_class_conv(im) = std(max_log_cond_class_conv);
   std_max_log_cond_class_nonconv(im) = std(max_log_cond_class_nonconv);
   
end

% flag
mean_max_zeri_ort_perc_conv = mean_max_zeri_ort_conv*(100/n_species);
mean_max_zeri_ort_perc_nonconv = mean_max_zeri_ort_nonconv*(100/n_species);
mean_max_zeri_class_perc_conv = mean_max_zeri_class_conv*(100/n_species);
mean_max_zeri_class_perc_nonconv = mean_max_zeri_class_nonconv*(100/n_species);

std_max_zeri_ort_perc_conv = std_max_zeri_ort_conv*(100/n_species);
std_max_zeri_ort_perc_nonconv = std_max_zeri_ort_nonconv*(100/n_species);
std_max_zeri_class_perc_conv = std_max_zeri_class_conv*(100/n_species);
std_max_zeri_class_perc_nonconv = std_max_zeri_class_nonconv*(100/n_species);



% figure;
% semilogy(mean_max_cond_ort, '*');
% hold on
% semilogy(mean_max_cond_class, '*');
% hold off;
% title '% of mean of max conditioning number during PNG iterations';
% legend('PNG ortogonale', 'PNG classico');
% set(legend, 'Location', 'East')


%% Create Latex table 1


string_table1 = "\begin{table}[ht]" + newline + " \centering " + newline + "	\begin{tabular}{c*{11}{c|}}	\cline{3-12}" ...
    + "	& & \textbf{phys} & \textbf{k-Ras}  & \textbf{Raf} " + newline + " & \textbf{PI3K} & " ...
    + "\textbf{Betacatenin} & \textbf{APC} & \textbf{AKT} &  \textbf{SMAD4} & \textbf{PTEN}" + newline + "" ...
    + "& \textbf{p53} \\ " + newline +	"\hline	" + newline + "\multicolumn{1}{|c|}{\% of} 	&";
string_table2 = " \end{tabular}	" + newline + "  \end{table}";
string_results = "\textbf{Operator} $\mathcal{P}$";


for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf('%0.2f', mean_max_zeri_class_perc(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%0.2f', std_max_zeri_class_perc(i));
end
string_results = string_results + "\\ " + newline + " \cline{2-12}  " + newline + "\multicolumn{1}{|c|}{Zeros}  & \textbf{Ort. proj.}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf('%0.2f', mean_max_zeri_ort_perc(i)) ...
         + ' ' + char(177) + ' '  + sprintf('%0.2f', std_max_zeri_ort_perc(i));
end
string_results = string_results + "\\ " + newline + " \hline " + newline + " \multicolumn{1}{|c|}{Cond.} & \textbf{Operator $\mathcal{P}$}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf("%2.2e", mean_max_cond_class(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%2.2e', std_max_cond_class(i));
end
string_results = string_results + "\\ " + newline + " \cline{2-12}  " + newline + "\multicolumn{1}{|c|}{Number}  & \textbf{Ort. proj.}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf("%2.2e", mean_max_cond_ort(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%2.2e', std_max_cond_ort(i));
end
string_results = string_results + "\\  " + newline + " \hline  " + newline + "";

string_table = string_table1 + string_results + string_table2;

disp(string_table)



disp(newline)
disp(newline)
disp(newline)



%% Create Latex table 2


string_table1 = "\begin{table}[ht]" + newline + " \centering " + newline + "	\begin{tabular}{c*{11}{c|}}	\cline{3-12}" ...
    + "	& & \textbf{phys} & \textbf{k-Ras}  & \textbf{Raf} " + newline + " & \textbf{PI3K} & " ...
    + "\textbf{Betacatenin} & \textbf{APC} & \textbf{AKT} &  \textbf{SMAD4} & \textbf{PTEN} " + newline + "" ...
    + "& \textbf{p53} \\ " + newline +	"\hline	" + newline + "\multicolumn{1}{|c|}{\% of} 	&";
string_table2 = " \end{tabular}	" + newline + "  \end{table}";
string_results = "\textbf{Operator} $\mathcal{P}$";


for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf('%0.2f', mean_max_zeri_class_perc(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%0.2f', std_max_zeri_class_perc(i));
end
string_results = string_results + "\\ " + newline + " \cline{2-12}  " + newline + "\multicolumn{1}{|c|}{Zeros}  & \textbf{Ort. proj.}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf('%0.2f', mean_max_zeri_ort_perc(i)) ...
         + ' ' + char(177) + ' '  + sprintf('%0.2f', std_max_zeri_ort_perc(i));
end
string_results = string_results + "\\ " + newline + " \hline " + newline + " \multicolumn{1}{|c|}{Cond.} & \textbf{Operator $\mathcal{P}$}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf("%0.2f", mean_max_log_cond_class(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%0.2f', std_max_log_cond_class(i));
end
string_results = string_results + "\\ " + newline + " \cline{2-12}  " + newline + "\multicolumn{1}{|c|}{Number}  & \textbf{Ort. proj.}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf("%0.2f", mean_max_log_cond_ort(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%0.2f', std_max_log_cond_ort(i));
end
string_results = string_results + "\\  " + newline + " \hline  " + newline + "";

string_table = string_table1 + string_results + string_table2;

disp(string_table)





disp(newline)
disp(newline)
disp(newline)



%% Create Latex table 3 - phys only

string_table1 = "\begin{table}[ht]" + newline + " \centering " + newline + "	\begin{tabular}{c|*{2}{c|}}	\cline{2-3} " ...
    + "	& \textbf{Operator} $\mathcal{P}$ & Ort.Proj. \\ " + newline +	"\hline	" + newline + "\multicolumn{1}{|c|}{\% of zeros} 	&";
string_table2 = " \end{tabular}	" + newline + "  \end{table}";
string_results = sprintf('%0.2f', mean_max_zeri_class_perc(1)) +" " + char(177) + " " + sprintf('%0.2f', std_max_zeri_class_perc(1)) ...
    + "& " + sprintf('%0.2f', mean_max_zeri_ort_perc(1)) +" " + char(177) + " " + sprintf('%0.2f', std_max_zeri_ort_perc(1)) + " \\" ...
    + "\hline" + newline + "\multicolumn{1}{|c|}{Log. cond number} & " + sprintf('%d', mean_max_log_cond_class(1)) +" " ...
    + char(177) + " " + sprintf('%d', std_max_log_cond_class(1)) ...
    + "& " + sprintf('%d', mean_max_log_cond_ort(1)) +" " + char(177) + " " + sprintf('%d', std_max_log_cond_ort(1));


string_results = string_results + "\\  " + newline + " \hline  " + newline + "";

string_table = string_table1 + string_results + string_table2;

disp(string_table)



%% Create Latex table 4 - only convergence

string_table1 = "\begin{table}[ht]" + newline + " \centering " + newline + "	\begin{tabular}{c*{11}{c|}}	\cline{3-12}" ...
    + "	& & \textbf{phys} & \textbf{k-Ras}  & \textbf{Raf} " + newline + " & \textbf{PI3K} & " ...
    + "\textbf{Betacatenin} & \textbf{APC} & \textbf{AKT} &  \textbf{SMAD4} & \textbf{PTEN} " + newline + "" ...
    + "& \textbf{p53} \\ " + newline +	"\hline	" + newline + "\multicolumn{1}{|c|}{\% of} 	&";
string_table2 = " \end{tabular}	" + newline + "  \end{table}";
string_results = "\textbf{Operator} $\mathcal{P}$";


for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf('%0.2f', mean_max_zeri_class_perc_conv(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%0.2f', std_max_zeri_class_perc_conv(i));
end
string_results = string_results + "\\ " + newline + " \cline{2-12}  " + newline + "\multicolumn{1}{|c|}{Zeros}  & \textbf{Ort. proj.}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf('%0.2f', mean_max_zeri_ort_perc_conv(i)) ...
         + ' ' + char(177) + ' '  + sprintf('%0.2f', std_max_zeri_ort_perc_conv(i));
end
string_results = string_results + "\\ " + newline + " \hline " + newline + " \multicolumn{1}{|c|}{Cond.} & \textbf{Operator $\mathcal{P}$}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf("%0.2f", mean_max_log_cond_class_conv(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%0.2f', std_max_log_cond_class_conv(i));
end
string_results = string_results + "\\ " + newline + " \cline{2-12}  " + newline + "\multicolumn{1}{|c|}{Number}  & \textbf{Ort. proj.}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf("%0.2f", mean_max_log_cond_ort_conv(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%0.2f', std_max_log_cond_ort_conv(i));
end
string_results = string_results + "\\  " + newline + " \hline  " + newline + "";

string_table = string_table1 + string_results + string_table2;

disp(string_table)


disp(newline)
disp(newline)
disp(newline)



%% Create Latex table 5 - only non convergence

string_table1 = "\begin{table}[ht]" + newline + " \centering " + newline + "	\begin{tabular}{c*{11}{c|}}	\cline{3-12}" ...
    + "	& & \textbf{phys} & \textbf{k-Ras}  & \textbf{Raf} " + newline + " & \textbf{PI3K} & " ...
    + "\textbf{Betacatenin} & \textbf{APC} & \textbf{AKT} &  \textbf{SMAD4} & \textbf{PTEN} " + newline + "" ...
    + "& \textbf{p53} \\ " + newline +	"\hline	" + newline + "\multicolumn{1}{|c|}{\% of} 	&";
string_table2 = " \end{tabular}	" + newline + "  \end{table}";
string_results = "\textbf{Operator} $\mathcal{P}$";


for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf('%0.2f', mean_max_zeri_class_perc_nonconv(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%0.2f', std_max_zeri_class_perc_nonconv(i));
end
string_results = string_results + "\\ " + newline + " \cline{2-12}  " + newline + "\multicolumn{1}{|c|}{Zeros}  & \textbf{Ort. proj.}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf('%0.2f', mean_max_zeri_ort_perc_nonconv(i)) ...
         + ' ' + char(177) + ' '  + sprintf('%0.2f', std_max_zeri_ort_perc_nonconv(i));
end
string_results = string_results + "\\ " + newline + " \hline " + newline + " \multicolumn{1}{|c|}{Cond.} & \textbf{Operator $\mathcal{P}$}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf("%0.2f", mean_max_log_cond_class_nonconv(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%0.2f', std_max_log_cond_class_nonconv(i));
end
string_results = string_results + "\\ " + newline + " \cline{2-12}  " + newline + "\multicolumn{1}{|c|}{Number}  & \textbf{Ort. proj.}";
for i = 1:numel(all_mutations)
    string_results = string_results + " & " + sprintf("%0.2f", mean_max_log_cond_ort_nonconv(i)) ...
         + ' ' + char(177) + ' ' + sprintf('%0.2f', std_max_log_cond_ort_nonconv(i));
end
string_results = string_results + "\\  " + newline + " \hline  " + newline + "";

string_table = string_table1 + string_results + string_table2;

disp(string_table)



disp(newline)
disp(newline)
disp(newline)







% %% Create graphs - zeros at first attempt
% 
% 
% for i=1:n_runs
%     f = figure('units','normalized','outerposition',[0 0 1 0.5]);
%     subplot(1,2,1);
%     xlabel('Iteration k');
%     ylabel('Null components of x_k (%)');
%     
%     hold on;
%     plot(100*(zeri_ort(i).n/n_species),  'Markersize', 3, 'LineWidth', 2, 'color', [0 0.5 0]);
%     plot(100*(zeri_class(i).n/n_species), 'Markersize', 3, 'LineWidth', 3, 'color', [1 0 0]);
% %     symlog('y');
% %     semilogy(100*(zeri_ort(i).n/n_species), 'LineWidth', 1.8, 'color', [0 0.5 0]);
% %     hold on;
% %     semilogy(100*(zeri_class(i).n/n_species), 'LineWidth', 1.8, 'color', [1 0 0]);
%     legend('Orthogonal projector','Non-projector', 'Location', 'east');
%     xlim([0 n_runs0]);
% %     ylim([-1, n_runs])
%     grid on
%     hold off;
%     
%     subplot(1,2,2);
%     semilogy(cond_number_ort(i).n,  'Markersize', 3, 'LineWidth', 2, 'color', [0 0.5 0]);
%     hold on;
%     semilogy(cond_number_class(i).n, 'Markersize', 3, 'LineWidth', 3, 'color', [1 0 0]);
%     xlabel('Iteration k');
%     ylabel('Condition number of J_F(x_k)');
%     legend('Orthogonal projector','Non-projector', 'Location', 'east');
%     xlim([0 n_runs0]);
%     grid on
%     hold off;
%     
%     %saveas(f, sprintf('./figures/fig_phys_%d', i), 'png');
%     %close all
% end


%% Create graphs - jacobian condition number at first attempt

% 
% for i=1:n_runs
%     f = figure('units','normalized','outerposition',[0 0 1 0.5]);
%     semilogy(cond_number_ort(i).n,  'Markersize', 3, 'LineWidth', 2, 'color', [0 0.5 0]);
%     hold on;
%     semilogy(cond_number_class(i).n, 'Markersize', 3, 'LineWidth', 3, 'color', [1 0 0]);
%     xlabel('Iteration k');
%     ylabel('Condition number of J(x_k)');
%     legend('Orthogonal projector','Non-projector', 'Location', 'east');
%     xlim([0 n_runs0]);
%     %title('First attempt');
%     hold off;
%     %saveas(f, sprintf('figures/fig_phys_condn_%d', i), 'png');
%     %close all
% end



%% Create graphs - jacobian condition number in convergence case


% for i=1:n_runs
%     f = figure;
%     semilogy(cond_number_ort(i).n, 'LineWidth', 1.8);
%     hold on;
%     semilogy(cond_number_class(i).n, 'LineWidth', 1.8);
%     xlabel('Iteration');
%     ylabel('Jacobian condition number');
%     legend('Orthogonal projector','Non-projector', 'Location', 'east');
%     xlim([0 n_runs0]);
%     title('First attempt');
%     hold off;
%     %saveas(f, sprintf('fig_phys_condn_%d', i), 'png');
%     %close all
% end


%% Create graphs - updated components in convergence case


% for i=1:n_runs
%     f = figure;
%     xlabel('Iteration');
%     ylabel('Number of updated components');
%     hold on;
%     plot(upd_comp_ort(i).n(end).n, 'b', 'LineWidth', 1.7);
%     plot(upd_comp_class(i).n(end).n, 'r', 'LineWidth', 1.7);
%     legend('Orthogonal projector','Non-projector');
%     xlim([0 n_runs0]);
%     title('Last attempt');
%     hold off;
%     %saveas(f, sprintf('fig_phys_end_%d', i), 'png');
%     %close all
% end


%% Create graphs - updated components at first attempt
% 
% close all
% for i=1:n_runs
%     f = figure;
%     xlabel('Iteration');
%     ylabel('Number of updated components');
%     hold on;
%     plot(upd_comp_ort(i).n(1).n, 'b', 'LineWidth', 1.7);
%     plot(upd_comp_class(i).n(1).n, 'r', 'LineWidth', 1.7);
%     legend('Orthogonal projector','Non-projector');
%     xlim([0 n_runs0]);
%     title('First attempt');
%     hold off;
%     saveas(f, sprintf('fig_phys_1_%d', i));
%     %close all
% end
% 

%% Save figures