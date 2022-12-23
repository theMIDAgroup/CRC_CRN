clear all
close all
clc


%% Path, Folders & Load

folder_results2 = './results/18_11';
folder_results1 = './results/18_11_ort';
folder_figures = './figures';

load(fullfile(folder_results1, 'png_ort_phys.mat'), 'png_ort_phys');
load(fullfile(folder_results2, 'png_phys.mat'), 'png_phys');


%% Plot results - updated components & conditioning number - convergence case


for i=1:50
   fig1 = figure; 
   hold on
   plot(png_ort_phys(i).upd_components(end).n, 'LineWidth', 1.5);
   plot(png_phys(i).upd_components(end).n, 'LineWidth', 1.5);
   hold off
%    fig2 = figure; 
%    hold on
%    plot(png_ort_phys(i).cond_number(end).n);
%    plot(png_phys(i).cond_number(end).n);
%    hold off
   %saveas(fig1, fullfile(folder_figures, sprintf('upd_comp_%d', i)), 'png');
   %saveas(fig2, fullfile(folder_figures, sprintf('cond_number_%d', i)), 'png');
   close all;
end


%% Plot results - updated components & conditioning number - non-convergence case


% for i=1:50
%    fig1 = figure; 
%    hold on
%    plot(png_ort_phys(i).upd_components(end).n);
%    plot(png_phys(i).upd_components(end).n);
%    hold off
% %    fig2 = figure; 
% %    hold on
% %    plot(png_ort_phys(i).cond_number(end).n);
% %    plot(png_phys(i).cond_number(end).n);
% %    hold off
%    saveas(fig1, fullfile(folder_figures, sprintf('upd_comp_%d', i)), 'png');
%    %saveas(fig2, fullfile(folder_figures, sprintf('cond_number_%d', i)), 'png');
%    close all;
% end