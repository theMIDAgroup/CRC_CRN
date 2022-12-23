clear all
close all
clc


%% Path and folders

addpath(fullfile('..', 'funcs'))

target = fullfile('..', 'data');
folder_results = './results';


%% Data

path_mim = fullfile(target, 'CRC_CRN_nodrug.mat');
load(path_mim, 'new_CMIM'); CRN = new_CMIM;
v = CRN.matrix.v;
Nl = CRN.matrix.Nl;
idx_basic_species = find(CRN.species.std_initial_values>0);
n_species = size(CRN.matrix.S, 1);
ind_one = n_species + 1;

% 3.1. Define network parameters specific to physiological cells
rate_constants = CRN.rates.std_values;
S = CRN.matrix.S;
rho = Nl*CRN.species.std_initial_values;

%% Parameters

max_t = [1 2 3 4 5] * 1e7;


%% Computing F norm in equilibrium points


for i=1:length(max_t)
    
    name_file = sprintf('dyn_phys_%s.mat', max_t(i));
    load(fullfile(folder_results, name_file), 'dyn_phys_t', 'x0_all');
    n_x0 = size(x0_all, 2);
    
    for j=1:n_x0
        
        x_eq = dyn_phys_t(j).x;
        F_x_eq = f_evaluate_mim(rate_constants, x_eq, ...
                    idx_basic_species, Nl, rho, S, v, ind_one);
        norm_F_x_eq = norm(F_x_eq);
        
        dyn_phys_t(j).norm_F = norm_F_x_eq;
        
        %Save
        %save(sprintf('dyn_phys_%s.mat', max_t(i)), 'dyn_phys_t', 'x0_all');
        
        dyn_phys(i).norm_F(j) = norm_F_x_eq;
        
    end
    
    
    
    clear dyn_phys_t x0_all
    
end


%% Plot results

figure;
hold on
xlabel('Starting points');
ylabel('||F(x\_eq)||');
for i=1:length(max_t)
    plot(dyn_phys(i).norm_F(1:10));
end
legend('max t = 1e7', 'max t = 2e7', 'max t = 3e7', 'max t = 4e7', 'max t = 5e7');



hold off

%% SARA
addpath('~/Documents/CRC_CRN/funcs')
all_x_eq = zeros(n_species, length(max_t));
all_norm_F = zeros(1, length(max_t));
for i=1:length(max_t)
    
    name_file = sprintf('dyn_phys_%s.mat', max_t(i));
    load(fullfile(folder_results, name_file), 'dyn_phys_t', 'x0_all');
    n_x0 = size(x0_all, 2);
    
    for j=30
        
        x_eq = dyn_phys_t(j).x;
        
        all_x_eq(:, i) = x_eq;
        
        F_x_eq = f_evaluate_mim(rate_constants, x_eq, ...
                    idx_basic_species, Nl, rho, S, v, ind_one);
        norm_F_x_eq = norm(F_x_eq);
        
        all_norm_F(i) = norm_F_x_eq;
        
    end
    
    clear dyn_phys_t x0_all
    
end

name_file = sprintf('png_phys_testproj.mat');
load(fullfile(folder_results, name_file), 'png_phys');

x_eq_ref = png_phys(j).x;

delta_1 = ( all_x_eq(:, 1) - x_eq_ref) ./ x_eq_ref;
delta_2 = ( all_x_eq(:, length(max_t)) - x_eq_ref) ./ x_eq_ref;


delta = (all_x_eq(:, length(max_t)) - all_x_eq(:, 1)) ./ all_x_eq(:, 1);

F_xeq_1 = f_evaluate_mim(rate_constants, all_x_eq(:, 1), ...
                    idx_basic_species, Nl, rho, S, v, ind_one);
F_xeq_2 = f_evaluate_mim(rate_constants, all_x_eq(:, length(max_t)), ...
                    idx_basic_species, Nl, rho, S, v, ind_one);
norm_1 = norm(F_xeq_1);
norm_2 = norm(F_xeq_2);
norm_png = norm(f_evaluate_mim(rate_constants, x_eq_ref, ...
                    idx_basic_species, Nl, rho, S, v, ind_one));


figure
subplot(2, 1, 1)
semilogy(all_x_eq(:, 1), 'k', 'linewidth', 2)
hold on
semilogy(all_x_eq(:, length(max_t)), 'r', 'linewidth', 1)
subplot(2, 1, 2)
plot(delta, 'k', 'linewidth', 2)

figure
plot(delta_1, 'k', 'linewidth', 2)
hold on
plot(delta_2, 'r--', 'linewidth', 2)

figure
plot(F_xeq_1, 'k', 'linewidth', 2)
hold on
plot(F_xeq_2, 'r--', 'linewidth', 2)




