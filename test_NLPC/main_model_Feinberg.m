clc
clear
close all

addpath(fullfile('..', 'funcs')) 
set(0,'defaulttextInterpreter','latex') 
set(0,'defaultAxesTickLabelInterpreter','latex');

folder_figures = './figures';

%% Step 1. Model definition
eps = 3; % free parameter for the rate constants
c = 1; % Free parameter for the conservation law.

% 1.1. Stoichiometric matrix
fein_S = [-1,  1; ...
           1, -1];
[n_species, n_reactions] = size(fein_S);

% 1.2. Flux vector
ind_one = n_species+1;
fein_v = [1, 1, 1; ...
         2, 1, 2];

% 1.3. Conservation laws
cons_laws.laws = f_compute_semipositive_conservations(fein_S);
fein_Nl = cons_laws.laws;
idx_basic_species = [1];
rho = c;

% 1.4. Rate constants
fein_k =[1, (eps-1)]';

% 1.5. Analytic solutions
x_eq_teo_1 = [0; c];
x_eq_teo_2 = c/eps*[(eps-1); 1]';

%% Step 2. Run NLPC
x0_a = 0.01*(1:70);
x0_all = [x0_a; c-x0_a];
% x0_all = [x0_a; rand(1, size(x0_a, 2))];
n_runs = size(x0_all, 2);

max_counter = 500;
proj = 0;

x_eq_fein = zeros(2, n_runs);
for ir = 1:n_runs
    x0 = x0_all(:, ir);
    ris = f_NLPC_restart(x0, fein_k, fein_S, fein_Nl, rho, idx_basic_species, ...
    fein_v, ind_one, max_counter, proj);
    x_eq_fein(:, ir) = ris.x;
end

%% Step 3. check results
solutions = zeros(n_runs, 1);
for ir = 1:n_runs
   [~, solutions(ir)] = min([norm(x_eq_teo_1 - x_eq_fein(:, ir)), ...
       norm(x_eq_teo_2 - x_eq_fein(:, ir))]); 
end

f_fein = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
plot(x0_a, solutions, '.', 'Markersize', 22, 'Color', [0, 100, 0]/256)
ylim([0.5, 2.5])
yticks([1, 2])
yticklabels({'Solution 1', 'Solution 2'})
xlim([0, 0.7])
xlabel('Initial concentration $x_{A,0}$')
set(gca, 'Fontsize', 22)
saveas(f_fein, fullfile(folder_figures, 'ex_multiple_sols.png'))


