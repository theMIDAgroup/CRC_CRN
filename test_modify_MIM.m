clc
clear
close all

addpath(fullfile('./funcs'))

%% Step 1. Define general parameters and load
% 1.a. File and folders
target_folder = './data';
file_mim = fullfile(target_folder, 'CRC_CRN.mat');

% 1.b. Species to be removed
death_species = {'DRUG', 'DRUG_Raf'};

% 1.c. Load
load(file_mim)

%% Step 2. Find species to be removed
n_ds = numel(death_species);
[aux, idx_ds] = ismember(death_species, CMIM.species.names);

%% Step 3. Find reactions to be removed
% --> Only reactions having one of the death species as reactant will be
%     removed
idx_dr = [];
for is = 1:n_ds
   idx_dr = [idx_dr, find(CMIM.reactions.details_ind(:, 2) == idx_ds(is)),...
       find(CMIM.reactions.details_ind(:, 3) == idx_ds(is))];
end

%% Step 4. Define the new variable
new_CMIM.species.names = CMIM.species.names;
new_CMIM.species.names(idx_ds) = [];

new_CMIM.species.is_constant = CMIM.species.is_constant;
new_CMIM.species.is_constant(idx_ds) = [];

new_CMIM.species.std_initial_values = CMIM.species.std_initial_values;
new_CMIM.species.std_initial_values(idx_ds) = [];

new_CMIM.rates.std_values = CMIM.rates.std_values;
new_CMIM.rates.std_values(idx_dr) = [];

new_CMIM.matrix.S = CMIM.matrix.S;
new_CMIM.matrix.S(idx_ds, :) = [];
new_CMIM.matrix.S(:, idx_dr) = [];

new_CMIM.matrix.ind_one = numel(new_CMIM.species.names)+1; % --> Metterlo a -1, 0.
new_CMIM.matrix.v = CMIM.matrix.v;
new_CMIM.matrix.v(idx_dr, :) = [];
temp = new_CMIM.matrix.v(:, 2:3);
temp(temp == CMIM.matrix.ind_one) = new_CMIM.matrix.ind_one;
new_CMIM.matrix.v(:, 2:3) = temp;


new_CMIM.matrix.Nl = f_compute_semipositive_conservations(new_CMIM.matrix.S);

save(fullfile('data', 'CRC_CRN_nodrug.mat'), 'new_CMIM')

%% Step 5. Sanity check
max_t = 2.5*10^7;
disp('Solve new CMIM')
x0 = new_CMIM.species.std_initial_values;
rates = new_CMIM.rates.std_values;
[~, sol_dyn] = ode15s(@(t_, x_) f_odefun_MIM(...
            t_, x_, rates, new_CMIM, 'Sv'), [0 max_t], x0);

clear x0

disp('Solve old CMIM')
x0 = CMIM.species.std_initial_values;
rates = CMIM.rates.std_values;
[~, sol_dyn_old] = ode15s(@(t_, x_) f_odefun_MIM(...
            t_, x_, rates, CMIM, 'Sv'), [0 max_t], x0);

x_eq = sol_dyn(end, :);
x_eq_old = sol_dyn_old(end, :);
x_eq_old = x_eq_old(1:end-2);

delta = (x_eq - x_eq_old) ./ x_eq_old;

figure
plot(delta)     