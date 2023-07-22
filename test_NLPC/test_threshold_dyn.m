clc
clear
close all

test_max_t = 0;

%% Path, folder & load
addpath(fullfile('..', 'funcs'))
target = fullfile('..', 'data');
folder_results = fullfile('.');
folder_results_paper = './results_paper';
aux_file_results = 'test_dyn_threshold_%d.mat';

path_mim = fullfile(target, 'CRC_CRN_nodrug.mat'); % Network
load(path_mim, 'new_CMIM'); CRN = new_CMIM;

path_dyn_phys = fullfile(folder_results_paper, 'dyn_phys.mat');
load(path_dyn_phys, 'x0_all')

%% Define general parameters of the network
rates_phys = CRN.rates.std_values;
n_species = size(x0_all, 1);
x0_phys = CRN.species.std_initial_values;
idx_basic_species = find(x0_phys>0);
Nl = CRN.matrix.Nl;
rho_phys = Nl*x0_phys;
S_phys = CRN.matrix.S;
v_phys = CRN.matrix.v;
ind_one = size(CRN.matrix.S, 1) + 1;

%% Solve the dynamical system for the physiological network with a 
%% a larger time interval
if test_max_t
    
max_t_all = [2.5*10^9]; n_runs = size(max_t_all, 2);
ix = 1; x0_dyn = x0_all(:, ix);

fprintf('Initial point = %d \n', ix)
elapse_times = zeros(1, n_runs);
sols = zeros(n_species, n_runs);
for ir = 1:n_runs
    fprintf('max_t n = %d \n', ir)
    max_t = max_t_all(ir);
    time_init = tic;
    [aux_time, aux_sol] = ode15s(@(t_, x_) f_odefun_MIM(...
            t_, x_, rates_phys, CRN, 'Sv'), [0 max_t], x0_dyn);
    elapse_times(ir) = toc(time_init);
    aux_sol = aux_sol';
    sols(:, ir) = aux_sol(:, end);

    if ir == n_runs
        ris.time = aux_time;
        ris.sol = aux_sol;
    end

    clear aux_time aux_sol
end

save(fullfile(folder_results, sprintf(aux_file_results, ix)), ...
    'elapse_times', 'sols', 'max_t_all', 'ris')

end

%% Plots

%% P1. Effects of the time-interval
ix = 1; aux_p = 244;

file_ris = fullfile(folder_results, sprintf(aux_file_results, ix));
load(file_ris)

n_times = numel(ris.time);
norm_sol = zeros(n_times, 1);
% Evaluate the norm of f on the trajectory of the ODE solution.
for it = 1:n_times
    norm_sol(it) = norm(f_evaluate_mim(rates_phys, ris.sol(:, it), ...
        idx_basic_species, Nl, rho_phys, S_phys, v_phys, ind_one));
end

f_max_t = figure('units','normalized','outerposition',[0 0 0.7 0.6]);
t_min = 10^-5;
subplot(2, 1, 1)
loglog(ris.time, norm_sol, 'k', 'linewidth', 3)
xlim([t_min, ris.time(end)])
ylabel({'Accuracy',  '||f(x(t))||'})
grid on
set(gca, 'Fontsize', 15)

subplot(2, 1, 2)
semilogx(ris.time, ris.sol(aux_p, :), 'k', 'linewidth', 3)
hold on
rectangle('Position', [0.75*10^8, 8.4*10^-5, ris.time(end)-0.8*10^8, 0.8*10^-5], ...
    'EdgeColor', 'red', 'Linewidth', 2)
xlim([t_min, ris.time(end)])
xlabel('Time t [s]')
ylabel({'Concentration', sprintf('%s', CRN.species.names{aux_p})}, ...
    'Interpreter', 'Latex')
set(gca, 'Fontsize', 15)
grid on
saveas(f_max_t, fullfile(folder_results, 'test_threshold_dyn.jpg'))


f_zoom = figure('units','normalized','outerposition',[0 0 0.35 0.15]);
semilogx(ris.time, ris.sol(aux_p, :), 'k', 'linewidth', 3)
xlim([0.8*10^8, ris.time(end)])
set(gca, 'Fontsize', 15)
saveas(f_zoom, fullfile(folder_results, 'test_threshold_dyn_zoom.jpg'))







