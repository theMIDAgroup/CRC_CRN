clc
clear
close all

rng(1003) % For reproducing the results of the manuscript

run_algo = 0;

%% Step 1. Load data
addpath(fullfile('..', 'funcs')) % CRC-CRN functions
target = fullfile('..', 'data');
file_ris = 'test_ipopt_CRN.mat';

path_mim = fullfile(target, 'CRC_CRN_nodrug.mat');
load(path_mim, 'new_CMIM'); CRN = new_CMIM;

path_dym = fullfile('./results_paper', 'dyn_phys.mat');
load(path_dym, 'dyn_phys') % Use resultf_NLPC_restart.ms from dynamic approach as reference
sol = dyn_phys(1).x;

[n_species, n_reactions] = size(CRN.matrix.S);
ind_one = n_species + 1;

rho = CRN.matrix.Nl * CRN.species.std_initial_values;
idx_basic_species = find(CRN.species.std_initial_values>0);

disp('Close to zero species')
disp(numel(find(sol<10^-3))/n_species)

switch run_algo
    
    case 1
        
    jacobian_v = f_compute_analytic_jacobian_v(CRN.matrix.v, n_species, ind_one);
    disp('mtfbwy')
    der_F =  my_f_compute_F_2der(CRN.rates.std_values, CRN.matrix.v, ...
        CRN.matrix.S, CRN.matrix.Nl, idx_basic_species);
    disp('Done')

%% Step 2. Define the initial point
x0 = sol' + 0.5*randn(1, n_species).*sol';
x0(x0 < 0) = 0;

%% Step 3. Run NLPC
max_counter = 500; proj = 0;

time_init = tic;
aux_ris = f_NLPC_restart(x0', CRN.rates.std_values, ...
        CRN.matrix.S, CRN.matrix.Nl, rho, idx_basic_species, ...
        CRN.matrix.v, ind_one, max_counter, proj);
time_nlpc = toc(time_init);

fprintf('Elsapsed time for the NLPC algorithm = %1.2f sec \n', time_nlpc)

%% Step 4. Define IPOPT options
% 4.a. Option related to the CR-CRN
rho = CRN.matrix.Nl * CRN.species.std_initial_values;
auxdata = {CRN.rates.std_values, ...     
           idx_basic_species, ... 
           CRN.matrix.Nl, ...
           rho, ...
           CRN.matrix.S, ...
           CRN.matrix.v, ...
           ind_one, ...
           jacobian_v, ...
           der_F};
options.auxdata = auxdata;

% 4.b. Variable bounds
n_vars = n_species;
options.lb = zeros(1, n_vars); 
options.ub = Inf*ones(1, n_vars);  

% 4.c. Options concerning the IPOPT algorithm
options.ipopt.jac_d_constant   = 'no';
options.ipopt.hessian_constant = 'no';
options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter         =  10000;
options.ipopt.tol              = 1e-3;

options.ipopt.linear_solver = 'mumps';

options.ipopt.hessian_approximation      = 'limited-memory';
options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
options.ipopt.bound_relax_factor = 0;
options.ipopt.print_level        = 0;

% 4.d. Objective functions
funcs.objective         = @objective;
funcs.gradient          = @gradient;

all_mu_max = [1, 10^-2, 10^-4, 10^-6, 10^-8];

%% Step 5. Run IPOPT
ris_ipopt.sols = zeros(n_species, numel(all_mu_max));
ris_ipopt.time = zeros(numel(all_mu_max), 1);

for im = 1:numel(all_mu_max)
    fprintf('IPOPT %d \n', im)
    options.ipopt.mu_max           = all_mu_max(im);
    time_init = tic;
    [aux_x, info] = ipopt_auxdata(x0,funcs,options);
    ris_ipopt.time(im) = toc(time_init);
    ris_ipopt.sols(:, im) = aux_x';
    clear aux_x time_init
end

%% Step 7. Save
 save(file_ris, 'ris_ipopt', 'all_mu_max')

    otherwise
        
%% Analyze results

%% R1. Load
load(file_ris, 'ris_ipopt', 'all_mu_max')

%% R2. Compute evaluation matrics
all_delta = zeros(size(ris_ipopt.sols));
all_norm_f = zeros(1, numel(all_mu_max));

for im = 1:numel(all_mu_max)
    
    tmp_sol = ris_ipopt.sols(:, im);
    
    all_delta(:, im) = abs(tmp_sol - sol) ./ sol;
    
    idx_basic_species = find(CRN.species.std_initial_values>0);
    rho = CRN.matrix.Nl * CRN.species.std_initial_values;
    f_vec = f_evaluate_mim(CRN.rates.std_values, tmp_sol, idx_basic_species, ... 
                           CRN.matrix.Nl, rho, CRN.matrix.S, ...
                           CRN.matrix.v, CRN.matrix.ind_one);
    all_norm_f(im) = norm(f_vec);  
    
end

%% R3. Plot
f_ipopt_crn = figure('units','normalized','outerposition',[0 0 0.7 0.8]);
subplot(2, 1, 1)
loglog(all_mu_max, all_norm_f, 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 10)
set(gca, 'Fontsize', 15)
xlabel('$\mu^*$', 'Interpreter', 'Latex', 'Fontsize', 20)
ylabel('Accuracy - $||\textbf{f}(\textbf{x})||$', 'Interpreter', 'Latex', 'Fontsize', 20)
grid on


subplot(2, 1, 2)
for im = 1:numel(all_mu_max)
    loglog(sol, all_delta(:, im), '.', 'Markersize', 25, 'Displayname',  sprintf('%1.0e', all_mu_max(im)))
    hold on
end
 set(gca, 'Fontsize', 15)
hleg = legend('show', 'Location', 'eastoutside', 'Fontsize', 18);
htitle = get(hleg,'Title');
set(htitle,'String','\mu^*')
xlabel('$\textbf{x}^{eq}$', 'Interpreter', 'Latex', 'Fontsize', 20)
ylabel('Relative error', 'Interpreter', 'Latex', 'Fontsize', 20)
grid on
saveas(f_ipopt_crn, 'ipopt_crn.png')



        
    

end


%%
function f = objective(x,auxdata)
  x = x';
  [rate_constants, idx_basic_species, Nl, rho, S, v, ind_one, ~, ~] = ...
      deal(auxdata{:});
  
  f_vec = f_evaluate_mim(rate_constants, x, idx_basic_species, ... 
                             Nl, rho, S, v, ind_one);
  f = 0.5*norm(f_vec)^2;
end

function g = gradient(x,auxdata)
  x = x';
  [rate_constants, idx_basic_species, Nl, rho, S, v, ind_one, jacobian_v, ~] = ...
      deal(auxdata{:});
  J_x = f_evaluate_jacobian(rate_constants, x, ...
                    S, idx_basic_species, jacobian_v, Nl);
  F_x = f_evaluate_mim(rate_constants, x, idx_basic_species, ... 
                             Nl, rho, S, v, ind_one);
  g = J_x'*F_x;
end

function h = hessian(x, sigma, lambda, auxdata)
    x = x';
    [rate_constants, idx_basic_species, Nl, rho, S, v, ind_one, ~, der_F] = ...
      deal(auxdata{:});
    n_species = size(S, 1);

    jac_v = f_compute_analytic_jacobian_v(v, n_species, ind_one);
    J_x = f_evaluate_jacobian(rate_constants, x, ...
                    S, idx_basic_species, jac_v, Nl);
    F_x = f_evaluate_mim(rate_constants, x, idx_basic_species, ... 
                             Nl, rho, S, v, ind_one);
    h = J_x'*J_x;
    for i=1:n_species
        for j=i:n_species
            h(i,j) = h(i,j) + F_x'*squeeze(der_F(:,i,j));
        end
    end
    h = sparse(tril(h));
end

function h_s = hessianstructure(auxdata)
    n_species = size(auxdata{5}, 1);
    h = tril(ones(n_species,n_species));
    h_s = sparse(h);
end