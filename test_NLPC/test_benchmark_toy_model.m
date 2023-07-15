% To do:
% 1) Referenze per k e c;
% 2) Provare una c con valori molto piccoli

clc
clear
close all

addpath(fullfile('..', 'funcs')) 

%% Step 1. Define the model
% (Ref.
% https://epubs.siam.org/doi/abs/10.1137/080719820?journalCode=smjmap)
% 1.a. Stoichiometric matrix
toy_S = [-2,  2,  0,  0,  0,  0; ... 
          1, -1, -1,  1,  0,  0; ... 
          0,  0, -1,  1,  0,  0; ...
          0,  0,  1, -1, -1,  1; ...
          0,  0,  0,  0, -1,  1; ...
          0,  0,  0,  0,  1, -1];
[n_species, n_reactions] = size(toy_S);
     
% 1.b. Flux vector
ind_one = size(toy_S, 1)+1;
toy_v = [1, 1, 1; ...
         2, 2, ind_one; ...
         3, 2, 3; ...
         4, 4, ind_one; ...
         5, 4, 5; ...
         6, 6, ind_one];

% 1.c. Conservation laws
cons_laws = f_compute_semipositive_conservations(toy_S);
idx_basic_species = [1; 3; 5];
Nl = zeros(size(cons_laws));
for ibs = 1:numel(idx_basic_species)
    aux = idx_basic_species(ibs);
    Nl(ibs, :) = cons_laws(cons_laws(:, aux)==1, :);
end

% 1.d. Rate constants
toy_k =[1, 1, 1, 1, 1, 2]';

% 1.e. Store information
toy_MIM.matrix.S = toy_S; 
toy_MIM.matrix.v = toy_v;
toy_MIM.matrix.ind_one = ind_one;
toy_MIM.matrix.Nl = Nl;
toy_MIM.rates.std_values = toy_k;
toy_MIM.species.is_constant = zeros(n_species, 1);

% 1.f. Stoichiometric compatibility class
% rho = [10; 7; 10^-7];
% rho = [10; 7; 10];
rho = [10^-7; 7; 10];

%% Step 1. Simulate the dynamics of the toy model
x0_dyn = zeros(1, size(toy_MIM.matrix.S, 1));
x0_dyn(idx_basic_species) = rho;

max_t = 100000;

[time, x_t] = ode15s(@(t_, x_) f_odefun_MIM(...
    t_, x_, toy_MIM.rates.std_values, toy_MIM, 'Sv'), [0 max_t], x0_dyn);
x_t = x_t';
x_eq = x_t(:, end)';

%% Step 2. Randomly choose 50 initial points
n_points = 50;
x0_all = zeros(n_points, size(toy_MIM.matrix.S, 1));
for ip = 1:n_points
    x0_all(ip, :) = f_draw_from_ssurf(toy_MIM.matrix.Nl, rho, ...
        idx_basic_species, [-3, 3]);
end

%% Step 3. Compute NLPC solutions
max_counter = 500; proj = 0;
x_nlpc_all = zeros(n_points, size(toy_MIM.matrix.S, 1));
err_rel_nlpc = zeros(n_points, size(toy_MIM.matrix.S, 1));
time_nlpc = zeros(n_points, 1);

for ip = 1:n_points
    time_init = tic;
    aux_ris = f_NLPC_restart(x0_all(ip, :)', toy_MIM.rates.std_values, ...
        toy_MIM.matrix.S, toy_MIM.matrix.Nl, rho, idx_basic_species, ...
        toy_MIM.matrix.v, ind_one, max_counter, proj);
    x_nlpc_all(ip, :) = aux_ris.x';
    time_nlpc(ip) = toc(time_init);
    err_rel_nlpc(ip, :) = (x_nlpc_all(ip, :) - x_eq) ./ x_eq;
    clear time_init
end

%% Step 4. Compute IPOPT solutions
jacobian_v = f_compute_analytic_jacobian_v(toy_MIM.matrix.v, n_species, ind_one);

auxdata = {toy_MIM.rates.std_values, ...     % rate constants
        idx_basic_species, ... % indeces denoting basic species
        toy_MIM.matrix.Nl, ...
        rho, ...
        toy_MIM.matrix.S, ...
        toy_MIM.matrix.v, ...
        ind_one, ...
        jacobian_v};
options.auxdata = auxdata;

n_vars = n_species;
options.lb = zeros(1, n_vars); 
options.ub = Inf*ones(1, n_vars);  % Upper bound on the variables.

% Set the IPOPT options.
options.ipopt.jac_d_constant   = 'no';
options.ipopt.hessian_constant = 'no';
options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter         =  100;
options.ipopt.tol              = 1e-5;

options.ipopt.linear_solver = 'mumps';
funcs.objective         = @objective;
funcs.gradient          = @gradient;
options.ipopt.hessian_approximation      = 'limited-memory';
options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
options.ipopt.print_level           = 0;

x_ipopt_all = zeros(n_points, size(toy_MIM.matrix.S, 1));
err_rel_ipopt = zeros(n_points, size(toy_MIM.matrix.S, 1));
time_ipopt = zeros(n_points, 1);

for ip = 1:n_points
    time_init = tic;
    [x_ipopt_all(ip, :), info] = ipopt_auxdata(x0_all(ip, :),funcs,options);
    time_ipopt(ip) = toc(time_init);
    err_rel_ipopt(ip, :) = (x_ipopt_all(ip, :) - x_eq) ./ x_eq;
    clear time_init
end

%% Step 5. Compute IPOPT+Hessian solutions
funcs.hessian           = @hessian;
funcs.hessianstructure  = @hessianstructure;

options.ipopt = rmfield(options.ipopt, 'hessian_approximation');

x_ipopt_hess_all = zeros(n_points, size(toy_MIM.matrix.S, 1));
time_ipopt_hess = zeros(n_points, 1);
err_rel_ipopt_hess = zeros(n_points, size(toy_MIM.matrix.S, 1));

for ip = 1:n_points
    time_init = tic;
    [x_ipopt_hess_all(ip, :), info] = ipopt_auxdata(x0_all(ip, :),funcs,options);
    time_ipopt_hess(ip) = toc(time_init);
    err_rel_ipopt_hess(ip, :) = (x_ipopt_hess_all(ip, :) - x_eq) ./ x_eq;
    clear time_init
end


%% Plots
fx_nlpc = zeros(n_points, 1);
fx_ipopt = zeros(n_points, 1);
fx_ipopt_hess = zeros(n_points, 1);

for ip = 1:n_points
    fx_nlpc(ip) = objective(x_nlpc_all(ip, :), auxdata);
    fx_ipopt(ip) = objective(x_ipopt_all(ip, :), auxdata);
    fx_ipopt_hess(ip) = objective(x_ipopt_hess_all(ip, :), auxdata);
end

figure
subplot(3, 1, 1)
plot(1:n_points, time_nlpc, 'b*', 'MarkerSize', 10)
hold on
plot(1:n_points, time_ipopt, 'r*', 'MarkerSize', 10)
plot(1:n_points, time_ipopt_hess, 'g*', 'MarkerSize', 10)
legend('NLPC',  'IPOPT', 'IPOPT+Hessian')
xlabel('Runs')
ylabel('Elapsed time')
subplot(3, 1, 2)
plot(1:n_points, mean(err_rel_nlpc, 2), 'b*', 'MarkerSize', 10)
hold on
plot(1:n_points,  mean(err_rel_ipopt, 2), 'r*', 'MarkerSize', 10)
plot(1:n_points,  mean(err_rel_ipopt_hess, 2), 'g*', 'MarkerSize', 10)
legend('NLPC',  'IPOPT', 'IPOPT+Hessian')
xlabel('Runs')
ylabel('Relative Error')
subplot(3, 1, 3)
semilogy(1:n_points, mean(fx_nlpc, 2), 'b*', 'MarkerSize', 10)
hold on
semilogy(1:n_points,  mean(fx_ipopt, 2), 'r*', 'MarkerSize', 10)
semilogy(1:n_points,  mean(fx_ipopt_hess, 2), 'g*', 'MarkerSize', 10)
legend('NLPC',  'IPOPT', 'IPOPT+Hessian')
xlabel('Runs')
ylabel('0.5 ||f(x)||')


%% Function for IPOPT
function f = objective(x,auxdata)
  x = x';
  [rate_constants, idx_basic_species, Nl, rho, S, v, ind_one, ~] = ...
      deal(auxdata{:});
  
  f_vec = f_evaluate_mim(rate_constants, x, idx_basic_species, ... 
                             Nl, rho, S, v, ind_one);
  f = 0.5*norm(f_vec)^2;
end

function g = gradient(x,auxdata)
  x = x';
  [rate_constants, idx_basic_species, Nl, rho, S, v, ind_one, jacobian_v] = ...
      deal(auxdata{:});
  J_x = f_evaluate_jacobian(rate_constants, x, ...
                    S, idx_basic_species, jacobian_v, Nl);
  F_x = f_evaluate_mim(rate_constants, x, idx_basic_species, ... 
                             Nl, rho, S, v, ind_one);
  g = J_x'*F_x;
end

function h = hessian(x, sigma, lambda, auxdata)
    x = x';
    [rate_constants, idx_basic_species, Nl, rho, S, v, ind_one, ~] = ...
      deal(auxdata{:});
    n_species = size(S, 1);

    jac_v = f_compute_analytic_jacobian_v(v, n_species, ind_one);
    J_x = f_evaluate_jacobian(rate_constants, x, ...
                    S, idx_basic_species, jac_v, Nl);
    F_x = f_evaluate_mim(rate_constants, x, idx_basic_species, ... 
                             Nl, rho, S, v, ind_one);

    der_F = f_compute_F_2der(rate_constants, x, ...
                    v, S, Nl);  
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
