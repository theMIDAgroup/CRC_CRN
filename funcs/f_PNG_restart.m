function struct = f_PNG_restart(MIM, num_run, max_counter)

% The function 'newtonGD_restart' takes the following inputs:
% - 'MIM' is a struct that contains all cell's features (in a physiological or mutated state)
% - 'ris_dynamic' is a struct that contains the results obtained working on the
% same cell, but using the dynamic approach
% - 'num_run' is an integer that indicates how many experiments we want to
% do on the same stoichiometric matrix
% - 'max_counter' is an integer that indicates how many iterations we want
% the algorithm to do every time we choose a new starting point


addpath(fullfile('..', 'functions'));

%% Step 1. Read data
% 1.a. Define path and load
target = fullfile('.', 'data');   
folder_figures = './figures';
folder_results = './results';

% 1.b. Define parameters for the analysis
[n_species, n_reactions] = size(MIM.matrix.S);
x_0_phys = MIM.species.std_initial_values;
idx_basic_species = find(x_0_phys>0);
Nl = MIM.matrix.Nl; num_cons_laws = size(Nl, 1);
rate_constants = MIM.rates.std_values;

%% Step 2. Fix stoichiometric compatibility class
%     Note: we reorder conservation vectors according to basic species 
%           (Really needed? Anyway I do not have the structure 
%                           N = [I, N2]
%           as I have not reorderd the species, i.e. the columns.)
new_order = zeros(num_cons_laws, 1);
for ii = 1:num_cons_laws
    new_order(ii) = find(Nl(:, idx_basic_species(ii)));
end
Nl = Nl(new_order, :);
rho = Nl*x_0_phys;

%% Step 3. Compute 'symbolic' form of the Jacobian matrix of v(x; k)
v = MIM.matrix.v;
ind_one = MIM.matrix.ind_one;
jacobian_v = f_compute_analytic_jacobian_v(v, n_species, ind_one);

%% Step 4. Run the algorithm 
%%      4.a. Parameters of the algorithm
toll_cond_init_point = 10^17;
tol = 1e-12; 
poss_alpha = logspace(0, -2, 20);
poss_alpha_2 = 1e3 * logspace(0, -2, 20);
poss_alpha_3 = 1e1 * logspace(0, -2, 20);
poss_alpha_2 = [poss_alpha_2  poss_alpha_3];
sigma = 10^-4;
sigma_2 = 10^-4;
FLAG = 0;
elapsed_time = 0;
z = 0;
ir = 1;

while ir < num_run + 1%for ir = 1:num_run
   
    fprintf('@@@@@@@@@@@  RUN num %d @@@@@@@@@@ \n', ir)  

%%      4.b. Initialize
%   We draw x_0 on the given stoichiometric compatibility class until
%   cond(J_F(x_0)) is small enough
disp('Initialization')
aux_cond = Inf; 
while aux_cond > toll_cond_init_point
    x_ss = f_draw_from_ssurf(Nl, rho, idx_basic_species, [-3, 3]);
    aux_cond = cond(f_evaluate_jacobian(rate_constants, x_ss, ...
                MIM.matrix.S, idx_basic_species, jacobian_v, Nl));
end
disp('Let''s start')

%%      4.c. Run newton method 
time_init = tic;

x_0_nm = x_ss;

% Define initial point
x =  x_0_nm; x(x<0) = 0;

xnew = x; counter = 0; 
F_x_new = f_evaluate_mim(rate_constants, xnew, idx_basic_species, Nl, rho, MIM);
norm_F_x_new = norm(F_x_new);

norm_F_x_nm = zeros(max_counter, 1);
diff_x_eq_nm = zeros(max_counter, 1);
est_err = zeros(max_counter, 1);
step_lengths = zeros(max_counter, 1);
det_F_x = zeros(max_counter, 1);
norm_grad = zeros(1,max_counter);
j = 1;

% 5.b. Iteration
while norm_F_x_new > tol
    
    counter = counter + 1;
    x = xnew; F_x = F_x_new; norm_F_x = norm_F_x_new;
    J_x = f_evaluate_jacobian(rate_constants, x, ...
                MIM.matrix.S, idx_basic_species, jacobian_v, Nl);
     
    if FLAG == 0 % Newton
        
        delta = -J_x \ F_x;
        ia = 1;
        
        while ia <= numel(poss_alpha) 
            
            alpha = poss_alpha(ia);
            xnew = x + alpha * delta;
            xnew(xnew<0) = x(xnew<0);
            F_x_new = f_evaluate_mim(rate_constants, xnew, idx_basic_species, Nl, rho, MIM);
            norm_F_x_new = norm(F_x_new);
        
            if norm_F_x_new <= sqrt(1-alpha*sigma)*norm_F_x
                ia = numel(poss_alpha)+1;
                FLAG = 0;        
            else
                FLAG = 1;
                ia = ia+1;
            end
        end
    
    if (ia == numel(poss_alpha)+1) && (FLAG == 1)
        xnew = x; F_x_new = F_x; norm_F_x_new = norm_F_x; counter = counter-1;
    else
        % Store some informations
        norm_F_x_nm(counter) = norm_F_x_new;
        est_err(counter) = norm(J_x*delta + F_x, 'Inf');
        step_lengths(counter) = alpha;
        det_F_x(counter) = det(J_x);
        fprintf('Iteration %d - f(x) = %2.3e  \n', ...
            counter, norm_F_x_nm(counter));

    end
    
    else % Gradient Descent
        
        disp('********************************************************************************')
        
        delta = - J_x' * F_x;
        delta_vers = delta / norm(delta);
        ia = 1;
        
        while ia < numel(poss_alpha_2)
            
            alpha = poss_alpha_2(ia);
            xnew = x + alpha * delta_vers;
            xnew(xnew<0) = x(xnew<0);
            F_x_new = f_evaluate_mim(rate_constants, xnew, idx_basic_species, Nl, rho, MIM);
            norm_F_x_new = norm(F_x_new);
            theta_x = 0.5 * norm_F_x^2;
            theta_x_new = 0.5 * norm_F_x_new^2;
            
            if theta_x_new <= theta_x + sigma_2 * (-delta_vers)' * (xnew - x)
                is = ia;
                ia = numel(poss_alpha_2);
                FLAG = 0;
                norm_grad(j) = norm(delta);
                j = j+1;
            else
                ia = ia+1;
                alpha = poss_alpha_2(ia);
                is = ia;
            end
        end
        
        % Store some informations
        norm_F_x_nm(counter) = norm_F_x_new;
        est_err(counter) = norm(J_x * delta + F_x, 'Inf');
        step_lengths(counter) = alpha;
        det_F_x(counter) = det(J_x);
        fprintf('Iteration %d - f(x) = %2.3e  ia = %d \n', ...
            counter, norm_F_x_nm(counter), is);        
        
    end
    
    diff_iter = norm(xnew - x);
    
    if counter > max_counter, break; end
end

% 5.c. Store results (for plotting)
x_m1 = xnew;

elapsed_time = toc(time_init) + elapsed_time;

% Restart if convergence hasn't been reached
if (counter == max_counter+1) && (norm_F_x_new > tol)
    z = z+1;
    
else
    z = z+1;
    % Step 6. Store results over run
    ris(ir).x0 = x_0_nm;
    ris(ir).x  = x_m1;
    ris(ir).num_iter = counter;
    ris(ir).step_lengths = step_lengths(1:counter);
    ris(ir).elapsed_time = elapsed_time;
    ris(ir).norm_F = norm_F_x_nm(counter);
    ris(ir).num_trials = z;
    elapsed_time = 0;
    ir = ir+1;
    z = 0;
end


end

struct.ris = ris;

%% Check: how many times dynamic and Newton-GD get to the same result?


%% Norm of F in the equilibrium point

mean_norm_F = 0;

for i = 1:num_run
    mean_norm_F = ris(i).norm_F + mean_norm_F;
end

struct.all.mean_norm_F = mean_norm_F / num_run;

%% Convergence rate

num_trials = 0;

for i=1:num_run
    num_trials = num_trials + ris(i).num_trials;
end

struct.all.effect_conv = num_run/num_trials;

