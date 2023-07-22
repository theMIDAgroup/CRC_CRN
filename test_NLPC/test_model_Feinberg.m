clc
clear
close all

addpath(fullfile('..', 'funcs')) 
set(0,'defaulttextInterpreter','latex') 
set(0,'defaultAxesTickLabelInterpreter','latex');

folder_figures = './figures';

%% Step 1. Model definition
eps = 3; % free parameter for the rate constants

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

% 1.4. Rate constants
fein_k =[1, (eps-1)]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 Case 1. Two distint solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 1; % Free parameter for the conservation law.
rho = c;

% Analytic solutions
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

%% Step 4. Plot.
f_fein = figure('units','normalized','outerposition',[0 0 0.7 0.5]);
plot(x0_a, solutions, '.', 'Markersize', 22, 'Color', [0, 100, 0]/256)
ylim([0.5, 2.5])
yticks([1, 2])
yticklabels({'Solution 1', 'Solution 2'})
xlim([0, 0.7])
xlabel('Initial concentration $x_{A,0}$')
set(gca, 'Fontsize', 22)
saveas(f_fein, fullfile(folder_figures, 'ex_multiple_sols.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Case 2. No positive solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = -3; % Free parameter for the conservation law.
rho = c;

% Analytic solutions
x_eq_teo_1 = [0; c];
x_eq_teo_2 = c/eps*[(eps-1); 1]';

%% Step 2. Run NLPC 
x0 = [1; 1];
% x0 = rand(2, 1)*100;
n_alpha_all = [60, 80, 100, 120, 140];
n_runs = numel(n_alpha_all);

max_counter = 500;
proj = 0;

x_eq_fein2 = zeros(2, n_runs);
min_alphas = zeros(1, n_runs);
for ir = 1:n_runs
    n_alpha = n_alpha_all(ir);
    [ris, gs] = f_NLPC_restart_test(x0, fein_k, fein_S, fein_Nl, rho, idx_basic_species, ...
            fein_v, ind_one, max_counter, proj, n_alpha);
    gs_all(ir) = gs;
    x_eq_fein2(:, ir) = gs.sols(:, end);
    min_alphas(ir) = gs.min_alpha;
end

%% Step 3. Collect results in a table
for ir = 1:n_runs
    fprintf('%1.2e  & (%1.2e, %1.2e) \\\\ \n', min_alphas(ir), x_eq_fein2(1, ir), x_eq_fein2(2, ir))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Modified NLPC function to perform the test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [ris, gs] = f_NLPC_restart_test(x_0, rate_constants, S, Nl, rho, idx_basic_species, ...
    v, ind_one, max_counter, proj, n_alpha)

    %% Step 1. Define additional parameters within NLPC
    toll_cond_init_point = 10^17;
    tol = 1e-12;
    poss_alpha_2 = ones(1,n_alpha);
        % User-tuned number of tested alpha in the gradient descent step.
    poss_alpha_2(2) = 0.79;
    for i=3:length(poss_alpha_2)
        poss_alpha_2(i) = poss_alpha_2(i-1) * poss_alpha_2(2);
    end
    poss_alpha = poss_alpha_2(1:20); 
    poss_alpha_2 = poss_alpha_2 * 1e3; 
    sigma = 10^-4;
    sigma_2 = 10^-4;
    rho_newcond = 10^-2;
    FLAG = 0;

    num_try = 2; % The algorithm will not restart

    ir = 1;

    %% Step 2. Compute 'symbolic' form of the Jacobian matrix of v(x; k)
    n_species = size(S, 1);
    jacobian_v = f_compute_analytic_jacobian_v(v, n_species, ind_one);
    
    gs.method = ones(max_counter, 1);
    gs.conv = zeros(max_counter, 1);
    gs.sols = zeros(n_species, max_counter);
    gs.min_alpha = poss_alpha_2(end);

    %% Step 3. Run the algorithm 
    while ir < num_try

        fprintf('@@@@@@@@@@@  RESTART num %d @@@@@@@@@@ \n', ir)
    %%      3.a. Define initial point
    %   We draw x_0 on the given stoichiometric compatibility class until
    %   cond(J_F(x_0)) is small enough
        if ir > 1
        disp('Initialization')
        aux_cond = Inf; 
        while aux_cond > toll_cond_init_point
            x_0 = f_draw_from_ssurf(Nl, rho, idx_basic_species, [-3, 3]);
            aux_cond = cond(f_evaluate_jacobian(rate_constants, x_0, ...
                        S, idx_basic_species, jacobian_v, Nl));
        end
        disp('Let''s start')
        end

    %%      3.b. Initialize storing variables
        norm_F_x_nm = zeros(max_counter, 1);
        step_lengths = zeros(max_counter, 1);
        det_F_x = zeros(max_counter, 1);

    %%      3.c. External while
        x =  x_0; x(x<0) = 0;
        xnew = x; counter = 0; 
        F_x_new = f_evaluate_mim(rate_constants, xnew, idx_basic_species, ... 
                                 Nl, rho, S, v, ind_one);
        norm_F_x_new = norm(F_x_new);

        while (norm_F_x_new > tol) && (counter < max_counter)

            counter = counter + 1;
            x = xnew; F_x = F_x_new; norm_F_x = norm_F_x_new;
            J_x = f_evaluate_jacobian(rate_constants, x, ...
                        S, idx_basic_species, jacobian_v, Nl);
            
            gs.sols(:, counter) = x; % Store the solution at each iteration

    % ************************* Projected Newton ******************************  
            if FLAG == 0 % Newton

                delta = -J_x \ F_x;
                ia = 1;
                while ia <= numel(poss_alpha) 
                    alpha = poss_alpha(ia);
                    xnew = x + alpha * delta;
                    if proj == 0
                        xnew(xnew<0) = x(xnew<0);
                    else
                        xnew(xnew<0) = 0;
                    end
                    F_x_new = f_evaluate_mim(rate_constants, xnew, ...
                        idx_basic_species, Nl, rho, S, v, ind_one);
                    norm_F_x_new = norm(F_x_new);
                    if norm_F_x_new <= sqrt(1-alpha*sigma)*norm_F_x
                        ia = numel(poss_alpha)+1;
                        flag_save(counter) = 0;
                        FLAG = 0;
                        
                        gs.method(counter) = 0;
                        
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
                    step_lengths(counter) = alpha;
                    det_F_x(counter) = det(J_x);
    %                 fprintf('Iteration %d - f(x) = %2.3e  \n', ...
    %                     counter, norm_F_x_nm(counter));

                    % Condition number of the Jacobian Matrix
                    cond_number(counter) = cond(J_x); 
                    rcond_number(counter) = rcond(J_x); 
                    upd_components(counter) = n_species - sum(xnew==x);
                    zeri(counter) = sum(xnew==0);

                end


    % ******************* Projected Gradient Descent **************************
            else 

                flag_save(counter) = 1;

                delta = - J_x' * F_x;  % delta = - gradiente
                delta_vers = delta / norm(delta);
                ia = 1;

                P_x_grad = x + delta_vers; 
                P_x_grad(P_x_grad < 0) = 0;
                diff_P_x = abs(x - P_x_grad); %
                while ia <= numel(poss_alpha_2)
                    alpha = poss_alpha_2(ia);
                    xnew = x + alpha * delta_vers;
                    diff_P_x_M = diff_P_x(xnew >= 0); %
                    diff_P_x_N = diff_P_x(xnew < 0); %

                    if proj == 0
                        xnew(xnew<0) = x(xnew<0);
                    else
                        xnew(xnew<0) = 0;
                    end
                    F_x_new = f_evaluate_mim(rate_constants, xnew, ...
                        idx_basic_species, Nl, rho, S, v, ind_one);
                    norm_F_x_new = norm(F_x_new);
                    theta_x = 0.5 * norm_F_x^2;
                    theta_x_new = 0.5 * norm_F_x_new^2;

                    if ((theta_x_new <= theta_x + sigma_2 * (-delta)' * (xnew - x)) ...
                            && (norm(diff_P_x_M) >= rho_newcond * norm(diff_P_x_N)))
                        ia = numel(poss_alpha_2)+1;
                        FLAG = 0;
                        
                        gs.conv(counter) = 1;
                        
                    else
                        ia = ia+1;
                    end
                end

            % Store some informations
            norm_F_x_nm(counter) = norm_F_x_new;
            step_lengths(counter) = alpha;
            det_F_x(counter) = det(J_x);      

            % Condition number of the Jacobian Matrix
            zeri(counter) = sum(xnew==0);
            cond_number(counter) = cond(J_x); 
            rcond_number(counter) = rcond(J_x); 
            upd_components(counter) = n_species - sum(xnew==x);

            end 
        end

    % 5.c. Store results (for plotting)
    x_res = xnew;

    % Restart if convergence hasn't been reached
    if (counter == max_counter) && (norm_F_x_new > tol)
        ris.cond_number(ir).n = cond_number;
        ris.rcond_number(ir).n = rcond_number;
        ris.upd_components(ir).n = upd_components;
        ris.zeri(ir).n = zeri;
        ris.flag(ir).n = flag_save;
        ir = ir+1; 
        FLAG = 0;
        clear upd_components cond_number zeri flag_save
        if ir == num_try
            warning(['The algorithm did not converge within the allowed number of iterations and restarts. ' ...
                'Results may be inaccurate and the problem may not have a solution.'])
        end
    else
        % Step 6. Store results over run
        ris.x0 = x_0;
        ris.x  = x_res;
        ris.num_iter = counter;
        ris.step_lengths = step_lengths(1:counter);
        ris.norm_F = norm_F_x_nm(counter);
        ris.num_trials = ir;
        ris.upd_components(ir).n = upd_components;
        ris.cond_number(ir).n = cond_number;
        ris.rcond_number(ir).n = rcond_number;
        ris.zeri(ir).n = zeri;
        ris.flag(ir).n = flag_save;
        ir = num_try + 1;
        
        
       gs.method  = gs.method(1:ris.num_iter);
       gs.conv  = gs.conv(1:ris.num_iter);
       gs.sols = gs.sols(:, 1:ris.num_iter);

        clear upd_components cond_number zeri rcond_number flag_save
    end

    end

end




