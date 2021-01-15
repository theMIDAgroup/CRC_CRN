% This code is run by 'main_script_single_mutation.m'

%% Simulate MIM dynamics in mutated condition
fprintf('*****     MUTATION OF %s  perc = %1.2f   ***** \n', protein, perc)

%% Step 1. Define mutation
[~, idx_protein] = ismember(protein, CMIM.species.names);

%% Step 2. Define parameter in the mutated cell
% NOTE: A loss of function (LoF) mutation corresponds to a mutated initial
% condition, while a gain of function (GoF) mutation corresponds to 
% different values of the constant rates.

% 2.a. Initial condition from x_0
switch protein
    
    case lof_mutations % Loss of function mutation
        
    x_0_mutated = x_0; x_0_mutated(idx_protein) = perc*x_0_mutated(idx_protein);
    
    case gof_mutations % Gain of function mutations
        
    x_0_mutated = x_0;
    
    case lof_mutations_type2
    
   [x_0_mutated,  ~] = ...
    f_define_mutated_condition(protein, x_0, rate_constants, CMIM, perc);
            
end

% 2.b. Initial condition from x_eq and rate constants
[x_0_mutated_xe, rate_constants_mutated] = ...
    f_define_mutated_condition(protein, x_eq, rate_constants, CMIM, perc);

% 2.c. Maximum time instant
max_t_mutated = max_t;

%% Step 3. Simulate MIM dynamic 

disp('*****     Dynamic from mutated x_0   *****')
disp('Solving ODE system... ')
[time_mutated, x_t_mutated] = ode15s(@(t_, x_) f_odefun_MIM(...
    t_, x_, rate_constants_mutated, CMIM, 'Sv'), [0 max_t_mutated], x_0_mutated);
disp('Done!')

x_t_mutated = x_t_mutated'; % Size: n_species x n_times
n_times_mutated = numel(time_mutated);

t_eq_mutated = size(x_t_mutated, 2);


disp('*****     Dynamic from mutated x_e    *****')
disp('Solving ODE system... ')
[time_mutated_xe, x_t_mutated_xe] = ode15s(@(t_, x_) f_odefun_MIM(...
    t_, x_, rate_constants_mutated, CMIM, 'Sv'), [0 max_t_mutated], x_0_mutated_xe);
disp('Done!')

x_t_mutated_xe = x_t_mutated_xe'; % Size: n_species x n_times
n_times_mutated_xe = numel(time_mutated_xe);

t_eq_mutated_xe = size(x_t_mutated_xe, 2);

%% Save
ris_mutated.x_0_mut = x_0_mutated;
ris_mutated.x_mut_eq = x_t_mutated(:, end);
ris_mutated.t_eq_mut = t_eq_mutated;

ris_mutated.x_0_xemut = x_0_mutated_xe;
ris_mutated.x_xemut_eq = x_t_mutated_xe(:, end);
ris_mutated.t_eq_xemut = t_eq_mutated_xe;

ris_mutated.rate_constants_mutated = rate_constants_mutated;

switch protein 
    case [lof_mutations, lof_mutations_type2]
     save(fullfile(folder_results, sprintf('results_mutation_%s_perc_%1.1f.mat', protein, perc)), ...
        'ris_mutated')

    case gof_mutations
    save(fullfile(folder_results, sprintf('results_mutation_%s_perc_%1.1f.mat', protein, perc)), ...
        'ris_mutated')
    
end

%% Clean some memory
clear x_0_mutated x_0_mutated_xe x_t_mutated x_t_mutated_xe 