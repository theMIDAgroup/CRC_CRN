% This code is run by 'main_script_single_mutation.m'

%% Simulate MIM dynamics in mutated condition
fprintf('*****     MUTATION OF %s  perc = %1.2f   ***** \n', protein, perc)

%% Step 1. Define mutation
[~, idx_protein] = ismember(protein, CMIM.species.names);
name_protein = protein; 
if strcmp(protein, 'Ras')
    name_protein = 'KRAS';
end
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

%% Step 3. Compute equilibrium - NLPC
max_counter = 500;
rho_mutated = CMIM.matrix.Nl * x_0_mutated;
idx_basic_species = find(CMIM.species.std_initial_values>0);

disp('*****     Equilibrium with mutated x_0   *****')
disp('Solving through NLPC... ')
time_init = tic;
aux_mut = f_NLPC_restart(x_0_mutated, rate_constants_mutated, CMIM.matrix.S, CMIM.matrix.Nl, ...
       rho_mutated, idx_basic_species, CMIM.matrix.v, CMIM.matrix.ind_one, max_counter, 0);
aux_mut.elapse_time = toc(time_init);
nlpc_mut = aux_mut;
clear aux_mut tim_init
disp('Done!')

x_mut_eq = nlpc_mut.x(:,1);

disp('*****     Equilibrium with mutated x_e    *****')
disp('Solving through NLPC... ')
aux_mut = f_NLPC_restart(x_0_mutated_xe, rate_constants_mutated, CMIM.matrix.S, CMIM.matrix.Nl, ...
       rho_mutated, idx_basic_species, CMIM.matrix.v, CMIM.matrix.ind_one, max_counter, 0);
aux_mut.elapse_time = toc(time_init);
nlpc_mut = aux_mut;
clear aux_mut tim_init
disp('Done!')

x_xemut_eq = nlpc_mut.x(:,1);

%% Save
a = x_0_mutated;
nlpc_mut.x_0_mut = x_0_mutated;
nlpc_mut.x_mut_eq = x_mut_eq;

nlpc_mut.x_0_xemut = x_0_mutated_xe;
nlpc_mut.x_xemut_eq = x_xemut_eq;
nlpc_mut.rate_constants_mutated = rate_constants_mutated;

switch protein 
    case [lof_mutations, lof_mutations_type2]
     save(fullfile(folder_results, 'mutations', sprintf('nlpc_mut_%s_perc_%1.1f.mat', name_protein, perc)), ...
        'nlpc_mut')

    case gof_mutations
    save(fullfile(folder_results, 'mutations', sprintf('nlpc_mut_%s_perc_%1.1f.mat', name_protein, perc)), ...
        'nlpc_mut')
    
end

%% Clean some memory
clear x_0_mutated x_0_mutated_xe x_t_mutated x_t_mutated_xe 