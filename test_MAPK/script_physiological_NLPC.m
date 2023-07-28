% This code is run by 'main_script_single_mutation.m'

%% Simulate MIM dynamics in physiological condition
disp('*****     PHYSIOLOGICAL CONDITION     *****')

%% Step 1. Set parameters
x_0 = CMIM.species.std_initial_values;

%% Step 2. Compute x_t
max_counter = 500;
rho = CMIM.matrix.Nl * x_0;
idx_basic_species = find(CMIM.species.std_initial_values>0);

disp('Solving through NLPC... ')
time_init = tic;
aux_phys = f_NLPC_restart(x_0, rate_constants, CMIM.matrix.S, CMIM.matrix.Nl, ...
       rho, idx_basic_species, CMIM.matrix.v, CMIM.matrix.ind_one, max_counter, 0);
aux_phys.elapse_time = toc(time_init);
nlpc_phys = aux_phys;
clear aux_mut tim_init
disp('Done!')

x_eq = nlpc_phys.x(:,1);

%% Step 5. Save

nlpc_phys.x_0 = x_0;
nlpc_phys.x_mut_eq = x_eq;
nlpc_phys.rate_constants = rate_constants;

save(fullfile(folder_results, 'nlpc_phys.mat'), ...
        'nlpc_phys')

%% Clean some memory
clear x_0 x_eq
