% This code is run by 'main_script_single_mutation.m'

%% Simulate MIM dynamics in physiological condition
disp('*****     PHYSIOLOGICAL CONDITION     *****')

%% Step 1. Set parameters
x_0 = CMIM.species.std_initial_values;

%% Step 2. Compute x_t
disp('Solving ODE system... ')
[time, x_t] = ode15s(@(t_, x_) f_odefun_MIM(...
    t_, x_, rate_constants, CMIM, 'Sv'), [0 max_t], x_0);

disp('Done!')

x_t = x_t'; % Size: n_species x n_times
n_times = numel(time);

%% Step 3. Compute equilibrium time
disp('Computing equilibrium time...')
t_eq = size(x_t, 2);

%% Step 5. Save
disp('Saving..')
ris_phys.x_0 = x_0;
ris_phys.x_eq = x_t(:, end);
ris_phys.t_eq = t_eq;

save(file_ris_phys, 'ris_phys')

clear time x_t

disp('*****     Done!   *****')