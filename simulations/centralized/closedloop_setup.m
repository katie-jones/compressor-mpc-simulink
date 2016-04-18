% Setup simulation parameters for controller
clear all

% Choose filename and directory for saving results
% Plotting function should take care of ensuring no results are overwritten
results_folder = '../results/centralized_barrier2/';
results_fname = 'out';
results_overwrite = 0;

if ~exist(results_folder,'dir')
    mkdir(results_folder)
end

% Reference output
yss = [0.2175 0.2175 0 1.12]';
yref = [0.2175 0.2175 0 1.12]';

saveplots = 1;

% Choose type of disturbance
% 1: output, 2: input, 3: asymmetric output, 4: asymmetric input, 5: big output
for n_disturbance=1:6
    % Run MpcSetup script, perform simulation and plot results
    MpcSetup;
    sim('closedloop');
    makeplots;
end
