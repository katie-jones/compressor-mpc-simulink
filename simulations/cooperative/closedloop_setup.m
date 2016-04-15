% Setup simulation parameters for controller
clear all
addpath('../decentralized_common/')

% Choose type of disturbance
% 1: output, 2: input, 3: asymmetric output, 4: asymmetric input, 5: big output
n_disturbance = 1;

% Choose number of controller iterations
n_controller_iterations = 2;

% Choose filename and directory for saving results
% Plotting function should take care of ensuring no results are overwritten
results_folder = '../results/cooperative/';
results_fname = [num2str(n_controller_iterations),'2it'];
results_overwrite = 0;

% Reference output
yss = [0.2175 0.2175 0 1.12]';
yref = [0.2175 0.2175 0 1.12]';

% Run MpcSetup script, perform simulation and plot results
MpcSetup;
sim('closedloop');
makeplots;

