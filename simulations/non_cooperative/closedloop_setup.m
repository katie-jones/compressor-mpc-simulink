% Setup simulation parameters for controller
% clear all
addpath('../decentralized_common/')

% Choose type of disturbance
% 1: output, 2: input, 3: asymmetric output, 4: asymmetric input, 5: big output

% Choose number of controller iterations
n_controller_iterations = 3;

% Choose filename and directory for saving results
% Plotting function should take care of ensuring no results are overwritten
results_folder = '../results1/';
results_label = 'Non-cooperative';
results_fname = ['noncoop',num2str(n_controller_iterations),'it'];
results_overwrite = 0;
saveplots = 0;



weights;

for n_disturbance=5
% Run MpcSetup script, perform simulation and plot results
MpcSetup;
% sim('decentralized_closedloop');
% makeplots;
end

