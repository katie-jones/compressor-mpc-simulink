% Setup simulation parameters for controller
clear all
addpath('../decentralized_common/')

% Choose type of disturbance
% 1: output, 2: input, 3: asymmetric output, 4: asymmetric input, 5: big output

% Choose number of controller iterations
n_controller_iterations = 3;

% Choose filename and directory for saving results
% Plotting function should take care of ensuring no results are overwritten
results_folder = '../results';
results_label = 'Cooperative';
results_fname = ['coop',num2str(n_controller_iterations),'it'];
results_overwrite = 0;
saveplots = 1;

if ~exist(results_folder,'dir')
    mkdir(results_folder)
end


weights;

% Run MpcSetup script, perform simulation and plot results
for n_disturbance=1:4
    MpcSetup;
    sim('decentralized_closedloop');
    makeplots;
end

