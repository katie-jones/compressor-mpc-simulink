% Setup simulation parameters for controller
% clear all
addpath('../decentralized_common/')

% Choose type of disturbance
% 1: output, 2: input, 3: asymmetric output, 4: asymmetric input, 5: big output

% Choose number of controller iterations
n_controller_iterations = 10;

% Choose filename and directory for saving results
% Plotting function should take care of ensuring no results are overwritten
results_folder = '../results/cooperative/';
results_fname = [num2str(n_controller_iterations),'it'];
results_overwrite = 0;
saveplots = 1;

if ~exist(results_folder,'dir')
    mkdir(results_folder)
end

% Reference state
xinit = [  0.86733   1.0308  0.17648   394.96        0  0.99889   1.1872  0.17648   394.96        0 ]';
yss = [   1.0308   8.1258   1.1872   8.1258 ]';
yref = [1.03 8.12 1.1872 8.12];

weights;

% Run MpcSetup script, perform simulation and plot results
for n_disturbance=1
    MpcSetup;
%     sim('decentralized_closedloop');
%     makeplots;
end

