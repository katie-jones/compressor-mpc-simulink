results_folder = 'results/output_dist_big/';
mkdir(results_folder)
results_fname = '02';
results_title = 'Output Disturbance';

tdist = [30 100]; % disturbance times

% disturbance at time 1
% in1, out1, outtank, in2, out2
udist1 = [0 0 -0.4 0 0];

% disturbance at time 2
% udist2 = udist1;
udist2 = [0 0 0 0 0];

% reference
yref = [0.2 0.2 0 1.12]';
yss = [0.2175 0.2175 0 1.12]';