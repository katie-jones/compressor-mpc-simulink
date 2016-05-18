% Return disturbance times, magnitudes and filename for saving results based on
% predetermined disturbance types.
% Disturbance times are: 30s, 100s
% Disturbance types are:
% 1: Output disturbance (-0.2, -0.2)
% 2: Input disturbance (-0.2, -0.1)
% 3: Asymmetric compressor output disturbance (-0.1, -0.1)
% 4: Asymmetric input disturbance (-0.2, -0.1)
function [tdist,udist1,udist2,fname] = disturbances(n_disturbance)
tdist = [50 100]; % disturbance times

switch n_disturbance
    case 1 % output disturbance
        % in1, out1, outtank, out2
        udist1 = [0 0 -0.1];
        udist2 = 0*udist1;
        fname = 'output_dist';
    case 2 % input disturbance
        udist1 = [-0.1 0 0];
        udist2 = 0*udist1;
        fname = 'input_dist';
    case 3 % asymmetric output disturbance
        udist1 = [0 -0.1 0 0 0];
        udist2 = udist1;
        fname = 'asymm_output';
    case 4 % asymmetric input disturbance
        udist1 = [-0.2 0 0 0 0];
        udist2 = 0.5*udist1;
        fname = 'asymm_input';
    case 5
        udist1 = [0 0 -0.3 0 0];
        udist2 = [0 0 0 0 0];
        fname = 'big_output';
    otherwise
        udist1 = [0 0 -0.1 0 0];
        udist2 = udist1;
        fname = 'small_output';
end

