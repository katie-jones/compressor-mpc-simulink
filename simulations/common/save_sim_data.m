%% Save simulation data to mat file

[n_delay,dsize,ucontrolsize,p,m,UWT,YWT] = const_mpc();
[~,Vtank] = const_tank();
[~, ~, ~, ~, ysize] = const_sim();

Info.delay = n_delay(2);
Info.p = p;
Info.m = m;
Info.UWT = diag(UWT(1:2*ucontrolsize,1:2*ucontrolsize));
Info.YWT = diag(YWT(1:ysize,1:ysize));
Info.Vtank = Vtank;

if ~exist('dist_desc','var')
    Info.dist_desc = input('Description of disturbance(s):\n');
else
    Info.disturbance = dist_desc;
end

if ~exist('controltype','var')
    Info.controltype = input('Type of controller:\n');
else
    Info.controltype = controltype;
end

if ~exist('yref','var')
    yref = input('Reference signal: ');
end

fname = input('Filename to save as (no .mat extension):');

save([fname,'.mat'],'Info','PD','pout','SD','u_rec','Td','yref')