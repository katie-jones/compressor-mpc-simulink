function print_results_info(results_folder,results_label)
if results_folder(end)~='/'
    strcat(results_folder,'/');
end
ofile = fopen([results_folder,results_label,'settings.txt'],'w+');
weight_file = fopen('weights.m','r+');
mpc_file = fopen('const_mpc.m','r+');
tank_file = fopen('../parallel_common/const_tank.m','r+');
barrier_file = fopen('../common/const_barrier.m','r+');
setup_file = fopen(which('MpcSetup.m'),'r+');

find_print(weight_file,ofile,{'UW ','YW '});
find_print(mpc_file,ofile,{'m =','p =','n_delay'});
find_print(tank_file,ofile,{'Vtank'});
find_print(barrier_file,ofile,{'n =','delta'});
find_print(setup_file,ofile,{'x_init_lin =','yref','yss'});

fclose(ofile);
fclose(weight_file);
fclose(mpc_file);
fclose(tank_file);
fclose(barrier_file);
fclose(setup_file);

end


function find_print(in_file, out_file, strings)
tline = fgets(in_file);
while (ischar(tline))
    for i=1:length(strings)
        if strncmp(tline,strings{i},length(strings{i}))
            if ~strcmp(tline(end),sprintf('\n'))
                fprintf(out_file,[tline,'\n']);
            else
                fprintf(out_file,tline);
            end
        end
    end
    tline = fgets(in_file);
end

end