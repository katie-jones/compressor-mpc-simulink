% generate init file for constant matrices:
delete get_matrices.m;
delete get_matAB.m;
disp(' ');
disp(['All auto-generated functions were deleted!']);

file_name = 'get_matrices.m';
file = fopen(file_name, 'w');
fprintf(file, 'function [H,C,M,UB,LB,b,Ain,Ga,Gb,Gc,xinit,upast] = get_matrices(option)\n');
fprintf(file, '%%#eml\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'H = [');
fclose(file);
dlmwrite(file_name, H ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'C = [');
fclose(file);
dlmwrite(file_name, C ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'M = [');
fclose(file);
dlmwrite(file_name, M ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'UB = [');
fclose(file);
dlmwrite(file_name, UB ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'LB = [');
fclose(file);
dlmwrite(file_name, LB ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'b = [');
fclose(file);
dlmwrite(file_name, b ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'Ain = [');
fclose(file);
dlmwrite(file_name, Ain ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'Ga = [');
fclose(file);
dlmwrite(file_name, Ga ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'Gb = [');
fclose(file);
dlmwrite(file_name, Gb ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'Gc = [');
fclose(file);
dlmwrite(file_name, Gc ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'xinit = [');
fclose(file);
dlmwrite(file_name, xinit ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'upast = [');
fclose(file);
dlmwrite(file_name, upast ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

disp(['The file ' file_name ' was created!']);
%%%%%%%%%%%%%%%%%%%%%

file_name = 'get_matAB.m';
file = fopen(file_name, 'w');
fprintf(file, 'function [A,B] = get_matAB(option)\n');
fprintf(file, '%%#eml\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'A = [');
fclose(file);
dlmwrite(file_name, A ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

file = fopen(file_name, 'a');
fprintf(file, 'B = [');
fclose(file);
dlmwrite(file_name, B ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

disp(['The file ' file_name ' was created!']);