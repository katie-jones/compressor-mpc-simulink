% generate init file for constant matrices:
delete get_constant_matrices.m;
disp(' ');
disp(['All auto-generated functions were deleted!']);

file_name = 'get_constant_matrices.m';
file = fopen(file_name, 'w');
fprintf(file, 'function [C,M,UB,LB,b,Ain,xinit,upast,deltax] = get_constant_matrices(option)\n');
fprintf(file, '%%#eml\n');
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

file = fopen(file_name, 'a');
fprintf(file, 'deltax = [');
fclose(file);
dlmwrite(file_name, deltax ,'-append', 'newline', 'pc');
file = fopen(file_name, 'a');
fprintf(file, '];\n');
fclose(file);

disp(['The file ' file_name ' was created!']);
