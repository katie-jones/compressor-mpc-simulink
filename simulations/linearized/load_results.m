folder = '/home/katie/school/MasterThesis/cpp/build/';

file = fopen([folder,'output.txt'],'r');
A = fscanf(file, '%f',[1+xsize_orig+ysize+usize,inf])';
fclose(file);

results.t = A(:,1);
results.x = A(:,2:1+xsize_orig);
results.y = A(:,2+xsize_orig:1+xsize_orig+ysize);
results.u = A(:,2+xsize_orig+ysize:end);