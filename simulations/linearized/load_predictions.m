folder = '/home/katie/school/MasterThesis/cpp/build/';

file = fopen([folder,'Su.txt'],'r');
Su2 = fscanf(file, '%f',[m*usize,p*ysize])';
fclose(file);

file = fopen([folder,'Sf.txt'],'r');
Sf2 = fscanf(file, '%f',[orig_xsize,p*ysize])';
fclose(file);

file = fopen([folder,'Sx.txt'],'r');
Sx2 = fscanf(file, '%f',[xsize,p*ysize])';
fclose(file);

Sx2 = [Sx2(:,1:orig_xsize),Sx2(:,orig_xsize+dsize+1:end),Sx2(:,orig_xsize+1:orig_xsize+dsize)];