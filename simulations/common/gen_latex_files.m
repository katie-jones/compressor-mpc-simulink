fnames = {...
    'parallel_p';
    'parallel_sd';
    'parallel_td';
    'parallel_ur';
    'serial_p1';
    'serial_sd1';
    'serial_td1';
    'serial_ur1';
    'serial_p2';
    'serial_sd2';
    'serial_td2';
    'serial_ur2';
    };

for i=1:length(fnames)
    fig1 = open([fnames{i},'.fig']);
    matlab2tikz(['/home/katie/school/MasterThesis/katie-thesis/report/src/results/figs/',...
        fnames{i},'.tex'],'width','6cm','figurehandle',fig1,'showInfo',false);
end

close all