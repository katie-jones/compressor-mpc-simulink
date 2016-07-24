function figs = PlotResults(results,fname,inds,n_skip)
if nargin < 3 || isempty(inds)
    inds = 1:size(results{1}.(fname),2);
end

if nargin < 4
    n_skip = 5;
end

figs = cell(size(inds));
linestyles={'-','--',':'};

for n=1:length(inds)
    
    fig = figure;
    hold on; grid on;
    
    for i=1:length(results)
        plottime = results{i}.t;
        plotvars = results{i}.(fname)(:,inds(n));
        p = plot(plottime(1:n_skip:end),plotvars(1:n_skip:end));%,'linestyle',linestyles{i});
        p.LineStyle = linestyles{i};
    end
    
    figs{n} = printplot(fig);
    
end
end