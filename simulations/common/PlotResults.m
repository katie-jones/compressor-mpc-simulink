function figs = PlotResults(results,fname,inds)
if nargin < 3
    inds = 1:size(results{1}.(fname),2);
end

n_skip = 10;

figs = cell(size(inds));

for n=1:length(inds)
    
    fig = figure;
    hold on; grid on;
    
    for i=1:length(results)
        plottime = results{i}.t;
        plotvars = results{i}.(fname)(:,inds(n));
%         plot(results{i}.t-40,results{i}.(fname)(:,inds(n)))
        plot(plottime(1:n_skip:end),plotvars(1:n_skip:end));
    end
    
    figs{n} = printplot(fig);
    
end
end