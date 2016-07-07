function figs = PlotResults(results,fname,inds)
if nargin < 3
    inds = 1:size(results{1}.(fname),2);
end

figs = cell(size(inds));

for n=1:length(inds)
    
    fig = figure;
    hold on; grid on;
    
    for i=1:length(results)
        plot(results{i}.t-40,results{i}.(fname)(:,inds(n)))
    end
    
    figs{n} = printplot(fig);
    
end
end