function fig = PlotCostFunction(results, fname,title_string,flag)
if nargin < 4
    flag = 0;
end

n_res = length(results);
max_value = -1;
for i=1:n_res
    if flag==1 && i~= n_res
        max_value = max(max_value, max(results{i}.(fname)-results{n_res}.(fname)));
    else
        max_value = max(max_value, max(results{i}.(fname)));
    end
end

fig = figure;
if flag==1
end
%     set(gca,'yscale','log')
hold on
grid on

for i=1:1:n_res
    inds = (results{i}.t >= 50);

    if flag==0
        plot(results{i}.t(inds),results{i}.(fname)(inds)/max_value);
    else
        if i~=n_res
            semilogy(results{i}.t(inds), (results{i}.(fname)(inds)-results{n_res}.(fname)(inds)).*(results{i}.(fname)(inds)>results{n_res}.(fname)(inds)) + 1e-15*(results{i}.(fname)(inds)<=results{n_res}.(fname)(inds)));
        end
    end
    xlabel('Time [s]')
    ylabel('Normalized penalty term')
    title(title_string)
    
    fig=printplot(fig);
    
end
