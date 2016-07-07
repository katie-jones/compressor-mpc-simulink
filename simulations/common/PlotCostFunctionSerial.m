function fig = PlotCostFunctionSerial(results, fname,title_string,flag)
if nargin < 4
    flag = 0;
end

n_res = length(results);
max_value = -1*ones(size(results{1}.(fname),2),1);
for i=1:n_res
    if flag==1 && i~= n_res
        max_value = max(max_value, max(results{i}.(fname)-results{n_res}.(fname)));
    else
%         max_value = max(max_value,max(results{i}.(fname)));
%         for j=1:size(results{1}.(fname),2)
%         max_value(j) = max(max_value(j), max(results{i}.(fname)(:,j)));
%         end
    end
end

for j=1:4
fig1 = figure;
ax1 = axes;
    title(title_string{1,j});

hold on
grid on


fig2 = figure;
ax2 = axes;
    title(title_string{2,j});


hold on
grid on

for i=1:1:n_res
    inds = (results{i}.t >= 50);

    if flag==0
    plot(ax1,results{i}.t,results{i}.Ju(:,j));
        xlabel('Time [s]')
    ylabel('Normalized penalty term')
    plot(ax2,results{i}.t,results{i}.Jy(:,j));
    else
        if i~=n_res
            semilogy(results{i}.t(inds), (results{i}.(fname)(inds)-results{n_res}.(fname)(inds)).*(results{i}.(fname)(inds)>results{n_res}.(fname)(inds)) + 1e-15*(results{i}.(fname)(inds)<=results{n_res}.(fname)(inds)));
        end
    end
    xlabel('Time [s]')
    ylabel('Normalized penalty term')
    
    fig=printplot(fig1);
    
end
end
