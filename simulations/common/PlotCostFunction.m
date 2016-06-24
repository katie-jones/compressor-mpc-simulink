function fig = PlotCostFunction(results, fname,title_string)
n_res = length(results);
max_value = -1;
for i=1:n_res
    max_value = max(max_value, max(results{i}.(fname)));
end

fig = figure; 
hold on
grid on

for i=1:n_res

inds = results{i}.t >= 50;
plot(results{i}.t(inds),results{i}.(fname)(inds)/max_value);

xlabel('Time [s]')
ylabel('Normalized Penalty term')
title(title_string)

fig=printplot(fig);

end

