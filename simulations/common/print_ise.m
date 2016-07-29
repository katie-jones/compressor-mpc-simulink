res_cent = add_ise(res_cent);
res_coop = add_ise(res_coop);
res_ncoop = add_ise(res_ncoop);

for i=1:4
    fprintf(' & %8.2g',[res_cent.yse(i),res_cent.yae(i),res_coop.yse(i),res_coop.yae(i),res_ncoop.yse(i),res_ncoop.yae(i)]);
    fprintf('\n\n');
end
fprintf('\nInputs\n\n')
for i=1:4
    fprintf(' & %8.2g',[res_cent.use(i),res_cent.uae(i),res_coop.use(i),res_coop.uae(i),res_ncoop.use(i),res_ncoop.uae(i)]);
    fprintf('\n\n');
end
