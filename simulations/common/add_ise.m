function results = add_ise(results)

results.dy = (results.y - repmat(results.yref',length(results.t),1));
results.yse = sum(results.dy .^ 2)/length(results.t);
results.yae = sum(abs(results.dy))/length(results.t);

inds = results.t >= 50;
tshort = results.t(inds);

results.du = results.u(inds,:)-repmat(results.u(end,:),length(tshort),1);
results.use = sum(results.du .^ 2)/length(tshort);
results.uae = sum(abs(results.du))/length(tshort);

end