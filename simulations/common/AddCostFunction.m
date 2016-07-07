function res = AddCostFunction(res,uwt,ywt)

res.Jtd = (res.u(:,[1,3])-(repmat(res.u(1,[1,3]),size(res.u,1),1).*repmat(res.t<=50,1,2) + ...
    repmat(res.u(end,[1,3]),size(res.u,1),1).*repmat(res.t>50,1,2))).^2; 
res.Jur = (res.u(:,[2,4])-(repmat(res.u(1,[2,4]),size(res.u,1),1).*repmat(res.t<=50,1,2) + ...
    repmat(res.u(end,[2,4]),size(res.u,1),1).*repmat(res.t>50,1,2))).^2; 

d_sd = (res.y(:,1:2)-repmat(res.yref(1:2)',size(res.y,1),1));
res.Jsd = (d_sd.^2.*(d_sd<0));
res.Jp = (res.y(:,4)-repmat(res.yref(4),size(res.y,1),1)).^2;

res.J = [res.Jsd, res.Jp]*ywt + [res.Jtd, res.Jur]*uwt;

end