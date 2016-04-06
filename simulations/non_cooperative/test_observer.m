%%
if (~exist('M','var'))
    MpcSetup;
end

N = 200;
u_rec = 0.0 *ones(N,1);
Td = 0.1*ones(N,1);
Inflow_opening = 0.405;
Outflow_opening = 0.55;

Ts = 0.05;

u_init = [0.304+Td(1),Inflow_opening,Outflow_opening,u_rec(1)]';

if ~(exist('p_out','var'))
    sim('compressor.mdl')
end

t = 0:Ts:(N-1)*Ts;
y = [SD.signals.values,pout.signals.values(:,1)-pout.signals.values(:,2),PD.signals.values];
yr = interp1(pout.time,y,t)';

u = [Td,u_rec,Td,u_rec]';

n_delay = [0;20];

dxaug = zeros(xtotalsize,N+1);
xaug = zeros(xtotalsize,N+1);
yout = zeros(4,N);
yout(:,1) = yr(:,1);

xaug(:,1) = xinit;



for i=2:length(t)
    dxaug(:,i) = dxaug(:,i) + M*(yr(:,i)-yr(:,i-1) - C*(dxaug(:,i)));
    xaug(1:xsize,i) = xaug(1:xsize,i-1) + dxaug(1:xsize,i);
    xaug(xsize+1:end,i) = dxaug(xsize+1:end,i);
    
    yout(:,i) = C*dxaug(:,i) + yout(:,i-1);

    [A,B,C,H1,H2,Ga1,Ga2,Gb1,Gb2,Gc1,Gc2,fd,Sx,Gd1,Gd2,Sf,Su1,Su2] = get_qp_matrices(xaug(:,i),u(:,i-1));
    
    xlin = dxaug(:,i);
    xlin(xsize+1) = u(2,i-1);
    xlin(xsize+2:end) = 0;
    xlin(xsize+1+n_delay(2)) = u(4,i-1);
    
    dxaug(:,i+1) = B*([u(1,i)-u(1,i-1); u(2,i); u(3,i)-u(3,i-1); u(4,i)]) + A*(dxaug(:,i)-xlin) + fd;

end

% 
% xr = x_out.signals.values;
% for i=2:length(x_out.time)
%     fr(i,:) = get_comp_deriv(xr(i,:)',[0.304+u(min(i-1,N),1),Inflow_opening,Outflow_opening,u(min(i-1,N),2)]);
%     fr2(i,:) = (xr(i,:)-xr(i-1,:))/(x_out.time(i)-x_out.time(i-1));
% end


figure; plot(t,[yr(1,:);yout(1,:)]')
figure; plot(t,[yr(2,:);yout(2,:)]')
figure; plot(t,[yr(3,:);yout(3,:)]')
figure; plot(t,[yr(4,:);yout(4,:)]')
