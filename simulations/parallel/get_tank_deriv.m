function pdot = get_tank_deriv(x,u_d)

xsize = 5;

% States
P_D = x(end);

x1 = x(1:xsize);
x2 = x(xsize+1:end-1);


[~,~,~,~,~,~,D2,m_out_c] = comp_coeffs();

[SpeedSound,~,~,~,~,~,ud1,ud2] = flow_params();


% Calculate m_out1 from tank of compressor 1

p2 = x1(2);%
m_out1 = get_mass_flow(p2,P_D,ud1,D2,m_out_c);




% Calculate m_out2 from tank of compressor 2

p2 = x2(2);%
m_out2 = get_mass_flow(p2,P_D,ud2,D2,m_out_c);



% Calculate m_out from large tank

[Out_pres,VolumeT,D2,m_out_c] = tank_params();
m_out_tank = get_mass_flow(P_D,Out_pres,u_d,D2,m_out_c);


% m_in to large tank
m_in = m_out1+m_out2;


pdot = SpeedSound * SpeedSound / VolumeT * (m_in - m_out_tank) * 1e-5; % p1


end


% Get mass flow out of a tank given pressures, valve setting, valve
% coefficients
function m_out = get_mass_flow(p,Out_pres,u_out,D2,m_out_c)

dp_sqrt2 = sqrt(abs(p*100 - Out_pres*100)) * sign(p*100 - Out_pres*100);      

M5 = [dp_sqrt2*u_out^3 dp_sqrt2*u_out^2 dp_sqrt2*u_out dp_sqrt2 ...
       u_out^3 u_out^2 u_out 1]';

m_out = D2 * M5 + m_out_c;

end