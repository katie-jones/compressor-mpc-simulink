function pdot = get_tank_deriv(x,u)

% States
P_D = x(1);

m_in = u(1);
u_d = u(end);

% Calculate m_out from large tank

[Out_pres,VolumeT,D2,m_out_c] = const_tank();
m_out_tank = get_mass_flow(P_D,Out_pres,u_d,D2,m_out_c);

SpeedSound = const_flow();
pdot = SpeedSound * SpeedSound / VolumeT * (m_in - m_out_tank) * 1e-5; % p1


end


% Get mass flow out of a tank given pressures, valve setting, valve
% coefficients
function m_out = get_mass_flow(p,Out_pres,u_out,D2,m_out_c)
%#eml

dp_sqrt2 = sqrt(abs(p*100 - Out_pres*100)) * sign(p*100 - Out_pres*100);      

M5 = [dp_sqrt2*u_out^3 dp_sqrt2*u_out^2 dp_sqrt2*u_out dp_sqrt2 ...
       u_out^3 u_out^2 u_out 1]';

m_out = D2 * M5 + m_out_c;

end

