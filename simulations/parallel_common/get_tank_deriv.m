function pdot = get_tank_deriv(x,u)

% States
P_D = x(1);

% Inputs
m_in = u(1);
Outflow_opening = u(2);

[Out_pres,VolumeT,D2,m_out_c] = const_tank();

SpeedSound = const_flow();


dp_sqrt2 = sqrt(abs(P_D*100 - Out_pres*100)) * sign(P_D*100 - Out_pres*100);      

M5 = [dp_sqrt2*Outflow_opening^3 dp_sqrt2*Outflow_opening^2 dp_sqrt2*Outflow_opening dp_sqrt2 ...
       Outflow_opening^3 Outflow_opening^2 Outflow_opening 1]';

m_out = D2 * M5 + m_out_c;

pdot = SpeedSound * SpeedSound / VolumeT * (m_in - m_out) * 1e-5; % p1
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

