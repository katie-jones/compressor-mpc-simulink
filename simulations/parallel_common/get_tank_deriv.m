function pdot = get_tank_deriv(x,u)

m_in = u(1); % m_in to large tank
u_d = u(2);
Out_pres = u(3);

% States
P_D = x(1);

% Calculate m_out from large tank

VolumeT = const_tank();
[~,~,~,D2,m_out_c] = comp_coeffs();
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

