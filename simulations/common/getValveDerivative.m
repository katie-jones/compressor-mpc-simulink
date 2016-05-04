function dp2 = getValveDerivative(p1,p2,SpeedSound,V,u,D)
%#eml
dp2 = SpeedSound*SpeedSound/V * 1e-5 * 100/2/sqrt(abs(p1*100-p2*100)) * [u^3, u^2, u, 1] * D(1:4)';

end