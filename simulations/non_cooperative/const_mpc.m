function [n_delay,dsize,ucontrolsize,p,m] = const_mpc()
%#eml

dsize = 2;

ucontrolsize = 2;

n_delay = [0, 40];

p = 50;
m = 1;



end