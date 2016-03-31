% Function fig = printplot(fig,delta)
%
function [fig] = printplot(fig,delta)
set(fig,'units','centimeters');
position = get(fig,'position');
if nargin < 2
    delta = 0.5;
end

set(fig,'position',position);
set(fig,'PaperUnits','centimeters');
set(fig,'PaperSize', [position(3)+2*delta,position(4)+2*delta]);
set(fig,'PaperPosition',[delta, delta, position(3)+delta,position(4)+delta]);
set(fig,'PaperPositionMode','Manual');
