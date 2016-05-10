formatstring='%8.5g ';
fprintf('xinit = [ ');
for i=1:size(x.signals.values,2)
    fprintf(formatstring,x.signals.values(end,i));
end

fprintf(']'';\n\n');

fprintf('yss = [ ');
fprintf(formatstring,pout.signals.values(end,1),SD.signals.values(end,1),pout.signals.values(end,2),SD.signals.values(end,2),PD.signals.values(end,1));
fprintf(']'';\n\n');