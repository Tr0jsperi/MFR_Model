function [cp] = heatcp(tp,yelem)
% c this program calculates the heat capacity of a particle from Merrick's
% c correlations.
% c calculates daf coal heat capacity
% c tp particle temperature (k)
% c cp particle heat capacity (cal/g/k)
% c u(i) atomic weights of elements
% c y(i) mass fractions of elements (c,h,n,o,s)
% c rgas gas constant (1.987 cal/gmole/k)
% c rgasj gas constant (8314.3 j/kg/k)
u = [12.0; 1.0; 14.0; 16.0; 32.0];
rgas = 1.987;
% commented out because unused:
% rgasj = 8314.3;
% c calculate mean atomic weight, a
a = 0.0;
for i = 1:5
a = a+(yelem(i)/u(i));
end
a = 1/a;
f1 = 380/tp;
f2 = 1800./tp;
cp = rgas/a*((exp(f1)/((exp(f1)-1)/f1)^2)+2*(exp(f2)/((exp(f2)-1)/f2)^2));
% there used to be a separate function called g1:
% g1=exp(z)/((exp(z)-1)/z)^2;
% and cp was calculated as cp = rgas/a *(g1(f1)+2.*g1(f2));
% this was changed to the version above to remove the requirement for the
% g1 function