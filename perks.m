function [yp] = perks(y,yp,t,L0,c0,ab,eb0,ebsig,ac,ec0,ag,eg0,egsig,...
rg,fnca0,an,en0,ensig)
% USAGE:
% [yp] = perks(y,yp,t,L0,c0,ab,eb0,ebsig,ac,ec0,ag,eg0,egsig,...
% rg,fnca0,an,en0,ensig);
% c This subroutine is the meat of the devolatilization model
% y A four element array:
% y(1) = l labile bridges
% y(2) = del ends
% y(3) = c char links
% y(4) = fnca mass fraction of nitrogen in site (aromatic)
% c yp(i) = derivative of y(i) in time
% t = temperature
fx = 0.0;
L = y(1);
del = y(2);
c = y(3);
fnca = y(4);
p = L+c;
g1 = 2.0*(1-p)-del;
g2 = 2.0*(c-c0);
g0 = g1+g2;
% c calculate current activation energy using error function solution
if (c0 < 1.0)
fx = (g0/(1.0-c0))/2.0;
% originally this was fx = g0/(1.0-c0)/2.0
% parentheses added to avoid ambiguity
end
x = inverf(fx);
eg = eg0 + x*egsig;
if (fnca0 > 0)
fx = 1.0 - fnca/fnca0;
end
x = inverf(fx);
en = en0 + x*ensig;
if (L0 > 0.0)
fx = 1.0 - L/L0;
end
x = inverf(fx);
eb = eb0 + x*ebsig;
ec = ec0;
% c calculate rate constants
rt = rg*t;
kb = ab*exp(-eb/rt);
rho = ac*exp(-ec/rt);
kg = ag*exp(-eg/rt);
kn = an*exp(-en/rt);
% c calculate rate of destruction of labile bridges
yp(1) = -(kb)*L;
% c calculate rate of formation of ends (danglers)
yp(2) = 2.0*rho*kb*L/(rho+1.0) - kg*del;
% c calculate rate of formation of char
yp(3) = kb*L/(rho+1.0);
% c calculate rate of high t (slow) nitrogen loss (as hcn)
yp(4) = -kn*fnca;