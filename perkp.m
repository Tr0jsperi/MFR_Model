function [ftar,ftart,fgas,ft,mt] = perkp(y,ftar,intar,ma,rba,c0,...
sig,siginv,nmax,pstar)
% USAGE:
% [ftar,ftart,fgas,ft,mt] = perkp(y,ftar,intar,ma,rba,c0,...
% sig,siginv,nmax,pstar);
%
% c calculates fractions of tar, and gas from p, sig, L, and c
L = y(1);
del = y(2);
c = y(3);
p = L + c;
if (intar)
if(p > 0.9999)
delfac = 1.0;
else
delfac = del/(1.0-p);
end
a = 1.0 + rba*(L/p + (sig-1.0)/4.0*delfac);
b = (delfac/2.0 - L/p);
% find pstar
pstar0 = pstar;
pinv = siginv+1.0e-4;
if(p >= 0.9999)
pstar = 0.0;
elseif (p >= pinv)
for i = 1:25
f = pstar*(1-pstar)^(sig-1) - p*(1-p)^(sig-1);
fp = (1-pstar)^(sig-1)-pstar*(sig-1)*(1-pstar)^(sig-2);
ppstar = pstar - f/fp;
err = abs(1.0 - ppstar/pstar);
if (err <= 1.0e-4)
break
end
pstar = ppstar;
end
if (err > 1.0e-4)
fprintf('\r!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!\r');
fprintf('pstar did not converge\r');
fprintf(' p = %d \r sig = %d \r pstar = %d \r',p,sig,pstar);
end
else
pstar = p;
end
% c check to see if pstar is in the right range
if ((pstar < 0.0)||((p ~= pstar) && (pstar >= siginv)))
fprintf('\r!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!\r');
fprintf('Error--pstar out of acceptable ranges!\r');
fprintf(' pstar = %d \r',pstar);
fprintf(' p = %d \r sig = %d \r pstar0 = %d \r',p,sig,pstar0);
return
end
sfac = (sig+1.0)/(sig-1.0);
fp = (pstar/p)^(sfac);
kp = (pstar/p)^sfac*(1.0-(sig+1)/2.0*pstar);
% c calculate wt fraction tar, gas, and char
ftart = 2.0*(a*fp + rba*b*kp)/(2.0+rba*(1.0-c0)*(sig+1.0));
else
ftart = 0;
end
tarfac = 1.0 - ftar;
g1 = (2.0*(1.0 - p) - del);
g2 = 2.0*(c - c0);
g0 = g1 + g2;
mgas = rba*ma*g0*(sig+1)/4.0*tarfac;
mtot = ma + rba*ma*(sig+1)/2.0 *(1.0 - c0);
fgas = mgas/mtot;
ft = zeros(1,nmax); % preallocated for speed
mt = zeros(1,nmax); % preallocated for speed
% c calculate tar molecular weight distribution
if (intar)
ftsum = 0.0;
for n = 1:nmax
tn = n*(sig - 1.0) + 2;
xm = n*sig+1.0;
yk = n-1.0;
xm1 = xm+1.0;
% c gamln is the solution to the gamma function in
% the sandia math library
fg1 = gamln(xm1);
if (fg1 <= 1.0e-10)
fgam = 0.0;
else
yk1 = yk + 1.0;
fg2 = gamln(yk1);
xmyk = xm - yk + 1.0;
fg3 = gamln(xmyk);
fgam = exp(fg1 - fg2 - fg3);
end
bnn = (sig+1.0)/(n*sig+1.0)*fgam;
qn = bnn*(p^(n-1))*((1-p)^tn)/n;
% c ft(n) = weight fraction of each tar bin
ft(n) = 2.0*(n*a*qn + rba*b*qn)/(2.0 + rba*(1.0-c0)*(sig+1.0));
ftsum = ftsum + ft(n);
% c check to not divide by zero
if (p <= 1.0e-9)
fac = 0;
else
fac = L/p;
end
tst = 1.0-p;
if (tst <= 1.0e-9)
fac1 = 0.0;
else
fac1 = del/(1.0-p);
end
% c mt(n) = molecular weight of each tar bin
mt(n) = n*ma + (n-1)*rba*ma*fac + tn*rba*ma/4.0*fac1;
end
end