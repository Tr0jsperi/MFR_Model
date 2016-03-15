function [ftar,fmet,metold,metold2,ftold,ftold2,...
tarold,tarold2,fgasold,fgasold2,ftarold,ftarold2] = ...
flash(fgas,gasmw,ft,mt,fracr,temp,press,nmax,zero,small,ipred,...
metold,metold2,ftold,ftold2,tarold,tarold2,fgasold,fgasold2,...
ftarold,ftarold2)
% c ipred = true only on predictor step, when old values are updated
% Initialize values of some variables:
ntot = nmax + 1;
a = 87058;
b = 299;
G = 0.5903;
x3 = 0.2;
x2 = 0.3;
x = zeros(1,ntot);
y = zeros(1,ntot);
k = zeros(1,ntot);
L = zeros(1,ntot);
v = zeros(1,ntot);
xmw = zeros(1,ntot);
f = zeros(1,ntot);
z = zeros(1,ntot);
pv = zeros(1,ntot);
converged = false;
% c xxxold2 is the value at the previous time step, while
% c xxxold is the value from the last time the subroutine was called
% c metold(i) = mass fraction of coal contained in metaplast of mer size i
% c fracr = fraction to account for reduction of metaplast by crosslinking
% c in latest time step
% c
% c renormalize in tar fragment data to molar basis
% c f(i) = moles of mer i
% c ft(i) = mass of mer i
ftot = 0.0;
for i = 1:nmax
i1=i+1;
xmw(i1) = mt(i);
if (ipred)
ftold2(i) = ftold(i);
metold2(i) = metold(i);
tarold2(i) = tarold(i);
fgasold2 = fgasold;
ftarold2 = ftarold;
end
dif = ft(i) - ftold2(i);
dif = max(dif,zero);
f(i1) = (dif+metold2(i)*fracr)/mt(i);
ftold(i) = ft(i);
ftot = ftot + f(i1);
end
f(1) = (fgas-fgasold2)/gasmw;
f(1) = max(f(1),0.);
fgasold = max(fgas,fgasold2);
xmw(1) = gasmw;
ftot = ftot + f(1);
% c Get mole fraction of components in the feed stream
% c and compute equilibrium constants k(i) from vapor pressure and
% c Raoults law expression
sum = 0.0;
for ii = 1:ntot
sum = sum + f(ii);
pv(ii) = a*exp(-b*xmw(ii)^G/temp);
k(ii) = pv(ii)/press;
if (k(ii) < 0.001)
k(ii) = 0.0;
end
end
if (sum <= 1.0e-8)
return
end
for ii = 1:ntot
z(ii) = f(ii)/sum;
end
% c use the Rachford-Rice formulation for flash distillation
% c x = v/f, first guess
x1 = x3;
% c calculate sum (eq. 11-24, Separation Processes, by King, 1971)
f1 = 0.0;
for ii = 1:ntot
f1 = f1 + z(ii)*(k(ii)-1)/((k(ii)-1)*(x1)+1);
end
% c secant method for convergence
if (x2 == x1)
x2 = x1 + 0.005;
end
for iter = 1:100
% c calculate sum (eq. 11-24, separation processes, by King, 1971)
f2 = 0.0;
for ii = 1:ntot
f2 = f2 + z(ii)*(k(ii)-1)/((k(ii)-1)*(x2)+1);
end
if ((abs(f2) <= small) || (abs(f2-f1) <= small^2))
converged = true;
break
end
x3 = x2 - f2*(x2-x1)/(f2-f1);
if (x3 > 1.0)
x3 = 1.0 - small^2;
end
if (x3 < 0.0)
x3 = 0.0 + small^2;
end
if (x3 == x2)
fprintf('\r Problem---f(v/f) not converged, but x3=x2 \r');
fprintf('x3 = %d \r', x3);
fprintf('x2 = %d \r', x2);
if (x2 >= small)
x3 = x2 - small;
else
x3 = x2 + small;
end
end
if ((x2 <= 1.0e-5) && (x1 <= 1.0e-5))
x2 = 1.0e-7;
converged = true;
break
end
if ((x2 >= 0.9999) && (x1 >= 0.9999))
x2 = 0.9999;
converged = true;
break
end
f1 = f2;
x1 = x2;
% c under-relax solution (especially near the v/f=1 limit
x2 = 0.2*x2 + 0.8*x3;
end
if (converged == false)
% fprintf('\r Convergence not achieved in flash distillation\r');
% fprintf('\r last two guesses for v/f were:\r');
% fprintf('x3 = %d \r', x3);
% fprintf('x2 = %d \r', x2);
ftar = ftarold2;
fmet = 0;
return
end
% c now calculate molecular weight distributions on a
% c light-gas free basis, wt fractions
vtot = ftot * x2;
ltot = ftot - vtot;
vol = vtot/ltot;
sumx = 0.0;
sumy = 0.0;
% xmwtot = 0.0;
ttot = 0.0;
for ii = 2:ntot
i = ii-1;
L(ii) = f(ii)/(1.0 + k(ii)*vol);
v(ii) = f(ii) - L(ii);
x(ii) = L(ii)*xmw(ii);
y(ii) = v(ii)*xmw(ii);
metold(i) = max(x(ii),zero);
tarold(i) = tarold2(i) + y(ii);
% xmwtot = xmwtot + tarold(i)*xmw(ii);
ttot = ttot + tarold(i);
sumx = sumx + x(ii);
sumy = sumy + y(ii);
end
% This code commented out because xmwtot is an unused variable elsewhere
% Same applies to other places where xmwtot appears
% if (ttot > 0.0)
% xmwtot = xmwtot/ttot;
% end
for ii = 2:ntot
if (sumx ~= 0.0)
x(ii) = x(ii)/sumx;
end
if (sumy ~= 0.0)
y(ii) = y(ii)/sumy;
end
end
ftar = ftarold2 + sumy;
% ftarold = ftar; % This should only be executed later in cpdcp_nlg.m
% because of the way ftarold is passed back to the
% calling program. This is different to the FORTRAN
% source code, but the results for fntar are all zero
% for the MATLAB version if this line is not commented
% out.
fmet = sumx;