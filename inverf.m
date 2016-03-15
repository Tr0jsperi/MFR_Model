function [x] = inverf(y)
% c this program calculates the inverse of the area under the normal curve.
% c if y=area(x), then given y, this program will calculate x.
% c a table lookup is performed.
% Calculates the number of standard deviations (x) from the mean
% corresponding to the area under the standard normal probability curve
% from -infinity to x.
xx = [3.4; 3.2; 3.; 2.8; 2.6; 2.4; 2.2; 2.0; 1.8;
1.6; 1.4; 1.2;1.0; 0.8; 0.6; 0.4; 0.2; 0.0];
yy = [0.9997; 0.9993; 0.9987; 0.9974; 0.9953; 0.9918; 0.9861; 0.9772;
0.9641; 0.9452; 0.9192; 0.8849; 0.8413; 0.7881; 0.7257; 0.6554;
0.5793; 0.5];
fac = 1.0;
% c check to see if y is within range
% c if(y<0.0228)then
% c x = -2.0
% c return
if (y < 0.0003)
x = -3.4;
return
elseif (y < 0.5)
yp = 1.0-y;
fac = -1.0;
elseif (y > 0.9997)
x = 3.5;
return
else
yp = y;
end
% c search for range
for i = 18:-1:1
if (yp <= yy(i-1))
x = xx(i) + (yp-yy(i))*(xx(i-1)-xx(i))/(yy(i-1)-yy(i));
x = fac*x;
return
end
end