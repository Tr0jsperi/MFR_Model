function [y] = gamln(x)
% c this is a program to calculate the ln of the gamma function,
% c taken from Abramowitz, p. 257, 6.1.41
PI = 3.14159;
y = (x-0.5)*log(x)-x+0.5*log(2.0*PI)+1.0/(12.0*x)-1.0/(360.0*x^3)+...
1.0/(1260.0*x^5)-1.0/(1680.0*x^7);
% original FORTRAN code:
% gamln = (x-.5)*alog(x)-x+.5*alog(2.*PI)+1./(12.*x)...
% -1./(360.*x**3)+1./(1260.*x**5)-1./(1680.*x**7)
% Note that in Fortran 77 "alog" is a function that returns the natural
% logarithm. The a at the beginning of alog is there to start the name of
% the function with a letter commonly associated with real rather than
% integer values.