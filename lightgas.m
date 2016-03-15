function [yygas,inside,lib] = lightgas(yf,xoc,yhc)
% c ------------------------------------------------------------------
% c --------------light gas distribution model------------------------
% c ------------------------------------------------------------------
% c this program calculates the distribution of light gas species
% c based on a look up table of reference data
% USAGE:
% [yygas,inside,lib] = lightgas(yf,xoc,yhc);
% c *******************************************************************
% c ****************light gas distribution reference library***********
% c *******************************************************************
%
% c this library can be moved outside of submodel as long as
% c it is linked to the light gas sub-model
%
% c xz = the index used in correlation. In the main program it
% c corresponds to variable "yf"
% c f*** = fraction of light gas that is species ***
% c The data is organized in the following order with 12 ordered
% c pairs for each species (ordered with xz)
% c Each table is organized in rows in the following order
% c 1 Lower Kittaning (Chen)
% c 2 Pocahontas #3 (ANL)
% c 3 Upper Freeport (ANL)
% c 4 Pittsburgh (Chen)
% c 5 Lewis Stockton (ANL)
% c 6 Utah Blind Canyon (ANL)
% c 7 Illinois #6 (ANL)
% c 8 Illinois #6 (Chen)
% c 9 Wyodak (ANL)
% c 10 Beulah Zap (ANL)
% c 11 Dietz (Chen)
% c 12 PSOC 1448 (Serio)
% When lib is 13 or 14, the library coal is Rhein Braun or Hongay
% respectively:
% From Genetti's MS Thesis pg 81: "The light gas composition of
% extremely high rank coals was estimated based on the measured light
% gas composition of the pyrolysis products of an anthracite coal,
% known as Hongay, reported by Xu and Tomita. The light gas composition
% of extremely low rank coals was estimated based on data on a lignite,
% known as Rhein Braun, also reported by Xu and Tomita."
% c reference data for xz = yf, the fractional light gas released
xz = [0,0.04,0.11,0.14,0.21,0.27,0.34,0.675,0.9,1.0,0,0;
0,0.161,0.442,0.663,0.777,0.874,0.921,0.967,1.0,0,0,0;
0,0.022,0.20,0.430,0.526,0.64,0.787,0.875,0.927,0.955,1.0,0;
0,0.04,0.12,0.15,0.23,0.29,0.36,0.68,0.9,1.0,0,0;
0,0.018,0.058,0.21,0.417,0.572,0.696,0.778,0.821,0.883,0.932,1.0;
0,0.052,0.144,0.291,0.498,0.639,0.746,0.859,0.925,0.949,0.966,1.0;
0,0.063,0.178,0.33,0.506,0.612,0.706,0.813,0.895,0.94,1.0,0;
0,0.04,0.12,0.15,0.23,0.29,0.36,0.68,0.9,1.0,0,0;
0,0.061,0.146,0.374,0.535,0.622,0.714,0.8,0.883,0.931,0.964,1.0;
0,0.034,0.087,0.179,0.316,0.472,0.585,0.694,0.777,0.872,0.935,1.0;
0,0.04,0.12,0.16,0.25,0.31,0.37,0.68,0.9,1.0,0,0;
0,0.02,0.055,0.17,0.313,0.434,0.546,0.716,0.874,0.935,0.973,1.0];
% c fraction of light gas that is water
fh2o = [0.772,0.772,0.738,0.455,0.371,0.304,0.290,0.273,0.218,0.218,0,0;
0.699,0.632,0.299,0.269,0.247,0.249,0.236,0.225,0.226,0,0,0;
0,0,0.35,0.297,0.301,0.299,0.284,0.291,0.306,0.297,0.283,0;
0.636,0.636,0.646,0.550,0.436,0.320,0.186,0.199,0.195,0.195,0,0;
1.0,0.983,0.754,0.488,0.413,0.385,0.373,0.382,0.377,0.362,0.367,0.348;
0.665,0.636,0.604,0.508,0.435,0.409,0.383,0.362,0.351,0.343,0.342,0.339;
0.763,0.737,0.698,0.572,0.527,0.470,0.438,0.411,0.411,0.396,0.378,0;
0.748,0.748,0.637,0.704,0.490,0.446,0.348,0.268,0.266,0.266,0,0;
0,0,0.385,0.461,0.396,0.369,0.344,0.323,0.292,0.277,0.266,0.257;
0,0,0.197,0.267,0.26,0.333,0.361,0.369,0.346,0.306,0.285,0.267;
0.521,0.521,0.55,0.523,0.511,0.46,0.414,0.388,0.313,0.313,0,0;
0,0,0.291,0.335,0.264,0.271,0.261,0.211,0.171,0.160,0.153,0.149];
% c fraction of light gas that is carbon dioxide
fco2 = [0,0,0,0.174,0.174,0.167,0.129,0.102,0.071,0.071,0,0;
0.259,0.234,0.113,0.086,0.097,0.109,0.116,0.118,0.122,0,0,0;
0.333,0.327,0.070,0.052,0.057,0.06,0.059,0.062,0.066,0.08,0.115,0;
0.194,0.194,0.152,0.117,0.116,0.122,0.081,0.092,0.065,0.065,0,0;
0,0,0,0.122,0.103,0.086,0.083,0.082,0.085,0.086,0.093,0.128;
0.332,0.318,0.165,0.141,0.120,0.108,0.105,0.119,0.120,0.122,0.125,0.130;
0.229,0.221,0.125,0.09,0.07,0.073,0.083,0.133,0.132,0.13,0.147,0;
0.111,0.111,0.142,0.175,0.149,0.155,0.136,0.122,0.133,0.133,0,0;
0.98,0.984,0.55,0.345,0.317,0.285,0.286,0.277,0.273,0.264,0.254,0.255;
0.993,0.989,0.786,0.572,0.519,0.416,0.375,0.345,0.335,0.32,0.303,0.299;
0.363,0.363,0.353,0.325,0.321,0.35,0.318,0.251,0.249,0.249,0,0;
1.0,0.983,0.448,0.179,0.104,0.09,0.104,0.151,0.166,0.160,0.158,0.154];
% c fraction of light gas that is methane
fch4 = [0.203,0.203,0.078,0.160,0.180,0.219,0.258,0.294,0.320,0.320,0,0;
0.041,0.037,0.388,0.389,0.359,0.332,0.323,0.307,0.299,0,0,0;
0.667,0.655,0.42,0.454,0.444,0.419,0.382,0.353,0.331,0.321,0.306,0;
0.055,0.055,0.073,0.088,0.116,0.124,0.170,0.15,0.189,0.189,0,0;
0,0,0.188,0.195,0.234,0.243,0.224,0.21,0.2,0.186,0.177,0.167;
0,0,0.11,0.155,0.176,0.172,0.185,0.173,0.163,0.159,0.156,0.151;
0,0,0.075,0.136,0.159,0.178,0.174,0.157,0.143,0.141,0.132,0;
0.02,0.02,0.026,0.042,0.045,0.049,0.064,0.1,0.128,0.128,0,0;
0,0,0,0.029,0.048,0.067,0.069,0.072,0.069,0.066,0.063,0.061;
0,0,0,0,0.035,0.05,0.061,0.058,0.057,0.053,0.049,0.046;
0.01,0.01,0.011,0.016,0.011,0.021,0.023,0.035,0.06,0.06,0,0;
0,0,0.216,0.262,0.362,0.327,0.307,0.25,0.203,0.189,0.182,0.177];
% c fraction of light gas that is carbon monoxide
fco = [0,0,0.157,0.121,0.141,0.112,0.139,0.085,0.145,0.145,0,0;
0,0,0,0.057,0.097,0.109,0.124,0.15,0.153,0,0,0;
0,0,0,0,0,0.024,0.078,0.097,0.099,0.104,0.099,0;
0.083,0.083,0.038,0.066,0.032,0.168,0.286,0.324,0.313,0.313,0,0;
0,0,0,0,0.055,0.091,0.124,0.131,0.142,0.171,0.168,0.162;
0,0,0,0.028,0.093,0.129,0.142,0.162,0.181,0.191,0.193,0.195;
0,0,0,0.075,0.099,0.122,0.139,0.133,0.148,0.167,0.177,0;
0.101,0.101,0.173,0.054,0.219,0.247,0.335,0.349,0.28,0.280,0,0;
0,0,0.055,0.115,0.151,0.168,0.172,0.2,0.236,0.264,0.287,0.298;
0,0,0,0.133,0.142,0.150,0.15,0.173,0.206,0.265,0.307,0.331;
0.096,0.096,0.066,0.113,0.123,0.13,0.2,0.281,0.334,0.334,0,0;
0,0,0,0.084,0.078,0.115,0.130,0.191,0.262,0.294,0.311,0.322];
% Transpose arrays above (necessary for FORTRAN to MATLAB translation
% without retyping all the data by hand as was done for the 'out'
% array below):
xz = xz';
fh2o = fh2o';
fco2 = fco2';
fch4 = fch4';
fco = fco';
% c ************************************************************
% c **********end of reference library**************************
% c ************************************************************
% c *********************************************************************
% c *******determine the appropriate triangle for interpolation**********
% c *********************************************************************
% c define the equations of line, distance and area
% These equations were moved into the function files yyy.m, xxx.m, d.m,
% and at.m:
% yyy(aa,xd,bb)=aa*xd+bb;
% xxx(aa,yd,bb)=(yd-bb)/aa;
% d(x2,x1,y2,y1)=((x2-x1)^2+(y2-y1)^2)^0.5;
% at(aa,bb,cc)=0.5*bb*cc*(1-((bb^2+cc^2-aa^2)/(2*bb*cc))^2)^0.5;
% c initialize variables
x = xoc;
y = yhc;
ind = yf;
% c look up table of the reference points of the triangular mesh
xx = [0.0177734,0.0203654,0.0659401,0.0773465,0.0893623,0.1077369,...
0.1305803,0.1357569,0.1803479,0.2093441,0.2603201,0.0687];
yy = [0.6717240,0.5810955,0.6550527,0.8088697,0.7575807,0.8506428,...
0.7671163,0.8523190,0.8499221,0.7890888,0.8572938,0.863];
% c look up table for a and b of equations of all triangle sides
a = [-34.965,1.6228,-0.34612,2.3021,1.7337,1.1993,4.3774,...
-4.2685,0.23134,5.0647,1.3746,-3.6565,0.059818,16.459,...
1.6638,-0.05375,0.27897,-2.0979,0.092179,1.3380,3.7559,...
-6.2604,-0.31655];
b = [1.2932,0.54805,0.67788,0.63081,0.54074,0.65041,0.36641,...
1.1390,0.73691,0.30499,0.70255,1.2446,0.84420,-1.3821,...
0.54985,0.85962,0.73069,1.2283,0.83330,0.50899,0.60497,...
1.2931,0.88475];
% c look up table for the three sides that correspond to each triangle
s1 = [1,3,4,7,8,10,12,14,15,18,21,22];
s2 = [2,7,6,5,10,9,14,15,17,20,4,11];
s3 = [3,6,8,9,11,12,13,16,18,19,22,23];
% c loop to find the appropriate triangle for interpolation
m=0;
inside = true;
for i = 1:12
z1=xxx(a(s1(i)),y,b(s1(i)));
z2=xxx(a(s2(i)),y,b(s2(i)));
z3=yyy(a(s3(i)),x,b(s3(i)));
if ((x >= z1) && (x <= z2) && (y <= z3))
m=i;
break
end
end
% this if statement moved outside of above loop for speed:
if((i == 12) && (m == 0.0))
% fprintf('\rOne or both ratios are out of bounds\r');
inside = false;
else
% fprintf('\rTriangle = %d',m);
end
% c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
% c ^^^^^^^^*triangular interpolation^^^^^^^^^^^^
% c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*
if (inside)
% c This interpolation scheme is taken from Zhao et al., 25th Symp.
% c on Comb. (Int.), 1994/pp. 553-560.
% c look up table for points 1,2, and 3 for each triangle
p1 = [2,3,1,3,5,5,7,7,7,10,1,4];
p2 = [1,1,4,5,4,6,6,8,9,9,12,12];
p3 = [3,5,5,7,6,7,8,9,10,11,4,6];
% c calculate the length of each side
ds1 = d(xx(p1(i)),xx(p2(i)),yy(p1(i)),yy(p2(i)));
ds2 = d(xx(p3(i)),xx(p1(i)),yy(p3(i)),yy(p1(i)));
ds3 = d(xx(p3(i)),xx(p2(i)),yy(p3(i)),yy(p2(i)));
ds4 = d(x,xx(p2(i)),y,yy(p2(i)));
ds5 = d(xx(p1(i)),x,yy(p1(i)),y);
ds6 = d(xx(p3(i)),x,yy(p3(i)),y);
% c calculate the area of each triangle used in interpolation scheme
a1 = at(ds1,ds2,ds3);
a2 = at(ds1,ds5,ds4);
a3 = at(ds5,ds2,ds6);
% c calculate s and r, the weighted fraction of two of the points
% c the weighted fraction of other point will be 1-r-s
s = a2/a1;
r = a3/a1;
% c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% c ^^^^^^*calculate light gas distribution^^^^^^^*
% c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% c n is the number of order pairs of data for coals 1-12
% n = [10,9,11,10,12,12,11,10,12,12,10,12];
% line above (n = [...) commented out because it is not used
% c Loop to calculate the light gas composition of each reference
% c coal (point) of triangle. This is accomplished using linear
% c interpolation between points in reference library
% c j specifies the point (1,2, or 3) of the triangle selected above
ygas = zeros(4,3);
for j = 1:3
if (j == 1)
lib=p1(i);
elseif (j == 2)
lib=p2(i);
else
lib=p3(i);
end
% c do loop to find the two xz points that surround ind
for ii = 1:12
if (ind >= xz(ii,lib) && ind <= xz(ii+1,lib))
break
end
end
% c for ygas(k,j)
% c k=1;fh2o
% c k=2;fco2
% c k=3;fch4
% c k=4;fco
% c k=5;other light gases (HC's, H2, parafins, olefins)
% c linear interpolation to find reference values as a function of ind
ygas(1,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*...
fh2o(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-...
xz(ii,lib)))*fh2o(ii+1,lib);
ygas(2,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*...
fco2(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-...
xz(ii,lib)))*fco2(ii+1,lib);
ygas(3,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*...
fch4(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-...
xz(ii,lib)))*fch4(ii+1,lib);
ygas(4,j)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*...
fco(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-...
xz(ii,lib)))*fco(ii+1,lib);
end
end
if (inside)
% c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% c ^^^*calculate gas composition from library coals^^^^*
% c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% c yygas(k) = the fraction of light gas that is k
% c 1=h20, 2=co2, 3=ch4, 4=co
yygas = zeros(5,1);
for k = 1:4
yygas(k) = (1-r-s)*ygas(k,1) + r*ygas(k,2) + s*ygas(k,3);
end
else
% c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% c ^^*estimate composition for coals outside mesh^^^^^^
% c ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
out = [0,0.085,0.12,0.155,0.222,0,0,0.05,0.066,0.175,0,0.222,0.285,0;
0.085,0.12,0.155,0.222,0.285,0.089,0.05,0.175,0.175,0.222,0.175,1,1,0.175;
0.835,0.835,0.835,0.835,0.75,0.75,0.63,0.63,0.69,0,0.55,0,0.75,0;
1.5,1.5,1.5,1.5,1.5,0.835,0.75,0.69,0.835,0.835,0.63,0.75,1.5,0.55;
12,6,8,9,11,4,1,3,7,10,2,10,13,14];
for kk = 1:14
if((x > out(1,kk)) && (x <= out(2,kk)))
if((y >= out(3,kk)) && (y < out(4,kk)))
lib = out(5,kk);
if (lib == 13)
yygas(1) = 0.24;
yygas(2) = 0.37;
yygas(3) = 0.06;
yygas(4) = 0.28;
break
end
if (lib == 14)
yygas(1) = 0.18;
yygas(2) = 0.08;
yygas(3) = 0.37;
yygas(4) = 0.18;
break
end
% c do loop to find the two xz points that surround ind
for ii = 1:12
if((ind >= xz(ii,lib)) && (ind <= xz(ii+1,lib)))
break
end
end
% c linear interpolation to find reference values
% as a function of ind
yygas(1)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*...
fh2o(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-...
xz(ii,lib)))*fh2o(ii+1,lib);
yygas(2)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*...
fco2(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-...
xz(ii,lib)))*fco2(ii+1,lib);
yygas(3)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*...
fch4(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-...
xz(ii,lib)))*fch4(ii+1,lib);
yygas(4)=((ind-xz(ii+1,lib))/(xz(ii,lib)-xz(ii+1,lib)))*...
fco(ii,lib)+((ind-xz(ii,lib))/(xz(ii+1,lib)-...
xz(ii,lib)))*fco(ii+1,lib);
break
end
end
end
% these lines are commented out to speed up execution:
% fprintf('\rLight gas distribution is based on ref. # %d\r',lib);
% fprintf('\rO/C = %d\rH/C = %d',x,y);
end
yygas(5) = 1 - yygas(1) - yygas(2) - yygas(3) - yygas(4);