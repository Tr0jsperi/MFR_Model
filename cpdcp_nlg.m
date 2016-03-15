function [tmsout,xmout,tpout,tgout,fvolout,fcharout,fcrossout,...
ftarout,fmetout,trateout,mwcharout,yNsiteout,fntout,fncharout,...
fntarout,fnhcnout,fntotout,fgasout,ffgasout,yygasout,...
yfout,waterout,convheatout,dpout]...
=cpdcp_nlg(twallvector,tbnr,texit,timax,yelem,mw1,p0,c0,sigp1,...
mdel,xt,tgc,zv,vpz,press,tg,vg,rhop,dp,swell,omegaw,omegaa,...
rhogas,xwbvector,ugvector,kgvector,...
diffwvector,cpgvector,prgas,emiss,rad,WallX)
% This file is the main function file for the cpdcp-nlg model for the
% Cantera-MATLAB model of the MFR combustion facility at BYU. See commments
% in the code for more information.
% Coded by Andrew Mackrory.
%--------------------------------------------------------------------------
% The cpdcp-nlg devolatilization model was originally written in FORTRAN
% 77. The function of this code is essentially identical to the FORTRAN
% version, but the following points should be noted regarding changes made
% to the code in translation:
%
% 1. In matlab .* is an operator and so all decimal points were preceded
% and followed by at least one numeral to avoid confusion:
% e.g. FORTRAN: 1. MATLAB: 1.0
% e.g. FORTRAN: .1 MATLAB: 0.1
% 2. MATLAB is case sensitive so all code was changed to lower case with a
% few exceptions such as L, L0, and PI. The L's were made uppercase
% to avoid confusion with the numeral 1 and PI was made uppercase to avoid
% confusion with the MATLAB function pi.
% 3. The C13 NMR parameter correlation of Genetti and Fletcher was added
% for use in cases where C13 NMR data are unavailable.
% (see cpdcp_nlg_input.m).
% 4. Single precision variables were the most common variable type used in
% the FORTRAN version of the code. MATLAB uses double precision variables
% almost exclusively resulting in slight differences in results as should
% be expected (rounding error is changed).
% 5. A list of variables with explanations of what most of them are was
% added to the comments.
% 6. Additional comments were added. Original FORTRAN comments can be
% recognized by the FORTRAN comment indicator "c" at the start of each
% comment line. Be sure to read the original file header below.
% 7. Parentheses were added to some expressions to enhance readability.
% 8. Some minor bugs were fixed and the style of some sections of code was
% changed to my own programming style (Andrew Mackrory).
%
% Translation from FORTRAN to the MATLAB m-file format was performed by
% Andrew Mackrory at Brigham Young University in 2007.
% andrew.mackrory@gmail.com
% No liability etc. is assumed for the use of the code as outlined in the
% original FORTRAN file header comments below.
%--------------------------------------------------------------------------
% Original file header comments from the FORTRAN source code (1999):
% c This is the CPDCP-NLG model
% c
% c This model was developed by Sandia National Laboratories under
% c FWP 0709 for the Department of Energy's Pittsburgh Energy
% c Technology Center and the DOE Division of Engineering and Geosciences
% c through the Office of Basic Energy Sciences;
% c and by the University of Utah through funding from
% c the Advanced Combustion Engineering Research Center (ACERC), which
% c is principally sponsored by the National Science Foundation, the
% c State of Utah, and by a consortium of industrial companies.
% c The code will not be formally licensed. Neither the U.S. or the
% c DOE, nor any of their employees, makes any warranty, express or
% c implied, or assumes any legal liability or responsibility for the
% c accuracy, completeness, or usefulness of any information, apparatus,
% c product, or process disclosed, or represents that its use would
% c infringe privately owned rights.
% c
% c The CPD model is intended to solve pyrolysis rate equations
% c based on a percolative bond-breaking scheme. This version includes
% c the flash distillation program to distinguish between tar and
% c metaplast. This program also includes a crosslinking scheme.
% c (January, 1991)
% c
% c Most recent modifications to this model include (a) the nitrogen
% c release model of Perry, (b) the model of Genetti to break the light
% c gas into species based on a correlation, and (c) slight modification
% c to mdel to account for the mass associated with c0.
% c These modifications were made at BYU by Dominic Genetti in his
% c M.S. thesis work (1999) and Steve Perry in his Ph.D. work (1999).
% c
% c This version is coupled with a solver for the particle energy and
% c momentum equations.
% c units are g,K,cm,s,cal
% c
% c A blowing correction to the heat transfer model is employed.
% c
% c Merrick heat capacity correlations are used
% c
%--------------------------------------------------------------------------
% SUBROUTINES
% The following functions are called by this main program or eachother:
%
% at.m
% - function required by lightgas.m to calculate area of triangle
% d.m
% - function required by lightgas.m to calculate distance (for
% interpolation)
% flash.m
% - the flash distillation program that distinguishes between tar and
% metaplast.
% gamln.m
% - a program to calculate the ln of the gamma function
% heatap.m
% - calculates ash heat capacity
% heatcp.m
% - calculates the heat capacity of a D.A.F. coal particle from
% Merrick's correlations
% inverf.m
% - calculates the number of standard deviations from the mean
% corresponding to the area under the standard normal probability curve
% lightgas.m
% - calculates the distribution of light gas species based on a look up
% table of reference data
% perkp.m
% - calculates fractions of tar, and gas
% perks.m
% - the meat of the devolatilization model
% xxx.m
% - function required by lightgas.m for interpolation
% yyy.m
% - function required by lightgas.m for interpolation
%
%--------------------------------------------------------------------------
% The variables used in this script and the related Cantera code that calls
% this function are defined here:
% These definitions do not necessarily apply to variables in the
% subroutines, fac being a good example of same name, different function in
% perkp than in this code.
% LIST OF VARIABLES in this script (alphabetic)
% ab pre-exponential factor (unitless) for labile bridge
% dissociation rate
% ac pre-exponential factor (unitless) for composite rate constant
% acr pre-exponential factor (unitless) for crosslinking rate
% ag pre-exponential factor (unitless) for gas release rate
% aind chemical structure parameter
% aind0 chemical structure parameter (initial value)
% alfa mass of ash per particle (stays constant during pyrolysis)
% alfc mass of D.A.F. portion of particle
% alfc0 initial value of alfc
% alfcp same as alfc, but used in predictor step
% alfw mass of water per particle
% alfwp same as alfw, but used in predictor step
% alpha particle mass (as received, i.e. moist)
% alphap this variable was the same as alpha and was removed from the
% code. (This comment is only here for reference).
% an pre-exponential factor (unitless) for nitrogen release by high
% temperature decomposition (slow)
% ans This variable will be generated by MATLAB when certain commands
% are executed. May be useful if you are debugging.
% ap surface area of spherical particle
% arad pre-exponential factor (unitless) for nitrogen attack by free
% radical (fast)
% ASTMvol dry, ash free ASTM volatile wt% (used for optional C13 NMR
% parameter correlation in cpdcp_nlg_input.m)
% b transfer number for effects of high mass transfer
% blow blowing factor = b/(exp(b)-1) for effects of high mass transfer
% bloww blowing factor for water evaporation (just like blow):
% bloww = bw/(exp(bw)-1)s
% bw transfer number for water evaporation
% c vector for storing empirical correlation coefficients for
% optional C13 NMR parameter correlation in cpdcp_nlg_input.m
% c0 char bridge population (C13 NMR parameter)
% cp specific heat capacity of particle
% cpa specific heat capacity of ash
% cpc specific heat capacity of the D.A.F. portion of the particle
% (char)
% cpg specific heat capacity of gas
% cpw specific heat capacity of water, 1 cal/g-K
% del2 chemical structure parameter
% delhv heat of pyrolysis (cal/g), negative indicates endothermic
% nominally -100.0 cal/g
% delhw heat of vaporization for water (cal/g) at 1 atm, negative
% indicates endothermic
% diffw diffusion coefficient of water in the gas
% dp particle diameter (cm)
% dp0 initial particle diameter (cm)
% dt time step (seconds)
% dtmax maximum time step (seconds)
% dvdt rate of mass loss of char due to volatile release
% dy1 change in y(1) over a single time step - used to determine if a
% change in time step is necessary
% eb0 activation energy (cal) for labile bridge dissociation rate
% ebsig standard deviation of activation energy for labile bridge
% dissociation rate
% ec0 activation energy (cal) for composite rate constant
% ecr activation energy (cal) for crosslinking rate
% eg0 activation energy (cal) for gas release rate
% egsig standard deviation of activation energy for gas release rate
% emiss emissivity of particle and other surfaces in reactor
% en0 activation energy (cal) for nitrogen release by high
% temperature decomposition (slow)
% ensig standard deviation of activation energy for nitrogen release by
% high temperature decomposition (slow)
% erad activation energy (cal) for nitrogen attack by free radical
% (fast)
% F(i,j) view factor for radiation heat transfer
% fchar equals 1 - fvol, so it's what's left of the D.A.F. coal after
% the light gas and tar leaves. NOTE that this includes
% metaplast (fmet) calculated by the flash subroutine. In some
% versions of CPD this is "fsolid" instead of "fchar".
% fcross fraction of original D.A.F. coal that was metaplast and
% crosslinked into the char matrix
% fcrossold same as fcross, but the last value of fcross calculated
% ffgas the fraction of total mass release (i.e. volatiles) that is
% h2o, co2, ch4, co, and other light gases.
% fgas mass fraction of D.A.F. coal evolved as light gas
% fgasold fgasold is the value from the last time the flash subroutine
% was called (see fgasold2). Coded in this main script to
% preserve its value between successive calls of the flash
% subroutine.
% fgasold2 fgasold2 is the value at the previous time step in the flash
% subroutine (see fgasold). Coded in this main script to
% preserve its value between successive calls of the flash
% subroutine.
% fid file identifier for output file
% fmet mass fraction of D.A.F. coal existing as metaplast
% fmetold same as fmet, but the value from earlier calculation/iteration
% fnca0 initial mass fraction of nitrogen in site (aromatic)
% fnchar fraction of original nitrogen remaining in char (and metaplast)
% (see comment above for fchar)
% fnhcn fraction of original nitrogen released as light gas (by
% difference: 1 - fnchar - fntar)
% fnit D.A.F. mass fraction of nitrogen in original coals
% fnt nitrogen content of char and metaplast
% fntar fraction of original nitrogen released as tar
% fntot total fractional release of nitrogen
% fracr fraction to account for reduction of metaplast by crosslinking
% in latest time step
% fstable initial fraction of mw decay with no radical n attack) as
% explained in this quote from Perry et al., "Modeling Nitrogen
% Evolution During Coal Pyrolysis", Energy & Fuels, vol. 14,
% no. 5, 2000 page 1099:
% "It was assumed that the radicals formed during the
% initial 3% of light gas release were stable (i.e.,
% fstable = 0.03). this means that Nsite was assumed to
% remain at the value in the parent coal until the
% molecular weight per cluster had decayed to 97% of the
% coal value. It is not clear whether this empiricism is
% really necessary, although it seems to fit the available
% data for high-rank coals somewhat better than using
% fstable = 0, consistent with the concept of the formation
% of a pool of free radicals before steady state is
% reached."
% In some documentation the variable "fst" is mentioned. fst is
% NOT the same thing as fstable. fst is explained in the MS
% Thesis of Genetti (BYU, April 1999). fstable is explained in
% the PhD Dissertation of Perry (BYU, December 1999).
% The value of 0.03 should be used for fstable.
% ft weight fraction of each tar bin
% ftold ftold is the value from the last time the flash subroutine
% was called (see ftold2). Coded in this main script to
% preserve its value between successive calls of the flash
% subroutine.
% ftold(i) = weight fraction of each tar bin
% ftold2 ftold2 is the value at the previous time step in the flash
% subroutine (see ftold). Coded in this main script to
% preserve its value between successive calls of the flash
% subroutine.
% ftar mass fraction of D.A.F. coal evolved as tar
% ftarold ftarold is the value from the last time the flash subroutine
% was called (see ftarold2). Coded in this main script to
% preserve its value between successive calls of the flash
% subroutine.
% ftarold2 ftarold2 is the value at the previous time step in the flash
% subroutine (see ftarold). Coded in this main script to
% preserve its value between successive calls of the flash
% subroutine.
% ftart wt fraction tar, gas, and char? - this definition from comment
% in perkp subroutine
% fvol mass fraction of D.A.F. coal evolved as volatiles
% volatiles = light gas + tar
% g 980.0 cm/s - gravitational acceleration constant
% g0 chemical structure coefficient
% gasmw presumably the molecular weight of gas given off - calculated
% from chemical structure coefficients
% h convective heat transfer coefficient
% heat heat (rate) required for pyrolysis and water evaporation
% i iteration or time step counter
% ik index used in a for loop
% inside (Logical) - when true, O/C and H/C ratios are inside the bounds
% of the library coals in the lightgas subroutine.
% intar (Logical) - when true, tar molecular weight distribution is
% calculated in subroutine perkp.
% ipred True (or 1) when on the predictor step (logical)
% ip array index for property interpolation
% iv array index for velocity interpolation
% ix array index for temperature interpolation
% j index used in for loops
% kg thermal conductivity of the gas (units: cal/cm/s/c)
% L number of labile bridges - see y
% (upper case to avoid confusion with the number one)
% L0 chemical structure coefficient = p0-c0
% (upper case to avoid confusion with the number ten)
% lib (integer) The number of the library coal (1-12) used in the
% lightgas subroutine. If lib = 13 or 14 this corresponds to
% the two extremes of the H/C vs O/C coalification diagram.
% ma chemical structure coefficient - see code
% machar a char nmr parameter calculated from machar = mwchar-sigp1*mdel
% mb chemical structure coefficient - see code
% mdel average molecular weight per side chain (C13 NMR parameter)
% metold metold is the value from the last time the flash subroutine
% was called (see metold2). Coded in this main script to
% preserve its value between successive calls of the flash
% subroutine.
% metold(i) = mass fraction of coal contained in metaplast of
% mer size i
% metold2 metold2 is the value at the previous time step in the flash
% subroutine (see metold). Coded in this main script to
% preserve its value between successive calls of the flash
% subroutine.
% mt molecular weight of each tar bin
% mw1 average molecular weight per aromatic cluster (includes side
% chains) (C13 NMR parameter)
% mwchar like mw1, but for char. Initially it is set to mw1
% mwcharold like mwchar, but used for... initially set to mwchar
% nmax number of terms in expansion for mol. wt. distribution
% nu Nusselt number for particle
% nv number of particle velocity data points
% nx number of gas temperature data points
% omegaa mass fraction of ash in the parent coal (as received)
% omegaw mass fraction of moisture in the parent coal
% (as received, i.e. including ash)
% output filename of file where results are saved
% p0 ratio of bridges to total attachments (C13 NMR parameter)
% PI stores value of pi (upper case signifies it is different to
% the MATLAB function pi, which returns the value of pi)
% pr Prandtl number of the gas
% press pressure (atm)
% pstar chemical structure coefficient
% qconv rate of convective particle heating
% qrad rate of radiative particle heating
% qbnr rate of radiative particle heating from the burner face
% qexit rate of radiative particle heating from the exhaust
% qwall rate of radiative particle heating from the wall
% r2 radiation view factor calculation parameter
% r3 radiation view factor calculation parameter
% rad reactor radius
% ratecr cross linking rate
% rba chemical structure coefficient
% re Reynolds number of flow around particle
% rg universal gas constant rg = 1.987 cal/gmole.K
% rhog gas density
% rhop initial particle apparent density (g/cm^3). See more detailed
% comment in cpdcp_nlg_input.m
% rtot Total rate of mass loss (volatiles and water)
% rtotp predicted total rate of mass loss (volatiles and water)
% rw water evaporation rate
% rwp same as rw, but in predictor step
% sig chemical structure coefficient = sigp1-1
% siginv chemical structure coefficient
% sigma Stefan-Boltzmann constant (for radiation heat transfer)
% sigp1 sig+1 is the coordination number (number of total attachments)
% small a constant small number for comparisons and making small
% adjustments to variables
% swell swelling factor (dpf/dp0 - 1) dpf = final/max diameter. See
% more detailed comment in cpdcp_nlg_input.m
% tarold tarold is the value from the last time the flash subroutine
% was called (see tarold2). Coded in this main script to
% preserve its value between successive calls of the flash
% subroutine.
% tarold2 tarold2 is the value at the previous time step in the flash
% subroutine (see tarold). Coded in this main script to
% preserve its value between successive calls of the flash
% subroutine.
% tg Gas temperature (K)
% tgc Gas temperatures (data as function of position in reactor) (K)
% tbnr burner temperature (K)
% timax maximum devolatilization time modeled (seconds)
% time time (in seconds). The independent variable of the main
% calculation loop.
% tms time converted to milliseconds for output
% tp Particle temperature
% tpred Predicted particle temperature
% trate Particle heating rate (K/s)
% tratep Predicted particle heating rate (K/s)
% texit exhaust tube temperature (K)
% twall reactor wall temperature (K)
% ug gas viscosity (units: g/(cm.s))
% vg gas velocity
% vp Particle velocity
% vpp Predicted particle velocity
% vpz velocities of particles along z axis (cm/s) (1-d flow of gas
% and particles assumed)
% x Particle position (cm)
% xm Particle position (x) converted to meters for output
% xoc O/C Molar Ratio
% xp Particle position (cm) on the predictor step
% xt z axis locations of gas temperatures (cm) and other properties
% xw0 particle surface water concentration (mole fraction)
% xwb bulk flow water concentration (mole fraction)
% y a four element array:
% y(1) = L labile bridges
% y(2) = del ends
% y(3) = c char links
% y(4) = mass fraction of nitrogen in site (aromatic)
% (initially fnca0)
% yelem 5 element array stores dry, ash-free mass fractions of
% C, H, N, O, S, in that order.
% per cluster) (C13 NMR parameter)
% yf A CPD indicator of the fraction of total light gas that has
% been released. The look up table on light gas composition is
% based on yf. Called Xgas in Genetti's MS thesis - see the
% thesis for more detail.
% yhc H/C Molar Ratio
% ynchar the nitrogen remaining in the char and metaplast
% yntar fraction of original nitrogen released as tar
% yp yp(i) = derivative of y(i) in time
% ypp same as yp, but used on the predictor step
% ypred same as y, but used on the predictor step
% yygas the fractions of light gas release that is h20, co2, ch4, co
% and other light gases.
% z Particle position (cm) in radiation calculations
% zero a constant (0) used for comparisons
% zv z axis locations of particle velocities (cm)
% define constants
g = 980;
delhw = -540.0;
cpw = 1.0;
rg = 1.987; % cal/gmole K
PI = 3.14159;
pr = prgas(1);
nmax = 20; % documentation says this is usually 20
zero = 0.0;
small = 1.0e-7;
dt = 10.0e-6; % initial time steps (seconds)
dtmax = 50.0e-6; % maximum time step (seconds)
% KINETIC PARAMETERS
% In these definitions, ? stands for a "subscript" letter:
% a?'s are pre-exponential factors (unitless)
% e?'s are activation energies (cal)
% (sometimes have a zero at end of variable name)
% e?sig's are standard deviations of activation energies (cal)
% subscript letters:
% b is for the labile bridge dissociation rate
% c is for composite rate constant
%(ec = 0 so ac = rho, the composite rate constant)
% g is for gas release rate
% cr is for crosslinking rate
% rad is for nitrogen attack by free radical (fast)
% n is for nitrogen release by high temperature decomposition (slow)
ab = 2.602e15;
eb0 = 55400;
ebsig = 1800;
%
ac = 0.9;
ec0 = 0;
%
ag = 3.0e15;
eg0 = 69000;
egsig = 8100;
%
acr = 3.0e15;
ecr = 65000;
%
arad = 18.4;
erad = 6000;
%
an = 5.5e7;
en0 = 90000;
ensig = 0;
%
fstable = 0.03; % fstable (initial fraction of mw decay with no radical N
% attack) see the comments in the variable list above for
% more detail. Set equal to 0.03
delhv = -100.0; % heat of pyrolysis (cal/g), negative indicates endothermic
% nominally -100.0 cal/g
% The next ten variables are used in the flash subroutine and are passed in
% and out to store them for the next time flash is called:
metold = zeros(1,nmax);
ftold = zeros(1,nmax);
tarold = zeros(1,nmax);
fgasold = 0.0;
ftarold = 0.0;
metold2 = zeros(1,nmax);
ftold2 = zeros(1,nmax);
tarold2 = zeros(1,nmax);
fgasold2 = 0.0;
ftarold2 = 0.0;
% initialization of variables
yp = zeros(4,1);
ypp = zeros(4,1);
ypred = zeros(4,1);
intar = false;
ftar = 0.0;
fgas = 0.0;
fchar = 1.0;
tms = 0.0;
blow = 1.0;
ix = 1.0;
iv = 1.0;
ip = 1.0;
ipT = 1.0;
x = 0.0;
xm = 0.01*x;
nx = length(tgc);
nv = length(vpz);
lib = uint8(0);
% Initialize output arrays and non-zero initial values
tmsout = zeros(size(xt));
xmout = zeros(size(xt));
tpout = zeros(size(xt));
tgout = zeros(size(xt));
fvolout = zeros(size(xt));
fcharout = zeros(size(xt));
fcrossout = zeros(size(xt));
ftarout = zeros(size(xt));
fmetout = zeros(size(xt));
trateout = zeros(size(xt));
mwcharout = zeros(size(xt));
yNsiteout = zeros(size(xt));
fntout = zeros(size(xt));
fncharout = zeros(size(xt));
fntarout = zeros(size(xt));
fnhcnout = zeros(size(xt));
fntotout = zeros(size(xt));
fgasout = zeros(size(xt));
ffgasout = zeros(5,length(xt));
yygasout = zeros(5,length(xt));
yfout = zeros(size(xt));
tpout(1) = tg;
tgout(1) = tg;
fcharout(1) = 1;
mwcharout(1) = mw1;
yNsiteout(1) = fstable;
fntout(1) = yelem(3);
fncharout(1) = 1;
waterout = zeros(size(xt));
% mass fraction of moisture in ash-containing coal
waterout(1) = omegaw;
convheatout = zeros(size(xt));
% set output "trigger" to 2
outputnow = 2;
% c Save initial char NMR parameters as coal NMR parameters
% c (char parameters calculated independent of those using
% c empirical correlation for mdel)
mwchar = mw1;
mwcharold = mwchar;
machar = mwchar-sigp1*mdel;
% c adjust mdel to correct for c0 (Steve Perry, May 1999)
mdel = mdel/(1.0-c0);
% c empirical correlation to allow a small portion of alpha-carbon to
% c stay with the aromatic cluster
mdel = mdel-7;
% c Now calculate other chemical structure coefficients
L0 = p0 - c0;
mb = 2.0*mdel;
ma = mw1-sigp1*mdel;
sig = sigp1-1;
rba = mb/ma;
fnit = yelem(3);
fnt = fnit;
fnca0 = fnit*mw1/machar;
dp0 = dp;
% c initialize variables
y(1) = L0;
y(2) = 2.0*(1.0-c0-L0);
y(3) = c0;
y(4) = fnca0;
aind0 = L0 + (1.0-c0-L0);
siginv = 1.0/sig;
pstar = 0.5*siginv;
yntar = 0.0;
yf = 0.0;
inside = true;
fcross = 0.0;
% c calculate initial particle velocity
% assumes particle is at the gas temperature at this point
tp = tg;
rhog = rhogas(1);
ug = ugvector(1);
kg = kgvector(1);
cpg = cpgvector(1);
diffw = diffwvector(1);
% also it is assumed that initial particle velocity is equal to initial gas
% velocity:
vp = vg;
fvol = 0.0;
fmet = 0.0;
rtot = 0.0;
xwb = xwbvector(1);
time = 0.0;
% c for now, assume that the apparent density is indicative of the as
% c received coal.
alpha = ((4/3)*PI*(dp/2)^3)*rhop;
alfa = alpha*omegaa;
alfw = alpha*omegaw;
alfc = alpha*(1-omegaa-omegaw);
alfc0 = alfc;
ap = PI*dp^2;
re = rhog*abs(vp-vg)*dp/ug;
nu = 2.0 + 0.6*re^0.5*pr^0.333;
h = nu*kg/dp;
sigma = 1.335e-12; % cal/s cm^2 K^4
heat = 0.0;
% c calculate O/C and H/C ratios for light gas model
xoc = (yelem(4)/16)/(yelem(1)/12);
yhc = yelem(2)/(yelem(1)/12);
% get daf coal heat capacity
[cpc] = heatcp(tp,yelem);
% get ash heat capacity
[cpa] = heatap(tp);
fntar = 0.0;
cp = (alfc*cpc + alfa*cpa + alfw*cpw)/(alpha);
% START OF MAIN CALCULATION LOOP
i = 0; % iteration counter
breakout = false;
while (time < timax)&&(breakout == false)
i = i+1;
% ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c
% c PREDICTOR STEP
% c
% ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
fvolold = fvol;
fcrossold = fcross;
fmetold = fmet;
xp = x + vp*dt;
% calculate gas temperature
if (xp <= xt(nx))
tg = (xt(ix+1)-xp)/(xt(ix+1)-xt(ix))*(tgc(ix)-tgc(ix+1))+tgc(ix+1);
else
% fprintf('\rReached end of gas temperature correlation\r');
breakout = true; % exits the while loop under this condition
end
if (inside == false)
% fprintf('\r!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!\r');
% fprintf('O/C and H/C ratios are outside the\r');
% fprintf('bounds of the library coals.\r');
% fprintf('Estimation of light gas distribution\r');
% fprintf('is based on library coal No. %u \r',lib);
end
if (tg > 4000.)
fprintf('\r!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!\r');
fprintf(' Gas temperature too high---- %d K\r',tg);
end
% c interpolate to get particle velocity from given velocities
if (xp <= zv(nv))
vpp = (zv(iv+1)-xp)/(zv(iv+1)-zv(iv))*(vpz(iv)-vpz(iv+1))...
+vpz(iv+1);
else
% fprintf('Reached end of particle velocity correlation\r');
breakout = true; % exits the while loop under this condition
end
% interpolate to get other properties
if (xp <= zv(nv))
rhog = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(rhogas(ip)-rhogas(ip+1))...
+rhogas(ip);
ug = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(ugvector(ip)-...
ugvector(ip+1))+ugvector(ip);
kg = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(kgvector(ip)-...
kgvector(ip+1))+kgvector(ip);
cpg = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(cpgvector(ip)-...
cpgvector(ip+1))+cpgvector(ip);
diffw = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(diffwvector(ip)-...
diffwvector(ip+1))+diffwvector(ip);
xwb = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(xwbvector(ip)-...
xwbvector(ip+1))+xwbvector(ip);
pr = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(prgas(ip)-...
prgas(ip+1))+prgas(ip);
twall = (WallX(ipT+1)-xp)/(WallX(ipT+1)-WallX(ipT))*...
(twallvector(ipT)-twallvector(ipT+1))+twallvector(ipT);
else
% fprintf('\rReached end of property correlations\r');
breakout = true; % exits the while loop under this condition
end
% Reynolds number set to zero here because as coded, this model
% assumes gas velocity = particle velocity:
% vlag = (vpp-vg);
% re = rhog*abs(vlag)*dp/ug;
re = 0.0;
% c energy equation
% c-- convection
nu = 2 + 0.6*re^0.5*pr^0.333;
b = cpg*(rtot)/(2.0*PI*dp*kg);
if (b >= 1.e-4)
blow = b/(exp(b)-1);
else
blow = 1.0;
end
h = blow*nu*kg/dp;
qconv = h*ap*(tg-tp);
% c-- mass transfer
if (alfw > 0)
bw = (rtot)/(2*PI*dp*diffw*rhog);
if (bw>=1.e-4)
bloww = bw/(exp(bw)-1);
else
bloww = 1.0;
end
else
bloww = 1.0;
end
% CHANGE THIS NEXT SECTION TO SUIT GEOMETRY OF SPECIFIC REACTORS:
% c-- radiation
z = xp; % distance from burner (cm)
% distance from exhaust will be (200-z) because reactor is 2m
% long
% set up areas, etc of radiation enclosure
Area(1) = PI*rad^2; % burner (Area vector = same order as emiss vector)
Area(2) = 2*PI*rad*200; % walls
Area(3) = Area(1); % exhaust
Area(4) = ap; % particle
Temp(1) = tbnr;
Temp(2) = twall;
Temp(3) = texit;
Temp(4) = tp;
% avoid divide by zero:
if z <= 0
F(4,1) = 0.5;
else
F(4,1) = 0.5*(1-(1/(1+(rad/z)^2)^0.5)); %disk to sphere view factor
end
F(4,3) = 0.5*(1-(1/(1+(rad/(200-z))^2)^0.5)); %same as above
F(4,2) = 1 - (F(4,1) + F(4,3)); % by summation rule
% by reciprocity:
F(1,4) = F(4,1)*Area(4)/Area(1);
F(2,4) = F(4,2)*Area(4)/Area(2);
F(3,4) = F(4,3)*Area(4)/Area(3);
qrad = emiss(4)*( F(1,4)*Area(1)*emiss(1)*sigma*Temp(1)^4+...
F(2,4)*Area(2)*emiss(2)*sigma*Temp(2)^4+...
F(3,4)*Area(3)*emiss(3)*sigma*Temp(3)^4-...
Area(4)*sigma*Temp(4)^4); % cal/s
% END OF SECTION TO CHANGE TO SUIT GEOMETRY OF SPECIFIC REACTORS
% THERE IS ANOTHER SECTION TO CHANGE BELOW
% c-- water evaporation rate
% using Antoine vapor pressure correlation
if (alfw > 0)
xw0 = exp(18.3036-3816.44/(tp-46.13))/(760*press);
% In above line, (760*press) was previously just 760 in the FORTRAN
% source code. The pressure was added to the equation to allow for
% pressures other than 1 atm.
xw0 = min(xw0,1.0);
xw0 = max(xw0,0.0);
rwp = bloww*2*rhog*diffw*PI*dp*(xw0-xwb)/(1.0-xw0);
else
rwp = 0;
end
% c-- coal pyrolysis rate
[ypp] = perks(y,ypp,tp,L0,c0,ab,eb0,ebsig,ac,ec0,ag,eg0,egsig,...
rg,fnca0,an,en0,ensig);
% c free radical light gas nitrogen release mechanism
if ((mw1-mwchar)/mw1 > fstable)
ypp(4) = ypp(4) - y(4)*arad*exp(-erad/rg/tp)*...
(mwcharold-mwchar)/mwchar*machar/mwchar/dt;
% This needs some parentheses to avoid ambiguity
end
% c component mass conservation
for j = 1:4
ypred(j) = y(j) + dt*ypp(j);
ypred(j) = max(ypred(j),zero);
end
% c crosslinking rate
fracr = 1.0;
if ((fmetold > small) && (acr > 0.0))
ratecr = acr*exp(-ecr/rg/tp)*fmetold*dt;
fracr = 1.0 - ratecr/fmetold;
fmet = fmetold - ratecr;
fcross = fcrossold + ratecr;
if (fmet < 0.0)
fcross = fcrossold + fmet;
fmet = 0.0;
fracr = 0.0;
end
end
% c get product distribution from ypred array
if(ypred(1) > small)
intar = true;
end
[ftar,ftart,fgas,ft,mt] = perkp(ypred,ftar,intar,ma,rba,c0,...
sig,siginv,nmax,pstar);
intar = false;
gasmw = rba*ma/2.0;
% c flash distillation
if(fgas >= 1.0e-5)
ipred = true;
[ftar,fmet,metold,metold2,ftold,ftold2,...
tarold,tarold2,fgasold,fgasold2,ftarold,ftarold2] = ...
flash(fgas,gasmw,ft,mt,fracr,tp,press,nmax,zero,small,ipred,...
metold,metold2,ftold,ftold2,tarold,tarold2,fgasold,fgasold2,...
ftarold,ftarold2);
elseif (fgas < 1.0e-5)
fmet = ftart;
ftar = 0.0;
end
fvol = fgas + ftar;
fchar = 1.0 - fvol;
dvdt = (fvol - fvolold)/dt*alfc0;
% c done with pyrolysis rate!
rtotp = dvdt + rwp;
heat = dvdt*delhv + rwp*delhw;
% c calculate heat capacity
% get daf coal heat capacity
[cpc] = heatcp(tp,yelem);
% get ash heat capacity
[cpa] = heatap(tp);
cp = (alfc*cpc + alfa*cpa + alfw*cpw)/(alpha);
tratep = (qconv + qrad + heat)/(alpha*cp);
tpred = tp + tratep*dt;
% c component mass conservation
alfwp = alfw - rwp*dt;
alfcp = fchar * alfc0;
alfcp = max(alfcp,0.0);
alfwp = max(alfwp,0.0);
alpha = (alfcp + alfa + alfwp);
omegaa = alfa/alpha;
% c particle swelling
L = ypred(1);
dp = dp0 * (1.0 + swell*(1.0-L/L0));
ap = PI * dp^2;
% c particle density changes during devolatilization
rhop = alpha/((4.0/3.0)*PI*(dp/2)^3);
% ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c
% c CORRECTOR STEP
% c
% ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
x = x + 0.5*(vpp+vp)*dt;
% c interpolate to get gas temperature
if (x <= xt(nx))
tg = (xt(ix+1)-x)/(xt(ix+1)-xt(ix))*(tgc(ix)-tgc(ix+1))+tgc(ix+1);
else
% fprintf('\rReached end of gas temperature correlation\r');
breakout = true;
end
% c interpolate to get particle velocity
if (x <= zv(nv))
vp = (zv(iv+1)-x)/(zv(iv+1)-zv(iv))*(vpz(iv)-vpz(iv+1))+vpz(iv+1);
else
% fprintf('\rReached end of particle velocity correlation\r');
breakout = true;
end
% interpolate to get other properties
if (x <= zv(nv))
rhog = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(rhogas(ip)-rhogas(ip+1))...
+rhogas(ip);
ug = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(ugvector(ip)-...
ugvector(ip+1))+ugvector(ip);
kg = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(kgvector(ip)-...
kgvector(ip+1))+kgvector(ip);
cpg = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(cpgvector(ip)-...
cpgvector(ip+1))+cpgvector(ip);
diffw = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(diffwvector(ip)-...
diffwvector(ip+1))+diffwvector(ip);
xwb = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(xwbvector(ip)-...
xwbvector(ip+1))+xwbvector(ip);
pr = (xt(ip+1)-xp)/(xt(ip+1)-xt(ip))*(prgas(ip)-...
prgas(ip+1))+prgas(ip);
twall = (WallX(ipT+1)-xp)/(WallX(ipT+1)-WallX(ipT))*...
(twallvector(ipT)-twallvector(ipT+1))+twallvector(ipT);
else
% fprintf('\rReached end of property correlations\r');
breakout = true; % exits the while loop under this condition
end
% If you had gas velocities, this is where you would use them:
% Reynolds number set to zero here because as coded, this model has no
% input fo gas velocity. It is assumed that the gas velocity is equal
% to the particle velocity:
% vlag = (vp-vg);
% re = rhog*abs(vlag)*dp/ug;
re = 0.0;
% c energy equation
% c
% c-- convection
nu = 2.0 + 0.6*re^0.5*pr^0.333;
b = cpg*(rtotp)/(2.0*PI*dp*kg);
if (b >= 1.0e-4)
blow = b/(exp(b)-1);
else
blow = 1.0;
end
h = blow*nu*kg/dp;
qconv = h*ap*(tg-tpred);
% c-- mass transfer
if (alfw > 0.0)
bw = (rtotp)/(2.0*PI*dp*diffw*rhog);
if (bw >= 1.0e-4)
bloww = bw/(exp(bw)-1);
else
bloww = 1.0;
end
else
bloww = 1.0;
end
% CHANGE THIS NEXT SECTION TO SUIT GEOMETRY OF SPECIFIC REACTORS:
% c-- radiation
z = x; % distance from burner (cm)
% distance from exhaust will be (200-z) because reactor is 2m
% long
% set up areas, etc of radiation enclosure
Area(1) = PI*rad^2; % burner (Area vector = same order as emiss vector)
Area(2) = 2*PI*rad*200; % walls
Area(3) = Area(1); % exhaust
Area(4) = ap; % particle
Temp(1) = tbnr;
Temp(2) = twall;
Temp(3) = texit;
Temp(4) = tpred;
% avoid divide by zero:
if z <= 0
F(4,1) = 0.5; % limiting value of view factor F(4,1)
else
F(4,1) = 0.5*(1-(1/(1+(rad/z)^2)^0.5)); %disk to sphere view factor
end
F(4,3) = 0.5*(1-(1/(1+(rad/(200-z))^2)^0.5)); %same as above
F(4,2) = 1 - (F(4,1) + F(4,3)); % by summation rule
% by reciprocity:
F(1,4) = F(4,1)*Area(4)/Area(1);
F(2,4) = F(4,2)*Area(4)/Area(2);
F(3,4) = F(4,3)*Area(4)/Area(3);
qrad = emiss(4)*( F(1,4)*Area(1)*emiss(1)*sigma*Temp(1)^4+...
F(2,4)*Area(2)*emiss(2)*sigma*Temp(2)^4+...
F(3,4)*Area(3)*emiss(3)*sigma*Temp(3)^4-...
Area(4)*sigma*Temp(4)^4); % cal/s
% END OF SECTION TO CHANGE TO SUIT GEOMETRY OF SPECIFIC REACTORS
% c-- water evaporation rate
% using Antoine vapor pressure correlation
if (alfw > 0)
xw0 = exp(18.3036-3816.44/(tpred-46.13))/(760*press);
% In above line, (760*press) was previously just 760 in the FORTRAN
% source code. The pressure was added to the equation to allow for
% pressures other than 1 atm.
xw0 = min(xw0,1.0);
xw0 = max(xw0,0.0);
rw = bloww*2.0*rhog*diffw*PI*dp*(xw0-xwb)/(1.0-xw0);
else
rw = 0;
end
% c-- coal pyrolysis rate
[yp] = perks(ypred,yp,tpred,L0,c0,ab,eb0,ebsig,ac,ec0,ag,eg0,egsig,...
rg,fnca0,an,en0,ensig);
% c time step control
if (y(1) > 5.0e-3)
dy1 = dt*0.5*(yp(1)+ypp(1));
else
dy1 = dt*0.5*(yp(3)+ypp(3));
end
if (abs(dy1) < 0.001)
dt = dt*2;
if (dt < dtmax)
% fprintf('\rAt time = %d dt changed to %d\r',time,dt);
end
elseif (abs(dy1) > 0.02)
dt = 0.01/abs(dy1)*dt;
% fprintf('\rAt time = %d dt changed to %d\r',time,dt);
end
dt = min(dt,dtmax);
% c free radical light gas nitrogen release mechanism
if ((mw1-mwchar)/mw1 > fstable)
yp(4)=yp(4)-y(4)*arad*exp(-erad/rg/tp)*...
(mwcharold-mwchar)/mwchar*machar/mwchar/dt;
end
% c component mass conservation
for j = 1:4
y(j) = y(j) + dt*0.5*(yp(j)+ypp(j));
y(j) = max(zero,y(j));
end
% c update current and old mwchar
mwcharold = mwchar;
g0 = 2.0*(1.0-y(1)-c0)-y(2);
mwchar = mw1-g0*mdel*sigp1/2.0;
% c crosslinking rate
fracr = 1.0;
if ((fmetold > small) && (acr > 0.0))
ratecr = acr*exp(-ecr/rg/tpred)*fmetold*dt;
fracr = 1.0-ratecr/fmetold;
fmet = fmetold-ratecr;
fcross = fcrossold+ratecr;
if (fmet<0.0)
fcross = fcrossold + fmet;
fmet = 0.0;
fracr = 0.0;
end
end
% c get product distribution from y array
if(y(1) > small)
intar = true;
end
[ftar,ftart,fgas,ft,mt] = perkp(y,ftar,intar,ma,rba,c0,...
sig,siginv,nmax,pstar);
gasmw = rba*ma/2.0;
% c flash distillation
if (fgas >= small)
ipred = false;
[ftar,fmet,metold,metold2,ftold,ftold2,...
tarold,tarold2,fgasold,fgasold2,ftarold,ftarold2] = ...
flash(fgas,gasmw,ft,mt,fracr,tpred,press,nmax,zero,small,...
ipred,metold,metold2,ftold,ftold2,tarold,tarold2,fgasold,...
fgasold2,ftarold,ftarold2);
elseif (fgas < 1.0e-5)
fmet = ftart;
ftar = 0.0;
end
intar = false;
fvol = fgas + ftar;
fchar = 1.0 - fvol;
dvdt = (fvol-fvolold)/dt*alfc0;
% c done with pyrolysis rate!
rtot = dvdt + rwp;
heat = dvdt*delhv + rwp*delhw;
% c calculate heat capacity
% get daf coal heat capacity
[cpc] = heatcp(tpred,yelem);
% get ash heat capacity
[cpa] = heatap(tpred);
cp = (alfcp*cpc + alfa*cpa + alfw*cpw)/(alpha);
trate = (qconv + qrad + heat)/(alpha*cp);
tp = tp + 0.5*(trate+tratep)*dt;
alfw = alfw - rw*dt;
alfc = fchar*alfc0;
alfc = max(alfc,0.0);
alfw = max(alfw,0.0);
alpha = (alfc + alfa + alfw);
omegaa = alfa/alpha;
% c particle diameter changes due to swelling during devolatilization
L = y(1);
del2 = y(2)/2.0;
aind = del2 + L;
dp = dp0*(1.0 + swell*(1.0-L/L0));
ap = PI*dp^2;
% c nitrogen release calculations
% c calculate tar nitrogen yield
% c (assumes tar is released before light gas and HCN)
yntar = yntar + (ftar-ftarold)*fnt;
ftarold = ftar;
% c nitrogen content of char and metaplast
fnt = y(4)*machar/mwchar;
% c nitrogen remaining in char
% since the way fchar is coded is fchar = char AND metaplast, ynchar is
% the nitrogen remaining in the char and metaplast
ynchar = fchar*fnt;
% c fraction of original nitrogen remaining in char
% (and metaplast) - see comment above
fnchar = ynchar/fnit;
% c fraction of original nitrogen released as tar
fntar = yntar/fnit;
% c fraction of original nitrogen released as light gas (diff.)
fnhcn = 1.0 - fnchar - fntar;
% c total fractional release of nitrogen
fntot = (fnit - fnt*fchar)/fnit;
% c distribute light gas into H2O, CO2, CO, CH4, & other HC's
% c yf is a CPD indicator of the fraction of total light gas
% c that has been released. The look up table on light gas
% c composition is based on yf.
% See Genetti's MS Thesis for more info. In his thesis, yf is called
% Xgas.
yf = 1 - aind/aind0;
[yygas,inside,lib] = lightgas(yf,xoc,yhc);
% c calculate fraction of total mass release that is h2o, co2, ch4,
% c co, and other light gases
% This comment is confusing. Here's the explanation:
% yygas stores the fractions of light gas release that is h20, co2,
% ch4, co and other light gases. ffgas stores the fraction of total
% mass release (i.e. volatiles) that is h2o, co2, ch4, co, and other
% light gases.
ffgas = zeros(5,1);
for ik = 1:5
ffgas(ik) = fgas*yygas(ik);
end
% c particle density changes during devolatilization
rhop = alpha/((4.0/3.0)*PI*(dp/2)^3);
time = time + dt;
tms = time*1000.0;
if (time >= timax)
breakout = true;
end
% c check to see if interpolation indices need update
if (x > xt(ix+1))
ix = ix + 1;
% if (ix >= 50)
% break
% end
end
if (x > zv(iv+1))
iv = iv + 1;
ip = ip + 1;
% if (iv >= 50)
% break
% end
end
if (x > WallX(ipT+1))
ipT = ipT + 1;
end
% store data in output arrays if appropriate
xm = 0.01*x;
yNsite = y(4);
if x >= xt(outputnow)
tmsout(outputnow) = tms;
xmout(outputnow) = xm;
tpout(outputnow) = tp;
tgout(outputnow) = tg;
fvolout(outputnow) = fvol;
fcharout(outputnow) = fchar;
fcrossout(outputnow) = fcross;
ftarout(outputnow) = ftar;
fmetout(outputnow) = fmet;
trateout(outputnow) = trate;
mwcharout(outputnow) = mwchar;
yNsiteout(outputnow) = yNsite;
fntout(outputnow) = fnt;
fncharout(outputnow) = fnchar;
fntarout(outputnow) = fntar;
fnhcnout(outputnow) = fnhcn;
fntotout(outputnow) = fntot;
fgasout(outputnow) = fgas;
ffgasout(:,outputnow) = ffgas;
yygasout(:,outputnow) = yygas;
yfout(outputnow) = yf;
waterout(outputnow) = alfw/alpha;
convheatout(outputnow) = qconv*(4.1868); % Watts/particle
dpout(outputnow) = dp;
outputnow = outputnow + 1;
end
end % end of while loop