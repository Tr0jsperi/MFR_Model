% This is a MATLAB + Cantera model of pulverized coal combustion in BYU's
% MFR combustion research facility. The model consists of four parts:
% 1. A series network of Cantera CSTR's with heat transfer between them for
% ignition of the reactants.
% 2. A Cantera CSTR in a loop that acts as a series of CSTR's without heat
% transfer between them to model the MFR post-ignition. Having heat
% transfer between CSTR's is too computationally expensive for a large
% number of CSTR's
% 3. The CPDCP-NLG coal devolatilization model to provide estimated
% devolatilization products to the CSTR's in 1 and 2 above.
% 4. A char oxidation and CO2 gasification model that begins once
% devolatilization is complete.
%
% (NOTE: This model mostly uses metric units, but the CPD model does
% not. Information is passed in terms of dimensionless quantities, or
% is converted where necessary). The char reactions model uses a
% mixture of units as detailed in the comments and the list of
% variables at the and of the code.
%
% For more information and instructions, see the comments in the code and:
%
% Mackrory, A. J. (2008) A MECHANISTIC INVESTIGATION OF NITROGEN EVOLUTION
% IN PULVERIZED COAL OXY-FUEL COMBUSTION, Ph.D. Dissertation, Brigham
% Young University, Mechanical Engineering Department, December 2008,
% Provo, UT, U.S.A
%
% Coded by Andrew Mackrory using:
% MATLAB Version: 7.4.0 (R2007a)
% Cantera Version: 1.7.1
%---ASSOCIATED SUBROUTINES (m-files)---------------------------------------
% heatap.m - calculates ash heat capacity
% heatcp.m - calculates the heat capacity of a DAF coal particle from
% Merrick's (1983) correlations
% cpdcp_nlg.m - this is the CPD model which uses these subroutines as
% described in the cpdcp_nlg comments:
% at.m
% d.m
% flash.m
% gamln.m
% heatap.m
% heatcp.m
% inverf.m
% lightgas.m
% perkp.m
% perks.m
% xxx.m
% yyy.m
% An input script such as Sub_bit_Air_Input.m is also required. It's name
% may be changed to anything provided the script is called in the code
% below where indicated.
%---CANTERA FUNCTIONS USED IN THE CODE-------------------------------------
% So that these functions may be differentitated from variable names, the
% following is a list of Cantera functions called by this model:
% NOTE: MATLAB is case-sensitive
%
% advance
% air
% cleanup
% cp_mass
% density
% elementIndex
% enthalpy_mass
% GRI30 - creates a gas object using the GRI30.cti mechanism file
% GRI30_B96 - creates a gas object using the GRI30_B96 mechanism file (not
% distributed with Cantera)
% insert
% install
% MassFlowController
% massFraction
% massFractions
% meanMolarMass
% mixDiffCoeffs
% molecularWeights
% moleFraction
% moleFractions
% nAtoms
% nSpecies
% oneatm
% pressure
% Reactor
% ReactorNet
% Reservoir
% setArea
% setInitialVolume
% setMassFlowRate
% setMultiplier
% setThermalResistance
% setValveCoeff
% SKG03 - creates a gas object using the SKG03.cti mechanism file (not
% distributed with Cantera)
% speciesIndex
% speciesName
% temperature
% thermalConductivity
% Valve
% viscosity
% Wall
%---IMPORTANT ASSUMPTIONS--------------------------------------------------
% Key assumptions made in the model are largely based on established
% practices in the literature and include the following:
% Coal particles are entrained (i.e. particle velocity is equal to gas
% velocity).
% All gas products from the coal consist of species in the gas-phase
% kinetic mechanism.
% Natural gas is modeled as 100% CH4
% All nitrogen in the volatiles is in the form of HCN.
% Char consists of C(s) and burns with a shrinking core of constant
% density and constant ash content with CO as the surface product.
% These assumptions were used in deriving the rate constants sourced
% from the literature.
% NO formation from char was not included in the model.
% Sulfur is ignored.
% CO from the char reactions was oxidized to CO2 by the gas phase
% kinetics.
% Fluid mechanics were not modeled as the focus of the model was the
% devolatilization and gas phase kinetics.
% Mixing of burnout oxidizer was assumed to occur in one CSTR
% (i.e. intense mixing).
% The coal particles were represented with one particle diameter based
% on the mean diameter for a Rosin-Rammler Distribution fit to the
% measured particle size distributions.
% More detailed information on these assumptions and references to the
% literature are available in the dissertation referenced above.
%---VARIABLES--------------------------------------------------------------
% An alphabetic list of variables appears at the end of the code.
%---BEGINNING OF CODE------------------------------------------------------
clear; % Clear MATLAB workspace
clc; % and command window
cleanup; % Clear Cantera objects in memory
%---RUN INPUT SCRIPT AND SPECIFY OUTPUT FILE-------------------------------
% Comment out all but one line:
% txt is the appropriate filename extension for the output file
% for use in Microsoft Excel - output is tab-delimited text
case03_Input; output = 'case03_Output.txt';
% Sub_bit_O25_Input; output = 'Sub_bit_O25_Output.txt';
% Sub_bit_O30_Input; output = 'Sub_bit_O30_Output.txt';
% Sub_bit_Air_Opt_Input; output = 'Sub_bit_Air_Opt_Output.txt';
% Sub_bit_O30_Opt_Input; output = 'Sub_bit_O30_Opt_Output.txt';
% Illinois6_Air_Input; output = 'Illinois6_Air_Output.txt';
% Illinois6_O30_Input; output = 'Illinois6_O30_Output.txt';
% Illinois6_O30_0ppm_Input; output = 'Illinois6_O30_0ppm_Output.txt';
% Illinois6_O30_525ppm_Input; output = 'Illinois6_O30_525ppm_Output.txt';
% Pitt8_Air_Input; output = 'Pitt8_Air_Output.txt';
% Pitt8_O30_Input; output = 'Pitt8_O30_Output.txt';
%---OPTIONAL C13 NMR PARAMETER ESTIMATION----------------------------------
% (see notes in input script)
if (mw1 == 0)
yelem = yelem.*100; % convert to percentages for this section of code
% Estimate C13 NMR parameters as follows:
% Declare c (vector of empirical coefficients):
% Estimate parameter using c, yelem, and ASTMvol
% Move on to next parameter and repeat
% Order: mdel, mw1, p0, sigp1, c0
c = [421.957;
-8.64692;
0.0463894;
-8.47272;
1.18173;
1.15366;
-0.0434024;
0.556772;
-0.00654575];
mdel = c(1)+c(2)*yelem(1)+c(3)*yelem(1)^2+c(4)*yelem(2)+...
c(5)*yelem(2)^2+c(6)*yelem(4)+c(7)*yelem(4)^2+c(8)*ASTMvol+...
c(9)*ASTMvol^2;
c = [1301.41;
16.3879;
-0.187493;
-454.773;
51.7109;
-10.072;
0.0760827;
1.36022;
-0.0313561];
mw1 = c(1)+c(2)*yelem(1)+c(3)*yelem(1)^2+c(4)*yelem(2)+...
c(5)*yelem(2)^2+c(6)*yelem(4)+c(7)*yelem(4)^2+c(8)*ASTMvol+...
c(9)*ASTMvol^2;
c = [0.489809;
-0.00981566;
0.000133046;
0.155483;
-0.0243873;
0.00705248;
0.000219163;
-0.0110498;
0.000100939];
p0 = c(1)+c(2)*yelem(1)+c(3)*yelem(1)^2+c(4)*yelem(2)+...
c(5)*yelem(2)^2+c(6)*yelem(4)+c(7)*yelem(4)^2+c(8)*ASTMvol+...
c(9)*ASTMvol^2;
c = [-52.1054;
1.63872;
-0.0107548;
-1.23688;
0.0931937;
-0.165673;
0.00409556;
0.00926097;
-8.26717E-05];
sigp1 = c(1)+c(2)*yelem(1)+c(3)*yelem(1)^2+c(4)*yelem(2)+...
c(5)*yelem(2)^2+c(6)*yelem(4)+c(7)*yelem(4)^2+c(8)*ASTMvol+...
c(9)*ASTMvol^2;
if yelem(1) > 85.9
c0 = min(0.1183*yelem(1)-10.16,0.36);
else
if yelem(4) > 12.5
c0 = min(0.014*yelem(4)-0.175,0.15);
else
c0 = 0;
end
end
yelem = yelem./100; % Undo percent conversion
clear c;
end
%---END OF OPTIONAL C13 NMR PARAMETER ESTIMATION---------------------------
% Setup for calculation of chemical equivalence ratio
% Positive Oxidation States for C,H,N,& O in that order
V_plus = [4 1 0 0];
% Negative Oxidation States for C,H,N,& O in that order
V_minus = [0 0 0 -2];
% Coefficients of element i in species j:
if mechanism == 1
gas = GRI30('Mix');
elseif mechanism == 2
gas = GRI30_B96('Mix');
elseif mechanism == 3
gas = SKG03;
end
CHNOIndex(1) = elementIndex(gas,'C');
CHNOIndex(2) = elementIndex(gas,'H');
CHNOIndex(3) = elementIndex(gas,'N');
CHNOIndex(4) = elementIndex(gas,'O');
% Locations of key species in Cantera gas mixture objects:
H2OIndex = speciesIndex(gas,'H2O');
CO2Index = speciesIndex(gas,'CO2');
CH4Index = speciesIndex(gas,'CH4');
COIndex = speciesIndex(gas,'CO');
C2H2Index = speciesIndex(gas,'C2H2');
HCNIndex = speciesIndex(gas,'HCN');
NOIndex = speciesIndex(gas,'NO');
NO2Index = speciesIndex(gas,'NO2');
O2Index = speciesIndex(gas,'O2');
N2Index = speciesIndex(gas,'N2');
for j = 1:nSpecies(gas)
for i = 1:4
aij(j,i) = nAtoms(gas,j,CHNOIndex(i));
end
end
clear gas
press = P/oneatm; % Note: This does affect the pressure used in the CPD
% sub-model too.
Area = 0.25*pi*d^2; % Cross section area of reactor tube (m^2)
for i = 1:number_reactors
Volume(i) = Area * Length; % volume of each CSTR (m^3)
end
% convert coal flow rate to kg/s
COAL = COAL * (1/3600);
% Convert other reactant input to molar flow rates
m_dot(1:number_reactors) = (NG_in + Air_in + O2_in + CO2_in + N2_in)/3600;
% kg/s mass flow - gas only (coal volatiles added later in code)
NG = NG_in/(12.011+4*1.0079); % moles/hr Natural Gas (assumed 100% CH4)
O2 = (O2_in/(2*15.999) + (1/(1+3.76))*(Air_in/28.851));
% moles/hr O2 from air and bottle
CO2 = (CO2_in/(12.011+2*15.999)); % moles/hr bottled CO2
NO = (NO_doping/1000000)*CO2; % moles/hr NO in the CO2
N2 = (N2_in/(2*14.007) + (3.76/(1+3.76))*(Air_in/28.851));
% moles/hr N2 from air and bottled sources
composition = ['CH4:',num2str(NG),',O2:',num2str(O2),...
',CO2:',num2str(CO2),',N2:',num2str(N2),...
',NO:',num2str(NO)];
particles = (COAL*1000)/(rhop*(4/3)*pi*(dp/2)^3);
% number of coal particles per second
%---BEGINNING OF IGNITION PART OF CODE-------------------------------------
% Create the Gas objects using the selected gas-phase mechanism:
for i = 1:number_reactors
if mechanism == 1
gas(i) = GRI30('Mix');
elseif mechanism == 2
gas(i) = GRI30_B96('Mix');
elseif mechanism ==3
gas(i) = SKG03;
end
set(gas(i),'T',T1,'P',P,'X',composition);
if mechanism == 3
% set multiplier for Thermal NOx Mechanism
setMultiplier(gas(i),74,thermal);
setMultiplier(gas(i),73,thermal);
setMultiplier(gas(i),72,thermal);
% set multiplier for Prompt NOx Fenimore Mechanism
setMultiplier(gas(i),507,prompt);
setMultiplier(gas(i),503,prompt);
else
% set multiplier for Thermal NOx Mechanism
setMultiplier(gas(i),178,thermal);
setMultiplier(gas(i),179,thermal);
setMultiplier(gas(i),180,thermal);
% set multiplier for Prompt NOx Fenimore Mechanism
setMultiplier(gas(i),239,prompt);
setMultiplier(gas(i),240,prompt);
end
end
% dummy gas used for property evaluation at film temperatures
if mechanism == 1
dummygas = GRI30('Mix');
elseif mechanism == 2
dummygas = GRI30_B96('Mix');
elseif mechanism == 3
dummygas = SKG03;
end
% Create upstream reservoirs that will supply the CSTR's with the products
% from the previous CSTR after each iteration (dt).
for i = 1:number_reactors
upstream(i) = Reservoir(gas(i));
end
% Now set the gases to the initial temperature of the CSTR's, and create
% the reactor objects.
% Set their volumes. In this model, the reactor volume is fixed, and
% pressure is mantained by a valve at the outlet of each CSTR.
for i = 1:number_reactors
set(gas(i),'T',TR1,'P',P);
cstr(i) = Reactor(gas(i));
setInitialVolume(cstr(i),Volume(i));
end
% Create a reservoir to represent the environment, and initialize its
% temperature.
ambient = air;
set(ambient,'T',300,'P',P);
env = Reservoir(ambient);
% Create heat-conducting walls between the CSTR's and the
% environment. Set their area, and overall heat transfer
% coefficients.
for i = 1:number_reactors
w(i) = Wall;
install(w(i),cstr(i),env);
setArea(w(i),pi*d*Length);
R = log(0.180/(d/2))/(2*pi*Length*k_wall);
setThermalResistance(w(i),R);
end
% Create heat-conducting walls between the CSTR's themselves.
for i = 1:(number_reactors-1)
gw(i) = Wall;
install(gw(i),cstr(i),cstr(i+1));
setArea(gw(i),Area);
k_gas = 0.5*(thermalConductivity(gas(i))+thermalConductivity(gas(i+1)));
R = Length/(k_gas*Area);
setThermalResistance(gw(i),R);
end
% Connect the upstream reservoirs to the CSTR's with mass flow
% controllers (constant mdot at each one). Set the mass flow rates.
for i = 1:number_reactors
mfc(i) = MassFlowController;
install(mfc(i),upstream(i), cstr(i));
setMassFlowRate(mfc(i),m_dot(i));
end
% Now create downstream reservoirs to exhaust into.
for i = 1:number_reactors
exhaust(i) = air;
set(exhaust(i),'T',300,'P',P);
downstream(i) = Reservoir(exhaust(i));
end
% Connect the CSTR's to the downstream reservoirs with valves, and
% set the coefficient sufficiently large to keep the reactor pressures
% close to the downstream pressure of exhaust.
for i = 1:number_reactors
v(i) = Valve;
install(v(i),cstr(i),downstream(i));
setValveCoeff(v(i), 1.0);
end
% create the network
for i = 1:number_reactors
network_cell_array(i) = {cstr(i)};
end
network = ReactorNet(network_cell_array);
% create vector of locations for variables of interest
position(1) = 0;
for i = 2:number_reactors+1
position(i) = Length+position(i-1);
end
% create arrays for estimation of volatiles species
Cgoal = zeros(1,number_reactors+1);
Hgoal = zeros(1,number_reactors+1);
C_nlg = Cgoal;
H_nlg = Hgoal;
Cdiff = Cgoal;
Hdiff = Hgoal;
% initialize figure for iteration control
figure; subplot(2,1,1); hold on;
title('Temperature Profile - use to check for correct model function');
ylabel('Gas Temperature (K)');
xlabel('Axial Position of CSTR''s (m)');
subplot(2,1,2)
text(0,1,'Choose "Yes" to keep iterating.');
text(0,0.75,'Iterate until max change in T is very small (<1e-5),');
text(0,0.5,'...then select No to continue calculations.');
text(0,0.25,'Select Cancel to end program if T is too low (no ignition)');
axis off;
subplot(2,1,1);
% now integrate in time to a steady state (in a while loop)
tme = 0.0;
n = 0;
iterate = true;
ButtonName = 'Yes';
CPU_time = 0;
delta_max = 1;
while iterate
n = n + 1; % n counts the iterations
tme = tme + dt;
t0 = cputime;
flagtime = 0;
errorcount = 0;
while flagtime < tme
try 
advance(network, tme);
flagtime = tme;
catch exception
errorcount = errorcount + 1;
if errorcount > 10
rethrow(exception);
end
end
end
CPU_time = CPU_time + (cputime - t0);
CPU_time_per_step = CPU_time/n;
% get variables for CPD input
tgas(1) = temperature(upstream(1));
for i = 2:number_reactors+1
tgas(i) = temperature(gas(i-1));
end
if n == 1
tp = tgas;
end
tfilm = 0.5*(tp + tgas);
velocity(1) = m_dot(1)/(density(upstream(1))*Area);
rhogas(1) = density(upstream(1));
set(dummygas,'T',tfilm(1),'P',P,'X',composition);
xwbvector(1) = moleFraction(dummygas,'H2O');
ugvector(1) = viscosity(dummygas);
kgvector(1) = thermalConductivity(dummygas);
DiffCoeffs = mixDiffCoeffs(dummygas);
diffwvector(1) = DiffCoeffs(H2OIndex); % Mixture-averaged diffusion
% coefficient (m^2/s)
clear DiffCoeffs;
cpgvector(1) = cp_mass(dummygas);
prgas(1) = ugvector(1)*cpgvector(1)/kgvector(1);
for i = 2:number_reactors+1
velocity(i) = m_dot(i-1)/(density(gas(i-1))*Area);
rhogas(i) = density(gas(i-1));
xwbvector(i) = moleFraction(gas(i-1),'H2O');
set(dummygas,'T',tfilm(i),'P',P,'Y',massFractions(gas(i-1)));
ugvector(i) = viscosity(dummygas);
kgvector(i) = thermalConductivity(dummygas);
DiffCoeffs = mixDiffCoeffs(dummygas);
diffwvector(i) = DiffCoeffs(H2OIndex); % Mixture-averaged diffusion
% coefficient (m^2/s)
clear DiffCoeffs;
cpgvector(i) = cp_mass(dummygas);
prgas(i) = ugvector(i)*cpgvector(i)/kgvector(i);
end
% call CPD model here (unit conversions are in this function call)
[tms,xm,tp,tg,fvol,fchar,fcross,ftar,fmet,trate,mwchar,yNsite,...
fnt,fnchar,fntar,fnhcn,fntot,fgas,ffgas,yygas,yf,water,convheat,...
dpout] = cpdcp_nlg(twallvector,tbnr,texit,timax,yelem,mw1,p0,c0,...
sigp1,mdel,position*100,tgas,position*100,velocity*100,press,...
tgas(1),velocity(1)*100,rhop,dp,swell,omegaw,omegaa,rhogas/1000,...
xwbvector,ugvector*10,kgvector/418.4,diffwvector*0.01^2,...
cpgvector*2.38846e-4,prgas,emiss,d*50,WallX*100);
y(n,1) = temperature(upstream(1));
for i = 1:number_reactors
y(n,i+1) = temperature(gas(i));
end
if n > 1
delta_max = max(abs((y(n,:)-y(n-1,:))));
plot(position,y(n,:),'bo');
ButtonName = questdlg(['Iterate Again? Max Temp Change = ',...
num2str(delta_max),' K ',...
'(CPU time so far: ',num2str(CPU_time),', = ',...
num2str(CPU_time_per_step),' seconds per step'], ...
'ITERATION CONTROL');
pause(0.1);
end
switch ButtonName,
case 'Yes',
% adjust the thermal conductivity between CSTR's
for i = 1:(number_reactors-1)
k_gas = 0.5*(thermalConductivity(gas(i))+...
thermalConductivity(gas(i+1)));
R = Length/(k_gas*Area);
setThermalResistance(gw(i),R);
end
% calculate mass flows of existing gases + volatiles + water
% (this is evaporated water, not light gas water)
for i = 2:number_reactors
m_dot(i) = m_dot(i-1) + ...
(fvol(i)-fvol(i-1))*COAL*(1-omegaa-omegaw)-...
COAL*(water(i)-water(i-1));
setMassFlowRate(mfc(i),m_dot(i));
end
i = number_reactors+1;
m_dot(i) = m_dot(i-1) + ...
(fvol(i)-fvol(i-1))*COAL*(1-omegaa-omegaw)-...
COAL*(water(i)-water(i-1));
% Estimate unknown species in volatiles (kg/s units)
% Step 1, work out what elemental mass release should be (C, H)
% Assumptions:
% C mass release is proportional to total mass release
% (Asay, 1982)
% H mass release is according to a curve fit to the data
% of Asay (1982).
% O mass release is entirely in the CPD predictions of CO
% H2O and CO2 (Niksa, 1996).
for i = 2:number_reactors+1
Cgoal(i) = COAL*(1-omegaa-omegaw)*yelem(1)*fvol(i)...
-Cgoal(i-1);
Hgoal(i) = COAL*(1-omegaa-omegaw)*yelem(2)*...
((-0.5597)*fvol(i)^2 + 1.5651*fvol(i))...
-Hgoal(i-1);
end
% Step 2, work out what elemental mass release is predicted by
% cpdcp_nlg in H2O, CO2, CH4, CO, and HCN (assume all N release
% is HCN)
for i = 2:number_reactors+1
C_nlg(i) = ffgas(2,i)*fvol(i)*COAL*(1-omegaa-omegaw)*...
(12.011/(12.011+2*15.999))+... % CO2
ffgas(3,i)*fvol(i)*COAL*(1-omegaa-omegaw)*...
(12.011/(12.011+4*1.0079))+... % CH4
ffgas(4,i)*fvol(i)*COAL*(1-omegaa-omegaw)*...
(12.011/(12.011+15.999))+... % CO
fntot(i)*yelem(3)*COAL*(1-omegaa-omegaw)*...
(12.011/14.007)-... % HCN
C_nlg(i-1);
H_nlg(i) = ffgas(1,i)*fvol(i)*COAL*(1-omegaa-omegaw)*...
(2*1.0079/(2*1.0079+15.999))+... % H2O
ffgas(3,i)*fvol(i)*COAL*(1-omegaa-omegaw)*...
(4*1.0079/(12.011+4*1.0079))+... % CH4
fntot(i)*yelem(3)*COAL*(1-omegaa-omegaw)*...
(1.0079/14.007)-... % HCN
H_nlg(i-1);
end
% Step 3, use the difference to get amounts of C and H to
% add (that presumably come from other light gases and cracked
% tars)
Cdiff = Cgoal - C_nlg;
Hdiff = Hgoal - H_nlg;
% Step 4, assume that the other light gas and cracked tar are
% made up of CH4 and C2H2. Calculate the C/H molar ratio and
% use this to choose the proportions of these two gases to make
% that ratio. Make a new array of all estimated volatiles.
CHratio = (Cdiff/12.011)./(Hdiff/1.0079);
ffgas2 = zeros(6,number_reactors+1);
for i = 2:number_reactors+1
MolarProportionCH4 = (CHratio(i)-1)/(0.25-1);
MWunknowns = MolarProportionCH4*(12.011+4*1.0079)+...
(1-MolarProportionCH4)*2*(12.011+1.0079);
MassProportionCH4 = MolarProportionCH4*(12.011+4*1.0079)/...
MWunknowns; % in other light gas and tar
% The values in ffgas are modifed and stored in ffgas2
% which has these species in its rows:
% H2O, CO2, CH4, CO, C2H2, and HCN in kg of species per kg
% of volatiles. After this there are no volatiles of unknown
% or unestimated composition.
ffgas2(1,i) = ffgas(1,i); % H2O
% (excludes evaporated moisture)
ffgas2(2,i) = ffgas(2,i); % CO2
ffgas2(4,i) = ffgas(4,i); % CO
if fvol(i) == 0
ffgas2(6,i) = 0; %HCN
ffgas2(3,i) = 0; %CH4
else
ffgas2(6,i) = (fntot(i)*COAL*(1-omegaa-omegaw)*...
yelem(3)*(1.0079+12.011+14.007)/(14.007))/...
(COAL*(1-omegaa-omegaw)*fvol(i));% HCN
ffgas2(3,i) = ffgas(3,i) + ...
(MassProportionCH4*(ffgas(5,i))*COAL*...
(1-omegaa-omegaw))/...
(fvol(i)*COAL*(1-omegaa-omegaw)); % CH4
end
ffgas2(5,i) = 1-(ffgas2(6,i)+sum(ffgas2(1:4,i)));
% C2H2 (by difference)
end
% ffgas & ffgas2 = fraction of total volatiles that is a
% particular species
% move products downstream, mixing in the volatiles calculated
% above
for i = 2:number_reactors
% get the existing gas composition from previous
% reactor into kg/s units:
ExistingGas = m_dot(i-1)*massFractions(gas(i-1));
% get the gases to be added into kg/s units:
AddedGas = zeros(size(ExistingGas));
AddedGas(H2OIndex) = ffgas2(1,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw)+...
(water(i-1)-water(i))*COAL;
% H2O (light gas) and evaporated moisture
AddedGas(CO2Index) = ffgas2(2,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % CO2
AddedGas(CH4Index) = ffgas2(3,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % CH4
AddedGas(COIndex) = ffgas2(4,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % CO
AddedGas(C2H2Index) = ffgas2(5,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % C2H2
AddedGas(HCNIndex) = ffgas2(6,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % HCN
% add the two compositions from above and insert the gas
% into the next reactor (Heat transfer for volatiles to be
% added to code later)
mixture = ExistingGas + AddedGas;
% lower enthalpy of gas to account for heating of particles
% by convection
H = enthalpy_mass(gas(i-1)); %J/kg
% convheat is in J/particle per second
H = H - (0.5*(convheat(i)+convheat(i-1)))*...
particles*(1/1000)*(tms(i)-tms(i-1))/(m_dot(i-1));
set(dummygas,'H',H,'P',P,'MassFractions',mixture);
insert(upstream(i),dummygas);
end
case 'No',
iterate = false;
close(gcf);
case 'Cancel',
close(gcf);
return; % stop program
end % switch
end % while
% get some variables ready for output
pressurevector(1) = pressure(upstream(1));
for i = 2:number_reactors+1
pressurevector(i) = pressure(gas(i-1));
end
%---END OF IGNITION PART OF CODE-------------------------------------------
Length = Length2;
% Calculate chemical equivalence ratios for ignition section and
% predictions for direct comparison to measurements (NOx, CO, O2, CO2)
MWmix(1) = 1/(sum(massFractions(upstream(1))./molecularWeights(gas(1))'));
nj = massFractions(upstream(1))./molecularWeights(gas(1))';
for k = 1:4
bi(k) = sum(aij(:,k).*nj');
end
V_p = sum(V_plus.*bi);
V_m = sum(V_minus.*bi);
r(1) = -V_p/V_m;
Yi = massFractions(upstream(1));
YWater = Yi(H2OIndex);
XWater = YWater*MWmix(1)/(2*1.0079+15.999);
YNO = Yi(NOIndex);
YNO2 = Yi(NO2Index);
YCO = Yi(COIndex);
YO2 = Yi(O2Index);
YCO2 = Yi(CO2Index);
XNO = YNO*MWmix(1)/(14.007+15.999);
XNO2 = YNO2*MWmix(1)/(14.007+2*15.999);
XCO = YCO*MWmix(1)/(12.011+15.999);
XO2 = YO2*MWmix(1)/(2*15.999);
XCO2 = YCO2*MWmix(1)/(12.011+2*15.999);
NOx_ppm_dry(1) = 1000000*(XNO + XNO2)/(1-XWater);
CO_ppm_dry(1) = 1000000*XCO/(1-XWater);
O2_vol_dry(1) = 100*XO2/(1-XWater);
CO2_vol_dry(1) = 100*XCO2/(1-XWater);
NCE(1) = (YNO*(14.007/(14.007+15.999))+YNO2*(14.007/(14.007+2*15.999)))*...
m_dot(1)/(COAL*(1-omegaa-omegaw)*yelem(3));
for j = 2:number_reactors+1
MWmix(j) = meanMolarMass(gas(j-1));
nj = (1/MWmix(j))*moleFractions(gas(j-1));
for k = 1:4
bi(k) = sum(aij(:,k).*nj);
end
V_p = sum(V_plus.*bi);
V_m = sum(V_minus.*bi);
r(j) = -V_p/V_m;
Xi = moleFractions(gas(j-1));
XWater = Xi(H2OIndex);
XNO = Xi(NOIndex);
XNO2 = Xi(NO2Index);
XCO = Xi(COIndex);
XO2 = Xi(O2Index);
XCO2 = Xi(CO2Index);
NOx_ppm_dry(j) = 1000000*(XNO + XNO2)/(1-XWater);
CO_ppm_dry(j) = 1000000*XCO/(1-XWater);
O2_vol_dry(j) = 100*XO2/(1-XWater);
CO2_vol_dry(j) = 100*XCO2/(1-XWater);
Yi = massFractions(gas(j-1));
YNO = Yi(NOIndex);
YNO2 = Yi(NO2Index);
NCE(j) = (YNO*(14.007/(14.007+15.999))+YNO2*...
(14.007/(14.007+2*15.999)))*m_dot(j-1)/...
(COAL*(1-omegaa-omegaw)*yelem(3));
end
%---START WRITING RESULTS TO FILE------------------------------------------
% File header
fid = fopen(output, 'wt');
fprintf(fid,'BYU MFR Coal Combustion Model\r');
fprintf(fid,'*****************************\r');
fprintf(fid,'* Coded by: Andrew Mackrory *\r');
timedata = fix(clock);
fprintf(fid,['Date: ',date,' \r']);
fprintf(fid,['Time: ',num2str(timedata(4)),':',...
num2str(timedata(5)),' \r\r']);
fprintf(fid,'Model Inputs\r');
fprintf(fid,'============\r');
fprintf(fid,['Methane Through Burner: ',num2str(NG_in),' kg/hr\r']);
fprintf(fid,['Coal: ',num2str(COAL*3600),' kg/hr\r']);
if COAL_Type == 1
fprintf(fid,'Coal Name: Wyoming Subbituminous\r');
elseif COAL_Type == 2
fprintf(fid,'Coal Name: Illinois #6\r');
elseif COAL_Type == 3
fprintf(fid,'Coal Name: Pittsburgh #8\r');
end
fprintf(fid,['Burner Air: ',num2str(Air_in),' kg/hr\r']);
fprintf(fid,['Burner Bottled O2: ',num2str(O2_in),' kg/hr\r']);
fprintf(fid,['Burner Bottled CO2: ',num2str(CO2_in),' kg/hr\r']);
fprintf(fid,['Burnout Air: ',num2str((Air_in/Primary)*(1-Primary)),...
' kg/hr\r']);
fprintf(fid,['Burnout Bottled O2: ',...
num2str((O2_in/Primary)*(1-Primary)),' kg/hr\r']);
fprintf(fid,['Burnout Bottled CO2: ',...
num2str((CO2_in/Primary)*(1-Primary)),' kg/hr\r']);
fprintf(fid,['NO in CO2: ',num2str(NO_doping),' ppm\r']);
fprintf(fid,['Initial Gas Temperature: ',num2str(T1),' K\r']);
fprintf(fid,['Burnout Oxidizer Temperature: ',num2str(T2),' K\r']);
fprintf(fid,['Pressure: ',num2str(press),' atm\r\r']);
if mechanism == 1
fprintf(fid,'Gas Phase Mechanism: GRI 3.0\r');
elseif mechanism == 2
fprintf(fid,'Gas Phase Mechanism: GRI 3.0 + B96\r');
elseif mechanism == 3
fprintf(fid,'Gas Phase Mechanism: SKG03\r');
end
if thermal == 0
fprintf(fid,'Thermal NOx Mechanism Disabled\r');
elseif thermal == 1
fprintf(fid,'Thermal NOx Mechanism Enabled\r');
end
if prompt == 0
fprintf(fid,'Prompt NOx Mechanism Disabled\r\r');
elseif prompt == 1
fprintf(fid,'Prompt NOx Mechanism Enabled\r\r');
end
if gasification == 0
fprintf(fid,'Char gasification by CO2 Disabled\r\r');
elseif gasification == 1
fprintf(fid,'Char gasification by CO2 Enabled\r\r');
end
fprintf(fid,['Percent of O2 Char Oxidation Energy to Char: ',...
num2str(Q_reactO2_x),'\r']);
fprintf(fid,['Percent of CO2 Char Gasification Energy from Char: ',...
num2str(Q_reactCO2_x),'\r\r']);
fprintf(fid,['Burner Emissivity: ',num2str(emiss(1)),'\r']);
fprintf(fid,['Reactor Wall Emissivity: ',num2str(emiss(2)),'\r']);
fprintf(fid,['Exhaust Tube Emissivity: ',num2str(emiss(3)),'\r']);
fprintf(fid,['Coal Particle Emissivity: ',num2str(emiss(4)),'\r']);
fprintf(fid,['Burner Temperature: ',num2str(tbnr),' K\r']);
fprintf(fid,['Exhaust Tube Temperature: ',num2str(texit),' K\r']);
fprintf(fid,'Reactor Wall Temperature Profile (m, K):\r');
for i = 1:length(WallX)
fprintf(fid,[num2str(WallX(i)),'\t',num2str(twallvector(i)),'\r']);
end
fprintf(fid,['\rk_wall ',...
'(Empirical Wall Thermal Conductivity Parameter): ',...
num2str(k_wall),'\r']);
fprintf(fid,['\rCoal Ultimate Analysis: (DAF wt%%)\r',...
'C: ',num2str(yelem(1)*100),...
'\rH: ',num2str(yelem(2)*100),...
'\rN: ',num2str(yelem(3)*100),...
'\rO: ',num2str(yelem(4)*100),...
'\rS: ',num2str(yelem(5)*100),...
' (Not used in model)\r\r']);
fprintf(fid,['Coal Moisture Content (as received): ',...
num2str(omegaw*100),'%%\r']);
fprintf(fid,['Coal Ash Content (as received): ',...
num2str(omegaa*100),'%%\r']);
if ASTMvol ~= 0
fprintf(fid,['ASTM Volatiles (DAF): ',num2str(ASTMvol),' %%\r\r']);
fprintf(fid,'Estimated C13 NMR Parameters for Coal:\r');
fprintf(fid,['mw1: ',num2str(mw1),'\r']);
fprintf(fid,['p0: ',num2str(p0),'\r']);
fprintf(fid,['c0: ',num2str(c0),'\r']);
fprintf(fid,['sigp1: ',num2str(sigp1),'\r']);
fprintf(fid,['mdel: ',num2str(mdel),'\r\r']);
else
fprintf(fid,'\r\r\rC13 NMR Parameters for Coal:\r');
fprintf(fid,['mw1: ',num2str(mw1),'\r']);
fprintf(fid,['p0: ',num2str(p0),'\r']);
fprintf(fid,['c0: ',num2str(c0),'\r']);
fprintf(fid,['sigp1: ',num2str(sigp1),'\r']);
fprintf(fid,['mdel: ',num2str(mdel),'\r\r']);
end
fprintf(fid,['Initial Particle Apparent Density: ',num2str(rhop),...
' g/cm^3\r']);
fprintf(fid,['Initial Particle Diameter: ',num2str(dp*1e4),' um\r']);
fprintf(fid,['Swelling Factor: ',num2str(swell),'\r\r']);
fprintf(fid,'Gas species are reported as mass fractions below.\r\r');
% write molecular weights of species to output file for easy conversion
% from mass fractions to mole fractions with a spreadsheet
fprintf(fid,'\t\tMolecular Weight of Species (kg/kmol):');
MW = molarMasses(gas(1));
for i = 1:nSpecies(gas(1))
fprintf(fid,strcat('\t%d'),MW(i));
end
fprintf(fid,'\r');
fprintf(fid,['Axial Position (m)\tGas Temperature (K)\t',...
'Gas Pressure (Pa)\t']);
for i = 1:nSpecies(gas(1))
text = speciesName(gas(1),i);
fprintf(fid,'%s\t',text{1,1});
end
fprintf(fid,['Gas Mass Flow (kg/s)\tResidence Time (ms)',...
'\tParticle Temperature (K)\tfvol\tfchar\tfcross\tftar\t',...
'fmet\t Trate\tMW\tNsite\tNchar\tfnchar\tfntar\tfnhcn\t',...
'fntot\tfgas\tfH2O\tfCO2\tfCH4\tfCO\tfOth',...
'\tyH2O\tyCO2\tyCH4\tyCO\tyOther\tXgas\tCoal Moisture %%',...
'\tO2 Consumptionby Char (kg/s)',...
'\tCO2 Consumption by Char (kg/s)\tDAF Char Flux (kg/s)',...
'\tGas Phase Chemical Equivalence Ratio',...
'\tNOx (ppm, excluding H2O)\tCO (ppm, excluding H2O)',...
'\tO2 (vol%%, excluding H2O)\tCO2 (vol%%, excluding H2O)',...
'\tNitrogen Conversion Efficiency\tMW of Gas mixture',...
'\tParticle Diameter (um)\r']);
% write model predictions from ignition section of MFR
i = 1;
fprintf(fid,strcat('%d\t %d\t %d\t'),position(i),tgas(i),...
pressurevector(i));
for j = 1:nSpecies(gas(1))
fprintf(fid,'%d\t',massFraction(upstream(1),char(speciesName(gas(1),j))));
end
fprintf(fid,strcat('%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t\t\t\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t %d\r'),...
m_dot(i),tms(i),tp(i),fvol(i),fchar(i),fcross(i),...
ftar(i),fmet(i),trate(i),mwchar(i),yNsite(i),...
fnt(i),fnchar(i),fntar(i),fnhcn(i),fntot(i),...
fgas(i),ffgas(1,i),ffgas(2,i),ffgas(3,i),...
ffgas(4,i),ffgas(5,i),...
yygas(1,i),yygas(2,i),yygas(3,i),yygas(4,i),...
yygas(5,i),yf(i),water(i)*100,r(i),NOx_ppm_dry(i),...
CO_ppm_dry(i),O2_vol_dry(i),CO2_vol_dry(i),NCE(i),...
MWmix(i),dpout(i)*10000);
for i = 2:number_reactors+1
fprintf(fid,strcat('%d\t %d\t %d\t'),position(i),tgas(i),...
pressurevector(i));
for j = 1:nSpecies(gas(1))
fprintf(fid,'%d\t',massFraction(cstr(i-1),char(speciesName(gas(1),j))));
end
fprintf(fid,strcat('%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t\t\t\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\r'),...
m_dot(i-1),tms(i),tp(i),fvol(i),fchar(i),...
fcross(i),ftar(i),fmet(i),trate(i),mwchar(i),...
yNsite(i),fnt(i),fnchar(i),fntar(i),fnhcn(i),...
fntot(i),fgas(i),ffgas(1,i),ffgas(2,i),...
ffgas(3,i),ffgas(4,i),ffgas(5,i),yygas(1,i),...
yygas(2,i),yygas(3,i),yygas(4,i),yygas(5,i),...
yf(i),water(i)*100,r(i),NOx_ppm_dry(i),...
CO_ppm_dry(i),O2_vol_dry(i),CO2_vol_dry(i),...
NCE(i),MWmix(i),dpout(i)*10000);
end
%---BEGINNING OF POST-IGNITION MODEL---------------------------------------
% set up CSTR and initial conditions
i = number_reactors + 1; % i is the variable that keeps track of the
% CSTR number from now on.
if mechanism == 1
maingas = GRI30('Mix');
elseif mechanism == 2
maingas = GRI30_B96('Mix');
elseif mechanism == 3
maingas = SKG03;
end
if mechanism == 3
% set multiplier for Thermal NOx Mechanism
setMultiplier(maingas,74,thermal);
setMultiplier(maingas,73,thermal);
setMultiplier(maingas,72,thermal);
% set multiplier for Prompt NOx Fenimore Mechanism
setMultiplier(maingas,507,prompt);
setMultiplier(maingas,503,prompt);
else
% set multiplier for Thermal NOx Mechanism
setMultiplier(maingas,178,thermal);
setMultiplier(maingas,179,thermal);
setMultiplier(maingas,180,thermal);
% set multiplier for Prompt NOx Fenimore Mechanism
setMultiplier(maingas,239,prompt);
setMultiplier(maingas,240,prompt);
end
ExistingGas = m_dot(i-1)*massFractions(gas(i-1));
% get the gases to be added into kg/s units:
AddedGas = zeros(size(ExistingGas));
AddedGas(H2OIndex) = ffgas2(1,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw)+...
(water(i-1)-water(i))*COAL;
% H2O (light gas) and evaporated moisture
AddedGas(CO2Index) = ffgas2(2,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % CO2
AddedGas(CH4Index) = ffgas2(3,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % CH4
AddedGas(COIndex) = ffgas2(4,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % CO
AddedGas(C2H2Index) = ffgas2(5,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % C2H2
AddedGas(HCNIndex) = ffgas2(6,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % HCN
m_dot(i) = m_dot(i-1) + ...
(fvol(i)-fvol(i-1))*COAL*(1-omegaa-omegaw) - ...
COAL*(water(i)-water(i-1));
% add the two compositions from above and insert the gas
% into the next reactor (Heat transfer for volatiles to be
% added to code later)
mixture = ExistingGas + AddedGas;
% lower enthalpy of gas to account for heating of particles
% by convection
H = enthalpy_mass(gas(i-1)); %J/kg
% convheat is in J/particle per second
H = H - (0.5*(convheat(i)+convheat(i-1)))*...
particles*(1/1000)*(tms(i)-tms(i-1))/(m_dot(i-1));
set(maingas,'H',H,'P',P,'MassFractions',mixture);
% clear old objects from ignition section:
clear upstream gas mfc cstr v downstream exhaust gw w network
% create new reactor network objects for this section of code:
% follows same procedure as above code
upstream = Reservoir(maingas);
insert(upstream,maingas);
cstr = Reactor(maingas);
setInitialVolume(cstr,Volume(end));
w = Wall;
install(w,cstr,env);
setArea(w,pi*d*Length);
R = log(0.180/(d/2))/(2*pi*Length*k_wall);
setThermalResistance(w,R);
mfc = MassFlowController;
install(mfc,upstream,cstr);
setMassFlowRate(mfc,m_dot(end));
exhaust = air;
set(exhaust,'T',300,'P',P);
downstream = Reservoir(exhaust);
v = Valve;
install(v,cstr,downstream);
setValveCoeff(v,1.0);
network_cell_array = {cstr};
network = ReactorNet(network_cell_array);
% start loop for reactors
tme = 0;
devol_incomplete = true;
while devol_incomplete
i = i + 1 %#ok<NOPTS>
position(i) = position(i-1) + Length;
% integrate the CSTR long enough to reach steady state
count = 1;
old_T = temperature(cstr);
delta_T = 1;
while (delta_T > tolerance) || (count < 3)
tme = tme + dt;
flagtime = 0;
errorcount = 0;
while flagtime < tme
try 
advance(network, tme);
flagtime = tme;
catch exception
errorcount = errorcount + 1;
if errorcount > 10
rethrow(exception);
end
end
end
new_T = temperature(cstr);
delta_T = abs(new_T-old_T);
old_T = new_T;
count = count + 1;
end
% update last point in output data variables
tfilm = 0.5*(tp + tgas);
tgas(i) = temperature(maingas);
pressurevector(i) = pressure(maingas);
% calculate/estimate properties required for CPD model
velocity(i) = m_dot(i-1)/(density(upstream)*Area);
rhogas(i) = density(upstream);
xwbvector(i) = moleFraction(maingas,'H2O');
set(dummygas,'T',tfilm(i-1),'P',P,'Y',massFractions(upstream));
ugvector(i) = viscosity(dummygas);
kgvector(i) = thermalConductivity(dummygas);
DiffCoeffs = mixDiffCoeffs(dummygas);
diffwvector(i) = DiffCoeffs(H2OIndex);
% Mixture-averaged diffusion coefficient (m^2/s) for water
clear DiffCoeffs;
cpgvector(i) = cp_mass(dummygas);
prgas(i) = ugvector(i)*cpgvector(i)/kgvector(i);
% Call the CPD model
[tms,xm,tp,tg,fvol,fchar,fcross,ftar,fmet,trate,mwchar,yNsite,...
fnt,fnchar,fntar,fnhcn,fntot,fgas,ffgas,yygas,yf,water,convheat,...
dpout] = cpdcp_nlg(twallvector,tbnr,texit,timax,yelem,mw1,p0,c0,...
sigp1,mdel,position*100,tgas,position*100,velocity*100,press,...
tgas(1),velocity(1)*100,rhop,dp,swell,omegaw,omegaa,...
rhogas/1000,xwbvector,ugvector*10,kgvector/418.4,...
diffwvector*0.01^2,cpgvector*2.38846e-4,prgas,emiss,d*50,WallX*100);
% check if devolatilization is complete
if ((fvol(i)- fvol(i-1)) <= 1e-4) && (fvol(i) > 0.4)
devol_incomplete = false;
end
% Add the products of devolatilization to the CSTR products
m_dot(i) = m_dot(i-1) + ...
(fvol(i)-fvol(i-1))*COAL*(1-omegaa-omegaw) - ...
COAL*(water(i)-water(i-1));
setMassFlowRate(mfc,m_dot(i));
Cgoal(i) = COAL*(1-omegaa-omegaw)*yelem(1)*fvol(i)...
-Cgoal(i-1);
Hgoal(i) = COAL*(1-omegaa-omegaw)*yelem(2)*...
((-0.5597)*fvol(i)^2 + 1.5651*fvol(i))...
-Hgoal(i-1);
C_nlg(i) = ffgas(2,i)*fvol(i)*COAL*(1-omegaa-omegaw)*...
(12.011/(12.011+2*15.999))+... % CO2
ffgas(3,i)*fvol(i)*COAL*(1-omegaa-omegaw)*...
(12.011/(12.011+4*1.0079))+... % CH4
ffgas(4,i)*fvol(i)*COAL*(1-omegaa-omegaw)*...
(12.011/(12.011+15.999))+... % CO
fntot(i)*yelem(3)*COAL*(1-omegaa-omegaw)*...
(12.011/14.007)-... % HCN
C_nlg(i-1);
H_nlg(i) = ffgas(1,i)*fvol(i)*COAL*(1-omegaa-omegaw)*...
(2*1.0079/(2*1.0079+15.999))+... % H2O
ffgas(3,i)*fvol(i)*COAL*(1-omegaa-omegaw)*...
(4*1.0079/(12.011+4*1.0079))+... % CH4
fntot(i)*yelem(3)*COAL*(1-omegaa-omegaw)*...
(1.0079/14.007)-... % HCN
H_nlg(i-1);
Cdiff = Cgoal - C_nlg;
Hdiff = Hgoal - H_nlg;
CHratio = (Cdiff/12.011)./(Hdiff/1.0079);
MolarProportionCH4 = (CHratio(i)-1)/(0.25-1);
MWunknowns = MolarProportionCH4*(12.011+4*1.0079)+...
(1-MolarProportionCH4)*2*(12.011+1.0079);
MassProportionCH4 = MolarProportionCH4*(12.011+4*1.0079)/...
MWunknowns; % in other light gas and tar
% The ffgas array is now modifed and becomes ffgas2
% which has these species in its rows:
% H2O, CO2, CH4, CO, C2H2, and HCN. After this there are no
% volatiles of unestimated composition.
ffgas2(1,i) = ffgas(1,i); % H2O
% (excludes evaporated moisture)
ffgas2(2,i) = ffgas(2,i); % CO2
ffgas2(4,i) = ffgas(4,i); % CO
ffgas2(6,i) = (fntot(i)*COAL*(1-omegaa-omegaw)*yelem(3)*...
(1.0079+12.011+14.007)/(14.007))/...
(COAL*(1-omegaa-omegaw)*fvol(i)); % HCN
ffgas2(3,i) = ffgas(3,i) + ...
(MassProportionCH4*(ffgas(5,i))*COAL*(1-omegaa-omegaw))/...
(fvol(i)*COAL*(1-omegaa-omegaw)); % CH4
ffgas2(5,i) = 1-(ffgas2(6,i)+sum(ffgas2(1:4,i)));
% C2H2 (by difference)
ExistingGas = m_dot(i-1)*massFractions(maingas);
% get the gases to be added into kg/s units:
AddedGas = zeros(size(ExistingGas));
AddedGas(H2OIndex) = ffgas2(1,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw)+...
(water(i-1)-water(i))*COAL;
% H2O (light gas) and evaporated moisture
AddedGas(CO2Index) = ffgas2(2,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % CO2
AddedGas(CH4Index) = ffgas2(3,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % CH4
AddedGas(COIndex) = ffgas2(4,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % CO
AddedGas(C2H2Index) = ffgas2(5,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % C2H2
AddedGas(HCNIndex) = ffgas2(6,i)*(fvol(i)-fvol(i-1))*...
COAL*(1-omegaa-omegaw); % HCN
% add the two compositions from above and insert the gas
% into the next reactor (Heat transfer for volatiles to be
% added to code later - ie heat brought into gas by volatiles
mixture = ExistingGas + AddedGas;
% lower enthalpy of gas to account for heating of particles
% by convection
H = enthalpy_mass(maingas); %J/kg
% convheat is in J/particle per second
H = H - (0.5*(convheat(i)+convheat(i-1)))*...
particles*(1/1000)*(tms(i)-tms(i-1))/(m_dot(i-1));
set(maingas,'H',H,'P',P,'MassFractions',mixture);
insert(upstream,maingas);
% Calculate chemical equivalence ratio for reactor i
MWmix(i) = meanMolarMass(maingas);
nj = (1/MWmix(i))*moleFractions(maingas);
for k = 1:4
bi(k) = sum(aij(:,k).*nj);
end
V_p = sum(V_plus.*bi);
V_m = sum(V_minus.*bi);
r(i) = -V_p/V_m;
Xi = moleFractions(maingas);
XWater = Xi(H2OIndex);
XNO = Xi(NOIndex);
XNO2 = Xi(NO2Index);
XCO = Xi(COIndex);
XO2 = Xi(O2Index);
XCO2 = Xi(CO2Index);
NOx_ppm_dry(i) = 1000000*(XNO + XNO2)/(1-XWater);
CO_ppm_dry(i) = 1000000*XCO/(1-XWater);
O2_vol_dry(i) = 100*XO2/(1-XWater);
CO2_vol_dry(i) = 100*XCO2/(1-XWater);
Yi = massFractions(maingas);
YNO = Yi(NOIndex);
YNO2 = Yi(NO2Index);
NCE(i) = (YNO*(14.007/(14.007+15.999))+...
YNO2*(14.007/(14.007+2*15.999)))*...
m_dot(i-1)/(COAL*(1-omegaa-omegaw)*yelem(3));
% Output data point to file
fprintf(fid,strcat('%d\t %d\t %d\t'),position(i),tgas(i),...
pressurevector(i));
for j = 1:nSpecies(maingas)
fprintf(fid,'%d\t',massFraction(cstr,char(speciesName(maingas,j))));
end
fprintf(fid,strcat('%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t\t\t\t %d\t %d\t %d\t %d\t %d\t %d\t',...
'%d\t %d\r'),...
m_dot(i-1),tms(i),tp(i),fvol(i),fchar(i),...
fcross(i),ftar(i),fmet(i),trate(i),mwchar(i),...
yNsite(i),fnt(i),fnchar(i),fntar(i),fnhcn(i),...
fntot(i),fgas(i),ffgas(1,i),ffgas(2,i),...
ffgas(3,i),ffgas(4,i),ffgas(5,i),yygas(1,i),...
yygas(2,i),yygas(3,i),yygas(4,i),yygas(5,i),...
yf(i),water(i)*100,r(i),NOx_ppm_dry(i),...
CO_ppm_dry(i),O2_vol_dry(i),CO2_vol_dry(i),...
NCE(i),MWmix(i),dpout(i)*10000);
end
% Update particle diameter to account for swelling during devolatilization
dp = dpout(end);
m_dot(i) = m_dot(i-1);
% Calculate mass flux of DAF CHAR (kg/s) and ASHratio (kg_ash/kg_DAF_CHAR)
% in the CHAR
% (Char is assumed C(s) for remainder of code)
% (Ash percentage in char is assumed to remain constant - shrinking core
% model - ash assumed shed from surface of char as char reacts)
CHAR = COAL*(1-omegaa-omegaw)*fchar(i);
ASHratio = COAL*omegaa/(CHAR);
% Calculate density of char (including ash) - this remains constant in the
% shrinking core model employed by Goetz et al. in their data reduction
% (according to Smith and Smoot)
% kg/m^3
rhoCHAR = (COAL*(1-omegaw)*fchar(i))/(particles*((4/3)*pi*(dp/200)^3));
% (dp in cm)
% Post-devolatilization model (prior to burnout oxidizer addition)
% Models char (C(s)) oxidation and gasification by O2 and CO2 respectively
% Set kinetic constants for the char (Values from Goetz et al. (1982) as
% presented in "Coal Combustion and Gasification" Chapter 4 by Smoot and
% Smith)
gas_const = 1.987; % cal/gmol.K
sigma = 1.335e-12; % cal/s cm^2 K^4 radiation constant
% Get vector of wall temps (interpolated from measurements in input file)
max_position = WallX(end);
max_index = floor(max_position/Length);
j = 1;
Current_position = 0;
twall = zeros(1,max_index);
for k = 1:max_index
twall(k) = ((Current_position-WallX(j))/(WallX(j+1)-WallX(j)))*...
(twallvector(j+1)-twallvector(j))+twallvector(j);
Current_position = Current_position + Length;
% Update interpolation index if necessary
if (Current_position > WallX(j+1))
j = j + 1;
end
end
if COAL_Type == 1
AO2 = 145; % g/(cm2s atmO2)
%Units are g/s of O2 consumed per unit char external surface
%area per atmosphere of O2 partial pressure
EAO2 = 19970; % cal/gmole (units to match R units)
ACO2 = 1040*gasification; % g/(cm2s atmCO2)
EACO2 = 42470; % cal/gmole
elseif COAL_Type == 2
AO2 = 60; % g/(cm2s atmO2)
%Units are g/s of O2 consumed per unit char external surface
%area per atmosphere of O2 partial pressure
EAO2 = 17150; % cal/gmole (units to match R units)
ACO2 = 12973*gasification; % g/(cm2s atmCO2)
EACO2 = 56368; % cal/gmole
elseif COAL_Type == 3
AO2 = 66; % g/(cm2s atmO2)
%Units are g/s of O2 consumed per unit char external surface
%area per atmosphere of O2 partial pressure
EAO2 = 20360; % cal/gmole (units to match R units)
ACO2 = 1390*gasification; % g/(cm2s atmCO2)
EACO2 = 53700; % cal/gmole
end
% C(s) + O2 -> CO
% C(s) + CO2 -> 2CO (This reaction may be disabled by setting
% gasification = 0)
if mechanism == 1
chargas = GRI30('Mix');
elseif mechanism == 2
chargas = GRI30_B96('Mix');
elseif mechanism == 3
chargas = SKG03;
end
while (position(end) < 0.667)
i = i + 1 %#ok<NOPTS>
position(i) = position(i-1) + Length;
% Calculate surface area of particles in previous reactor section
if CHAR < 0
A_char = 0;
else
A_char = particles*(tms(i-1)-tms(i-2))*(0.001)*4*pi*(dp/2)^2; %cm^2
end
% Calculate rate of oxidizer consumption in previous reactor section
% kp units: g of oxidizer consumed per second per unit char surface
% area (in cm^2) per atm of oxidizer
kpO2 = AO2*exp(-EAO2/(gas_const*tp(i-1)));
kpCO2 = ACO2*exp(-EACO2/(gas_const*tp(i-1)));
% rpo units: kg/s oxidizer consumed (in previous reactor section)
% Minimum of:
% 1. oxidizer consumption predicted by kinetic rate
% 2. oxidizer available
ExistingGas = m_dot(i-1)*massFractions(maingas); % kg/s
AddedGas = zeros(size(ExistingGas));
rpoO2 = min([0.001*kpO2*A_char*moleFraction(maingas,'O2')*P/101325;
ExistingGas(O2Index)]); % kg/s
rpoCO2 = min([0.001*kpCO2*A_char*moleFraction(maingas,'CO2')*P/101325;
ExistingGas(CO2Index)]); % kg/s
% Particle heating from convection and radiation (and some fraction of
% heat of char reaction)
% Estimated heat capacity of the char (using same methods as CPD model)
% get daf coal heat capacity
[cpc] = heatcp(tp(i-1),yelem); % cal/g/K
% get ash heat capacity
[cpa] = heatap(tp(i-1)); % cal/g/K
% combine heat capacities
cp = (CHAR*cpc + (ASHratio*CHAR)*cpa)/(CHAR*(1+ASHratio));
% convert to J/kg/K
cp = cp*4186.8;
% Radiation heat exchange with reactor walls (copied from CPD function
% and therefore has cm, cal, g units
z = position(i)*100; % distance from burner (cm)
% distance from exhaust will be (200-z) because reactor is 2m
% long
% set up areas, etc of radiation enclosure
A(1) = pi*(d*50)^2; % burner (Area vector = same order as emiss vector)
A(2) = 2*pi*(d*50)*200; % walls
A(3) = A(1); % exhaust
A(4) = A_char; % particles
Temp(1) = tbnr;
Temp(2) = twall(i-1);
Temp(3) = texit;
Temp(4) = tp(i-1);
F(4,1) = 0.5*(1-(1/(1+((d*50)/z)^2)^0.5)); %disk to sphere view factor
F(4,3) = 0.5*(1-(1/(1+((d*50)/(200-z))^2)^0.5)); %same as above
F(4,2) = 1 - (F(4,1) + F(4,3)); % by summation rule
% by reciprocity:
F(1,4) = F(4,1)*A(4)/A(1);
F(2,4) = F(4,2)*A(4)/A(2);
F(3,4) = F(4,3)*A(4)/A(3);
qrad = emiss(4)*( F(1,4)*A(1)*emiss(1)*sigma*Temp(1)^4+...
F(2,4)*A(2)*emiss(2)*sigma*Temp(2)^4+...
F(3,4)*A(3)*emiss(3)*sigma*Temp(3)^4-...
A(4)*sigma*Temp(4)^4); % cal/s
Q_rad = qrad*4.1868; % J/s
% Convection heat exchange with the gas
nu = 2; % Assumes entrained particles - Reynolds number is zero
rtot = 1000*(rpoO2*12.011/31.998 + rpoCO2*12.011/44.009);
b = cp_mass(maingas)*2.38846e-4*(rtot)/...
(2.0*pi*dp*thermalConductivity(maingas)/418.4);
if (b >= 1.e-4)
blow = b/(exp(b)-1);
else
blow = 1.0;
end
h = blow*nu*(thermalConductivity(maingas)/418.4)/dp;
qconv = h*A_char*(tgas(i-1)-tp(i-1)); % cal/s
Q_conv = qconv*4.1868; % J/s
% Some fraction (Q_react_x) of Heat released by reaction
% 6908557 J/kg_O2 for 2C(s) + O2 -> 2CO
% -3918744 J/kg_CO2 for C(s) + CO2 -> 2CO (negative means endothermic)
% Assume some fraction of heat goes into char, char temp will be
% adjusted and CO is released at the new particle temperature for
% mixing with the gas phase
Q_react = Q_reactO2_x*(6908557*rpoO2) -...
Q_reactCO2_x*(3918744*rpoCO2); % J/s
Q_total = Q_rad + Q_conv + Q_react; % J/s
tp(i) = tp(i-1) + (Q_total/(CHAR+COAL*omegaa))/cp;
% Remove consumed oxidizer from maingas and add products (CO) at
% particle temperature
% Heat released to the CO product gas by reaction
Q_react = (1-Q_reactO2_x)*(6908557*rpoO2) -...
(1-Q_reactCO2_x)*(3918744*rpoCO2); % J/s
AddedGas(COIndex) = 2*rpoO2*28.01/31.998+2*rpoCO2*28.01/44.009; %kgCO/s
set(chargas,'T',tp(i),'P',P,'MassFractions',AddedGas);
H_in = enthalpy_mass(chargas)*sum(AddedGas)+Q_react; % J/s
AddedGas(COIndex) = 0;
AddedGas(O2Index) = rpoO2;
AddedGas(CO2Index) = rpoCO2;
set(chargas,'T',tgas(i-1),'P',P,'MassFractions',AddedGas);
H_out = enthalpy_mass(chargas)*sum(AddedGas); % J/s
H_old = enthalpy_mass(maingas)*sum(ExistingGas); % J/s
AddedGas(COIndex) = 2*rpoO2*28.01/31.998 + 2*rpoCO2*28.01/44.009;
AddedGas(O2Index) = -rpoO2;
AddedGas(CO2Index) = -rpoCO2;
mixture = AddedGas + ExistingGas; % kg/s
H_new = (H_old + H_in - H_out)/sum(mixture); % J/kg
if CHAR < 0
% do nothing
else
set(maingas,'H',H_new,'P',P,'MassFractions',mixture);
end
insert(upstream,maingas);
% Remove C(s) mass from char and update gas mass flow rate and particle
% diameter
CHAR = CHAR - (rpoO2*2*12.011/31.998 + rpoCO2*12.011/44.009); % kgC/s
m_dot(i) = sum(mixture);
setMassFlowRate(mfc,m_dot(i));
dp = 200*(((3/(4*pi))*CHAR*(1+ASHratio)/(particles*rhoCHAR))^(1/3));
% integrate the CSTR long enough to reach steady state
count = 1;
old_T = temperature(cstr);
delta_T = 1;
while ((delta_T > tolerance) || (count < 3)) && (count < 10000)
tme = tme + dt;
flagtime = 0;
errorcount = 0;
while flagtime < tme
try 
advance(network, tme);
flagtime = tme;
catch exception
errorcount = errorcount + 1;
if errorcount > 10
rethrow(exception);
end
end
end
new_T = temperature(cstr);
delta_T = abs(new_T-old_T);
old_T = new_T;
count = count + 1;
end
count %#ok<NOPTS>
% update last point in output data variables
tgas(i) = temperature(maingas);
pressurevector(i) = pressure(maingas);
velocity(i) = m_dot(i-1)/(density(upstream)*Area);
tms(i) = tms(i-1) + 1000*(Length/velocity(i));
% Calculate chemical equivalence ratio for reactor i
MWmix(i) = meanMolarMass(maingas);
nj = (1/MWmix(i))*moleFractions(maingas);
for k = 1:4
bi(k) = sum(aij(:,k).*nj);
end
V_p = sum(V_plus.*bi);
V_m = sum(V_minus.*bi);
r(i) = -V_p/V_m;
Xi = moleFractions(maingas);
XWater = Xi(H2OIndex);
XNO = Xi(NOIndex);
XNO2 = Xi(NO2Index);
XCO = Xi(COIndex);
XO2 = Xi(O2Index);
XCO2 = Xi(CO2Index);
NOx_ppm_dry(i) = 1000000*(XNO + XNO2)/(1-XWater);
CO_ppm_dry(i) = 1000000*XCO/(1-XWater);
O2_vol_dry(i) = 100*XO2/(1-XWater);
CO2_vol_dry(i) = 100*XCO2/(1-XWater);
Yi = massFractions(maingas);
YNO = Yi(NOIndex);
YNO2 = Yi(NO2Index);
NCE(i) = (YNO*(14.007/(14.007+15.999))+YNO2*...
(14.007/(14.007+2*15.999)))*...
m_dot(i-1)/(COAL*(1-omegaa-omegaw)*yelem(3));
% Output data point to file
fprintf(fid,strcat('%d\t %d\t %d\t'),position(i),tgas(i),...
pressurevector(i));
for j = 1:nSpecies(maingas)
fprintf(fid,'%d\t',massFraction(cstr,char(speciesName(maingas,j))));
end
fprintf(fid,strcat('%d\t %d\t %d\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t',...
'\t\t\t\t\t\t\t\t\t\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t %d\t %d\t %d\t %d\t %d\r'),...
m_dot(i-1),tms(i),tp(i),rpoO2,rpoCO2,CHAR,r(i),...
NOx_ppm_dry(i),CO_ppm_dry(i),O2_vol_dry(i),...
CO2_vol_dry(i),NCE(i),MWmix(i),dp*10000);
end
% Final Code Section: Burnout Oxidizer Added for Char Burnout
% (instant or intense mixing is assumed)
ExistingGas = m_dot(i-1)*massFractions(maingas); % kg/s
H_Primary = enthalpy_mass(maingas)*sum(ExistingGas); % J/s
AddedGas = zeros(size(ExistingGas));
AddedGas(O2Index) = (O2_in/Primary)*(1-Primary)+...
0.232999715*(Air_in/Primary)*(1-Primary);
AddedGas(N2Index) = ((1-0.232999715)*(Air_in/Primary)+(N2_in/Primary))...
*(1-Primary);
AddedGas(CO2Index) = (CO2_in/Primary)*(1-Primary);
if AddedGas(CO2Index) > 0
AddedGas(NOIndex) = ((NO_doping/1000000)*(14.007+15.999)/...
(12.011+2*15.999))/AddedGas(CO2Index);
end
AddedGas = AddedGas/3600; % convert to kg/s
set(chargas,'T',T2,'P',P,'MassFractions',AddedGas);
H_Burnout = enthalpy_mass(chargas)*sum(AddedGas); % J/s
mixture = ExistingGas+AddedGas; % kg/s
H_new = (H_Primary + H_Burnout)/sum(mixture); % J/kg
set(maingas,'H',H_new,'P',P,'MassFractions',mixture);
m_dot(i) = sum(mixture);
while (position(end) < (WallX(end)-Length))
i = i + 1 %#ok<NOPTS>
position(i) = position(i-1) + Length;
% Calculate area of particles in previous reactor section
if CHAR < 0
A_char = 0;
else
A_char = particles*(tms(i-1)-tms(i-2))*(0.001)*4*pi*(dp/2)^2; %cm^2
end
% Calculate rate of oxidizer consumption in previous reactor section
% kp units: g of oxidizer consumed per second per unit char surface
% area (in cm^2) per atm of oxidizer
kpO2 = AO2*exp(-EAO2/(gas_const*tp(i-1)));
kpCO2 = ACO2*exp(-EACO2/(gas_const*tp(i-1)));
% rpo units: kg/s oxidizer consumed (in previous reactor section)
% Minimum of:
% 1. oxidizer consumption predicted by kinetic rate
% 2. oxidizer available
ExistingGas = m_dot(i-1)*massFractions(maingas); % kg/s
AddedGas = zeros(size(ExistingGas));
rpoO2 = min([0.001*kpO2*A_char*moleFraction(maingas,'O2')*P/101325;
ExistingGas(O2Index)]); % kg/s
rpoCO2 = min([0.001*kpCO2*A_char*moleFraction(maingas,'CO2')*P/101325;
ExistingGas(CO2Index)]); % kg/s
% Particle heating from convection and radiation (and some fraction of
% heat of char reaction)
% Estimated heat capacity of the char (using same methods as CPD model)
% get daf coal heat capacity
[cpc] = heatcp(tp(i-1),yelem); % cal/g/K
% get ash heat capacity
[cpa] = heatap(tp(i-1)); % cal/g/K
% combine heat capacities
cp = (CHAR*cpc + (ASHratio*CHAR)*cpa)/(CHAR*(1+ASHratio));
% convert to J/kg/K
cp = cp*4186.8;
% Radiation heat exchange with reactor walls (copied from CPD function
% and therefore has cm, cal, g units)
z = position(i)*100; % distance from burner (cm)
% distance from exhaust will be (200-z) because reactor is 2m
% long
% set up areas, etc of radiation enclosure
A(1) = pi*(d*50)^2; % burner (Area vector = same order as emiss vector)
A(2) = 2*pi*(d*50)*200; % walls
A(3) = A(1); % exhaust
A(4) = A_char; % particles
Temp(1) = tbnr;
Temp(2) = twall(i-1);
Temp(3) = texit;
Temp(4) = tp(i-1);
F(4,1) = 0.5*(1-(1/(1+((d*50)/z)^2)^0.5)); %disk to sphere view factor
F(4,3) = 0.5*(1-(1/(1+((d*50)/(200-z))^2)^0.5)); %same as above
F(4,2) = 1 - (F(4,1) + F(4,3)); % by summation rule
% by reciprocity:
F(1,4) = F(4,1)*A(4)/A(1);
F(2,4) = F(4,2)*A(4)/A(2);
F(3,4) = F(4,3)*A(4)/A(3);
qrad = emiss(4)*( F(1,4)*A(1)*emiss(1)*sigma*Temp(1)^4+...
F(2,4)*A(2)*emiss(2)*sigma*Temp(2)^4+...
F(3,4)*A(3)*emiss(3)*sigma*Temp(3)^4-...
A(4)*sigma*Temp(4)^4); % cal/s
Q_rad = qrad*4.1868; % J/s
% Convection heat exchange with the gas
nu = 2; % Assumes entrained particles - Reynolds number is zero
rtot = 1000*(rpoO2*12.011/31.998 + rpoCO2*12.011/44.009);
b = cp_mass(maingas)*2.38846e-4*(rtot)/...
(2.0*pi*dp*thermalConductivity(maingas)/418.4);
if (b >= 1.e-4)
blow = b/(exp(b)-1);
else
blow = 1.0;
end
h = blow*nu*(thermalConductivity(maingas)/418.4)/dp;
qconv = h*A_char*(tgas(i-1)-tp(i-1)); % cal/s
Q_conv = qconv*4.1868; % J/s
% Some fraction (Q_react_x) of Heat released by reaction
% 6908557 J/kg_O2 for 2C(s) + O2 -> 2CO
% -3918744 J/kg_CO2 for C(s) + CO2 -> 2CO (negative means endothermic)
% Assume some fraction of heat goes into char, char temp will be
% adjusted and CO is released at the new particle temperature for
% mixing with the gas phase
Q_react = Q_reactO2_x*(6908557*rpoO2) -...
Q_reactCO2_x*(3918744*rpoCO2); % J/s
Q_total = Q_rad + Q_conv + Q_react; % J/s
tp(i) = tp(i-1) + (Q_total/(CHAR+COAL*omegaa))/cp;
% Remove consumed oxidizer from maingas and add products (CO) at
% particle temperature
% Heat released to the CO product gas by reaction
Q_react = (1-Q_reactO2_x)*(6908557*rpoO2) -...
(1-Q_reactCO2_x)*(3918744*rpoCO2); % J/s
AddedGas(COIndex) = 2*rpoO2*28.01/31.998+2*rpoCO2*28.01/44.009; %kgCO/s
set(chargas,'T',tp(i),'P',P,'MassFractions',AddedGas);
H_in = enthalpy_mass(chargas)*sum(AddedGas)+Q_react; % J/s
AddedGas(COIndex) = 0;
AddedGas(O2Index) = rpoO2;
AddedGas(CO2Index) = rpoCO2;
set(chargas,'T',tgas(i-1),'P',P,'MassFractions',AddedGas);
H_out = enthalpy_mass(chargas)*sum(AddedGas); % J/s
H_old = enthalpy_mass(maingas)*sum(ExistingGas); % J/s
AddedGas(COIndex) = 2*rpoO2*28.01/31.998 + 2*rpoCO2*28.01/44.009;
AddedGas(O2Index) = -rpoO2;
AddedGas(CO2Index) = -rpoCO2;
mixture = AddedGas + ExistingGas; % kg/s
H_new = (H_old + H_in - H_out)/sum(mixture); % J/kg
if CHAR < 0
% do nothing
else
set(maingas,'H',H_new,'P',P,'MassFractions',mixture);
end
insert(upstream,maingas);
% Remove C(s) mass from char and update gas mass flow rate and particle
% diameter
CHAR = CHAR - (rpoO2*2*12.011/31.998 + rpoCO2*12.011/44.009); % kgC/s
m_dot(i) = sum(mixture);
setMassFlowRate(mfc,m_dot(i));
dp = 200*(((3/(4*pi))*CHAR*(1+ASHratio)/(particles*rhoCHAR))^(1/3));
% integrate the CSTR long enough to reach steady state
count = 1;
old_T = temperature(cstr);
delta_T = 1;
while ((delta_T > tolerance) || (count < 3)) && (count < 10000)
tme = tme + dt;
flagtime = 0;
errorcount = 0;
while flagtime < tme
try 
advance(network, tme);
flagtime = tme;
catch exception
errorcount = errorcount + 1;
if errorcount > 10
rethrow(exception);
end
end
end
new_T = temperature(cstr);
delta_T = abs(new_T-old_T);
old_T = new_T;
count = count + 1;
end
count %#ok<NOPTS>
% update last point in output data variables
tgas(i) = temperature(maingas);
pressurevector(i) = pressure(maingas);
velocity(i) = m_dot(i-1)/(density(upstream)*Area);
tms(i) = tms(i-1) + 1000*(Length/velocity(i));
% Calculate chemical equivalence ratio for reactor i
MWmix(i) = meanMolarMass(maingas);
nj = (1/MWmix(i))*moleFractions(maingas);
for k = 1:4
bi(k) = sum(aij(:,k).*nj);
end
V_p = sum(V_plus.*bi);
V_m = sum(V_minus.*bi);
r(i) = -V_p/V_m;
Xi = moleFractions(maingas);
XWater = Xi(H2OIndex);
XNO = Xi(NOIndex);
XNO2 = Xi(NO2Index);
XCO = Xi(COIndex);
XO2 = Xi(O2Index);
XCO2 = Xi(CO2Index);
NOx_ppm_dry(i) = 1000000*(XNO + XNO2)/(1-XWater);
CO_ppm_dry(i) = 1000000*XCO/(1-XWater);
O2_vol_dry(i) = 100*XO2/(1-XWater);
CO2_vol_dry(i) = 100*XCO2/(1-XWater);
Yi = massFractions(maingas);
YNO = Yi(NOIndex);
YNO2 = Yi(NO2Index);
NCE(i) = (YNO*(14.007/(14.007+15.999))+YNO2*...
(14.007/(14.007+2*15.999)))*...
m_dot(i-1)/(COAL*(1-omegaa-omegaw)*yelem(3));
% Output data point to file
fprintf(fid,strcat('%d\t %d\t %d\t'),position(i),tgas(i),...
pressurevector(i));
for j = 1:nSpecies(maingas)
fprintf(fid,'%d\t',massFraction(cstr,char(speciesName(maingas,j))));
end
fprintf(fid,strcat('%d\t %d\t %d\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t',...
'\t\t\t\t\t\t\t\t\t\t %d\t %d\t %d\t %d\t',...
'%d\t %d\t %d\t %d\t %d\t %d\t %d\r'),...
m_dot(i-1),tms(i),tp(i),rpoO2,rpoCO2,CHAR,r(i),...
NOx_ppm_dry(i),CO_ppm_dry(i),O2_vol_dry(i),...
CO2_vol_dry(i),NCE(i),MWmix(i),dp*10000);
end
% close the output file
fclose(fid);
%---LIST OF VARIABLES (alphabetical)---------------------------------------
% Note that MATLAB is case-sensitive.
%
% Units specified here are as at the first instance of the variable in the
% code (such as the user-input values). Unit conversions are noted in the
% code when they are changed.
%
% A - array of areas used in radiation heat transfer calculations (cm^2)
% A_char - surface area of char (cm^2)
% ACO2 - pre-exponential factor for char gasification by CO2
% (g/(cm2s atmCO2)
% AddedGas - vector of mass flow rates of species to be added to the gas
% mixture due to coal devolatilization and heterogeneous
% reactions (kg/s)
% aij - an array of the number of atoms of C, H, N, and O in the species in
% the gas objects - species in rows, number of CHNO atoms
% in columns
% Air_in - flow rate of air through the burner (kg/hr)
% ambient - a gas mixture object used to create the env reservoir
% AO2 - pre-exponential factor for char oxidiation (g/(cm2s atmO2))
% Area - cross sectional area of reactor (m^2)
% ASHratio - (kg_ash/kg_DAF_CHAR) in the CHAR
% ASTMvol - ASTM proximate analysis volatiles, DAF (0 < ASTMvol < 100)
% b - transfer number for effects of high mass transfer
% bi - intermediate variable in calculation of r
% blow - blowing factor = b/(exp(b)-1) for effects of high mass transfer
% ButtonName - yes, no, or cancel result of asking user whether to continue
% iterations in ignition section of code
% c - vector of empirical coefficients used when estimating C13 NMR
% parameters from proximate and ultimate analysis results
% c0 - char bridge population (C13 NMR parameter)
% C2H2Index - Index of the species in the gas mixture Cantera objects
% Cdiff - the difference between the carbon release predicted by CPD and
% the amount of carbon release expected from correlation
% with total mass loss
% Cgoal - the target amount of carbon release from volatiles correlated
% with total mass release
% CH4Index - Index of the species in the gas mixture Cantera objects
% CHAR - mass flux of DAF CHAR (kg/s)
% chargas - Cantera gas mixture object used for char reactions section of
% model
% CHNOIndex - indices for the elements in the Cantera objects
% CHratio - the CH ratio of the volatiles not predicted as specific species
% by the CPD model
% C_nlg - the carbon release predicted by the CPD model
% CO2 - moles/hr bottled CO2
% CO2_in - flow rate of bottled CO2 through the burner (kg/hr)
% CO2_vol_dry - CO2 expressed in units of vol% on a water-free basis
% COIndex - Index of the species in the gas mixture Cantera objects
% CO_ppm_dry - CO expressed in units of ppm on a water-free basis
% CO2Index - Index of the species in the gas mixture Cantera objects
% COAL - flow rate of pulverized coal through burner (kg/hr)
% COAL_Type - 1, 2, or 3 to determine which char reaction parameters will
% be used:
% 1 = Wyoming Sub-bituminous
% 2 = Illinois #6
% 3 = Pittsburgh #8
% composition - string containing details of initial gas composition
% convheat - convective heat transferred to particles (CPD model)
% count - counts iterations to get CSTR to steady state
% cp - particle heat capacity (J/kg-K and cal/g-K depending on context)
% cpa - ash heat capacity (cal/g-K)
% cpc - DAF char heat capacity (cal/g-K)
% cpgvector - mass-based heat capacity of gas (J/kg-K)
% CPU_time - stores the total time CPU is used in ignition part of model
% CPU_time_per_step - store the average time per iteration that the CPU is
% used during the ignition part of the model
% cstr - a Cantera reactor object
% Current_position - variable used in linear interpolation of wall
% temperature data
% d - reactor diameter (m)
% delhv - heat of pyrolysis (cal/g), negative indicates endothermic,
% nominally -100.0 cal/g
% delta_max - the maximum temperature change between iterations
% delta_T - change in temperature from 1 iteration
% devol_incomplete - logical variable that tracks whether devolatilization
% is complete or not
% DiffCoeffs - mixture-averaged diffusion coefficents (m^2/s)
% diffwvector - mixture-averaged diffusion coefficient (m^2/s) for water
% downstream - a Cantera reservoir object connected to the outlet of each
% CSTR
% dp - coal particle diameter (only 1 value is used to represent all
% particles)
% dpout - coal particle diameter from CPD model
% dt - time step for CSTR network integration (s)
% dummygas - a gas mixture object used for property evaluation at film
% temperatures
% EACO2 - activation energy for char gasification by CO2 (cal/gmole)
% EAO2 - activation energy for char oxidation (cal/gmole)
% emiss - emissivity (gray) values for radiation calculations as defined in
% input file comments
% env - a reservoir that represents the ambient environment that absorbs
% heat transferred through the reactor external walls
% exhaust - Cantera gas object used to create constant pressure reservoirs
% downstream of each CSTR
% ExistingGas - vector of mass flow rates of species in the gas mixture
% (kg/s)
% F - array of radiation view factors
% fchar - equals 1 - fvol. See cpdcp_nlg for more details
% fcross - fraction of original D.A.F. coal that was metaplast and
% crosslinked into the char matrix
% ffgas the fraction of total mass release (i.e. volatiles) that is
% h2o, co2, ch4, co, and other light gases.
% ffgas2 - ffgas modified to include unknown volatiles estimated as a
% mixture of CH4 and C2H2
% fgas - mass fraction of D.A.F. coal evolved as light gas
% fid - pointer to output file
% fmet - mass fraction of D.A.F. coal existing as metaplast
% fnchar - fraction of original nitrogen remaining in char (and metaplast)
% fnhcn - fraction of original nitrogen released as light gas (by
% difference: 1 - fnchar - fntar)
% fnt - nitrogen content of char and metaplast
% fntar - fraction of original nitrogen released as tar
% fntot - total fractional release of nitrogen
% ftar - mass fraction of D.A.F. coal evolved as tar
% fvol - mass fraction of D.A.F. coal evolved as volatiles (light gas +
% tar)
% gas - a gas mixture object
% gas_const - The universal gas constant in cal/gmol.K
% gasification - 1 or 0 values respectively determine whether char
% gasification by CO2 is modeled or not
% gw - a Cantera wall object between adjacent CSTR's to model conduction
% heat transfer between the gases in each CSTR
% h - convective heat transfer coefficient
% H - enthalpy of gas mixture (J/kg)
% H2OIndex - Index of the species in the gas mixture Cantera objects
% HCNIndex - Index of the species in the gas mixture Cantera objects
% Hdiff - the difference between the hydrogen release predicted by CPD and
% the amount of hydrogen release expected from correlation
% with total mass loss
% Hgoal - the target amount of hydrogen release from volatiles correlated
% with total mass release
% H_in - enthalpy from heterogeneous reactions and CO from the same (J/s)
% H_new - enthalpy of gas phase mixture after heterogeneous reactions (J/s)
% or after mixing of burnout and primary streamss
% H_nlg - the hydrogen release predicted by the CPD model
% H old - enthalpy of existing gas phase mixture before char reactions
% are calculated (J/s)
% H_out - enthalpy leaving the gas phase with O2 and CO2 consumed by the
% char (J/s)
% H_Primary - enthalpy of primary reactants just before burnout oxidizer
% is added (J/s)
% H_Burnout - enthalpy of burnout oxidizer stream (J/s)
% i - index used in loops
% iterate - a logical variable that causes iteration to continue while true
% j - index used in loops
% k_gas - gas thermal conductivity (W/m.K)
% k_wall - thermal conductivity for reactor wall (W/m.K) Value is empirical
% and linked to the value of Length and Length2
% kgvector - vector of gas thermal conductivity (W/m.K)
% kpCO2 - rate coefficient for CO2 gasification of char (g of oxidizer
% consumed per second per unit char surface area (in cm^2)
% per atm of oxidizer)
% kpO2 - rate coefficient for oxidation of char (g of oxidizer consumed per
% second per unit char surface area (in cm^2) per atm of
% oxidizer)
% Length - axial length of CSTR's in ignition section of model
% Length2 - axial length of CSTR's in post-ignition sections of model
% maingas - Cantera gas mixture object used for post-ignition model
% MassProportionCH4 - the mass proportion of CH4 in CH4 and C2H2 mixture
% that makes up the volatiles not predicted by CPD model
% max_index - number of points in linear iterpolation of wall temperature
% data
% max_position - the last axial position point in the wall temperature data
% mdel - average molecular weight per side chain (C13 NMR parameter)
% mechanism - 1, 2, or 3 to select from the available gas-phase mechanisms:
% 1 = GRI-Mech 3.0
% 2 = GRI-Mech 3.0 + B96 (includes advanced reburning)
% 3 = SKG03
% mfc - mass flow controller Cantera object
% mixture - vector of mass flow rates of species in the gas mixture passed
% to the next CSTR (kg/s)
% m_dot - gas phase mass flow (kg/s)
% MolarProportionCH4 - the molar proportion of CH4 in CH4 and C2H2 mixture
% that makes up the volatiles not predicted by CPD model
% mw1 - average molecular weight per aromatic cluster (C13 NMR parameter)
% mwchar - like mw1, but for char. Initially it is set to mw1
% MWmix - mixture molecular weight
% MWunknowns - the mixture molecular weight of the volatiles not predicted
% by CPD model
% n - iteration counter
% N2 - moles/hr N2 in the reactants
% N2_in - flow rate of bottled N2 (not air N2) through the burner (kg/hr)
% N2Index - Index of the species in the gas mixture Cantera objects
% NCE - conversion efficiency of fuel nitrogen converted to NO + NO2
% assuming all NO + NO2 originates from fuel nitrogen
% network - a Cantera object that holds the reactor networks
% network_cell_array - an array of CSTR reactor networks. Each network
% consists of the upstream reservoir, followed by the mass
% flow controller, the CSTR, the valve and the downstream
% reservoir
% new_T - temperature of CSTR after an iteration
% NG - moles/hr natural gas (CH4)
% NG_in - flow rate of natural gas through burner (kg/hr)
% nj - intermediate variable in calculation of r
% NO - moles/hr NO in the CO2 reactant stream
% NO_doping - ppm NO in the CO2 reactant streams
% NOIndex - Index of the species in the gas mixture Cantera objects
% NO2Index - Index of the species in the gas mixture Cantera objects
% NOx_ppm_dry - NO + NO2 expressed in units of ppm on a water-free basis
% nu - Nusselt number
% number_reactors - number of CSTR's in ignition section of model
% O2 - moles/hr O2 (from all reactants - air and bottled O2)
% O2_in - flow rate of bottled O2 (not air O2) through the burner (kg/hr)
% O2_vol_dry - O2 expressed in units of vol% on a water-free basis
% O2Index - Index of the species in the gas mixture Cantera objects
% old_T - temperature of CSTR before iteration
% omegaa - mass fraction of ash in the parent coal (as received)
% omegaw - mass fraction of moisture in the parent coal (as received)
% output - string containing name of output file
% P - pressure (Pa)
% p0 - ratio of bridges to total attachments (C13 NMR parameter)
% particles - number of coal particles per second
% position - vector of CSTR locations along the MFR axis
% press - pressure in atmospheres
% pressurevector - vector of CSTR pressures (Pa)
% prgas - Prandtl number of gas
% Primary - fraction of total oxidizer through the burner (0 < Primary < 1)
% prompt - multiplier for prompt NOx mechanism reactions (0 or 1)
% qconv - convective heat transer rate (cal/s)
% Q_conv - qconv converted to J/s
% qrad - radiation heat transfer rate (cal/s)
% Q_rad - qrad converted to J/s
% Q_react - rate of heat from heterogeneous reactions (J/s)
% Q_reactO2_x - Fraction (0 to 1) of heterogeneous O2 reaction heat to char
% Nominally 0 because 0.5 and 1.0 gave problems in testing
% Q_reactCO2_x - as for Q_reactO2_x, but for CO2 gasification
% Q_total - total heat rate to char from heterogeneous reactions, and heat
% transfer (radiatve and convective)
% r - chemical equivalence ratio (r>1 means fuel-rich, r=1 is
% stoichiometric, r<1 means fuel-lean)
% R - thermal resistance (calculated from k_wall or k_gas depending on
% context
% rhoCHAR - density of char (including ash) (kg/m^3)
% rhogas - gas density (kg/m^3)
% rhop - initial particle apparent density (g/cm^3) - see additional
% comments in input script
% rpoCO2 - rate of CO2 consumption by CO2 gasification (kg/s)
% rpoO2 - rate of O2 consumption by oxidation (kg/s)
% rtot - total rate of char consumption (kg/s) by char oxidation and
% gasification
% sigma - a radiation constant (cal/s cm^2 K^4)
% sigp1 - sigma + 1 is the number of total attachments per cluster (C13 NMR
% parameter)
% swell - swelling factor (dpf/dp0 - 1) where dpf = final/max diameter and
% dp0 = initial diameter. See additional comments in the
% input script
% t0 - the time at the start of iteration
% T1 - initial gaseous reactant temperature (K) (through the burner)
% T2 - temperature of preheated burnout oxidizer (K)
% tbnr - burner face temperature (K)
% Temp - temperatures used in radiation heat transfer calculations
% texit - exhaust pipe temperature (K)
% tfilm - film temperature (average of gas and particle temperatures) (K)
% tg - gas temperature (K)
% tgas - gas temperature (K)
% thermal - multiplier for thermal NOx mechanism reactions (0 or 1)
% timax - maximum allowable devolatilization time for CPD model
% timedata - time of day that output file is started
% tme - cumulative integration time for reactor networks
% tms - time in milliseconds from CPD model
% tolerance - maximum allowable change in temperature between time steps of
% CSTR integration to force steady state conditions
% tp - particle temperature (K)
% TR1 - initial guess of temperature for CSTR's in ignition section (K)
% trate - particle heating rate (K/s)
% twall - linearly interpolated wall temperatures
% twallvector - wall temperatures used as an input (K) - see WallX
% ugvector - vector of gas viscosity (g/(m.s) = Pa.s)
% upstream - reservoir upstream of a CSTR containing the gas mixture
% entering that CSTR
% v - a Cantera valve object between a CSTR and its downstream reservoir -
% set to maintain constant pressure in the CSTR
% V_m - intermediate variable in calculation of r
% V_minus - negative Oxidation States for C,H,N,& O in that order
% V_p - intermediate variable in calculation of r
% V_plus - positive oxidation states for C,H,N,& O in that order
% velocity - gas (and particle) velocity (m/s)
% Volume - volume of CSTR's (m^3)
% w - a Cantera wall object installed between the cstr's and the
% environment, env
% WallX - axial locations of wall temperatures used as an input (m)
% water - water content of particles
% XCO - mole fraction of CO
% XCO2 - mole fraction of CO2
% Xi - mole fractions of species in gas mixture
% xm - Particle position (m) from CPD model
% XNO - mole fraction of NO
% XNO2 - mole fraction of NO2
% XO2 - mole fraction of O2
% XWater - mole fraction of water
% xwbvector - vector of bulk gas water concentration (units: mole fraction)
% y - gas temperature (used in calculating temperature changes during
% iteration
% YCO - mass fraction of CO
% YCO2 - mass fraction of CO2
% yelem - DAF mass fractions of CHNOS for the coal in a 5 element vector,
% each element between 0 and 1
% yf - A CPD indicator of the fraction of total light gas that has been
% released. The look up table on light gas composition is
% based on yf. Called Xgas in Genetti's MS thesis - see the
% thesis for more detail.
% Yi - mass fractions of gas mixtures in upstream reservoirs
% YNO - mass fraction of NO
% YNO2 - mass fraction of NO2
% yNsite - variable from CPD model (undocumented)
% YO2 - mass fraction of O2
% Ywater - mass fraction of water in upstream reservoir
% yygas - the fractions of light gas release that is h20, co2, ch4, co and
% other light gases.
% z - distance from burner (cm)