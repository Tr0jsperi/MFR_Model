%---LIST OF __VARIABLES__ (alphabetical)---------------------------------------
% Note that MATLAB is case-sensitive.
%
% Units specified here are as at the first instance of the variable in the
% code (such as the user-input values). Unit conversions are noted in the
% code when they are changed.
%
A - array of areas used in radiation heat transfer calculations (cm^2)
A_char - surface area of char (cm^2)
ACO2 - pre-exponential factor for char gasification by CO2
% (g/(cm2s atmCO2)

AddedGas - vector of mass flow rates of species to be added to the gas
% mixture due to coal devolatilization and heterogeneous
% reactions (kg/s)
aij - an array of the number of atoms of C, H, N, and O in the species in
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

ButtonName - yes, no, or cancel result of asking user whether to continue

% iterations in ignition section of code
% c - vector of empirical coefficients used when estimating C13 NMR
% parameters from proximate and ultimate analysis results
% c0 - char bridge population (C13 NMR parameter)
% C2H2Index - Index of the species in the gas mixture Cantera objects

Cdiff - the difference between the carbon release predicted by CPD and the amount of carbon release expected from correlation with total mass loss

Cgoal - the target amount of carbon release from volatiles correlated with total mass release

% CH4Index - Index of the species in the gas mixture Cantera objects
% CHAR - mass flux of DAF CHAR (kg/s)
% chargas - Cantera gas mixture object used for char reactions section of
% model
% CHNOIndex - indices for the elements in the Cantera objects
% CHratio - the CH ratio of the volatiles not predicted as specific species
% by the CPD model

C_nlg - the carbon release predicted by the CPD model

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

cpgvector - mass-based heat capacity of gas (J/kg-K)

% CPU_time - stores the total time CPU is used in ignition part of model
% CPU_time_per_step - store the average time per iteration that the CPU is
% used during the ignition part of the model
% cstr - a Cantera reactor object
% Current_position - variable used in linear interpolation of wall
% temperature data
% d - reactor diameter (m)
% delhv - heat of pyrolysis (cal/g), negative indicates endothermic,
% nominally -100.0 cal/g

delta_max - the maximum temperature change between iteration

% delta_T - change in temperature from 1 iteration
% devol_incomplete - logical variable that tracks whether devolatilization
% is complete or not

DiffCoeffs - mixture-averaged diffusion coefficents (m^2/s)

diffwvector - mixture-averaged diffusion coefficient (m^2/s) for water

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

ffgas - the fraction of total mass release (i.e. volatiles) that is h2o, co2, ch4, co, and other light gases.

ffgas2 - ffgas modified to include unknown volatiles estimated as a mixture of CH4 and C2H2


% fgas - mass fraction of D.A.F. coal evolved as light gas
% fid - pointer to output file

flagtime - time to judge if the advance function works successfully

fmet - mass fraction of D.A.F. coal existing as metaplast

fnchar - fraction of original nitrogen remaining in char (and metaplast)

fnhcn - fraction of original nitrogen released as light gas (by difference: 1 - fnchar - fntar)

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

gw - a Cantera wall object between adjacent CSTR's to model conduction

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

kgvector - vector of gas thermal conductivity (W/m.K)

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

network - a Cantera object that holds the reactor networks

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

prgas - Prandtl number of gas

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

rhogas - gas density (kg/m^3)

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

tg - gas temperature (K)

tgas - gas temperature (K)

% thermal - multiplier for thermal NOx mechanism reactions (0 or 1)
% timax - maximum allowable devolatilization time for CPD model
% timedata - time of day that output file is started
% tme - cumulative integration time for reactor networks

tms - time in milliseconds from CPD model

% tolerance - maximum allowable change in temperature between time steps of
% CSTR integration to force steady state conditions

tp - particle temperature (K)

% TR1 - initial guess of temperature for CSTR's in ignition section (K)
% trate - particle heating rate (K/s)
% twall - linearly interpolated wall temperatures
% twallvector - wall temperatures used as an input (K) - see WallX

ugvector - vector of gas viscosity (g/(m.s) = Pa.s)

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

xm - Particle position (m) from CPD model

% XNO - mole fraction of NO
% XNO2 - mole fraction of NO2
% XO2 - mole fraction of O2
% XWater - mole fraction of water

xwbvector - vector of bulk gas water concentration (units: mole fraction) 

y - gas temperature (used in calculating temperature changes during iteration

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
