function param = parameters(varargin)
%param = parameters;
%param = parameters('name',value)
%
%This function creates a structure containing all the required parameters
%to solve the crack growth and serpentinisation kinetics problems.

p = inputParser;

p.addParamValue('T',300+273.15, @isnumeric);
p.addParamValue('P',50e6, @isnumeric);
p.addParamValue('spinf',0,@isnumeric);
p.addParamValue('Eprime',2.08e11,@isnumeric);
p.addParamValue('hsgrain',100e-6,@isnumeric);
p.addParamValue('hsgrainmin',10e-6,@isnumeric);
p.addParamValue('ac0',0.1,@isnumeric);
p.addParamValue('Gc',2,@isnumeric);
p.addParamValue('grainsizeprop',2,@isnumeric);
p.addParamValue('supcrtfile','SUPCRT/data_P',@ischar);
p.addParamValue('Nfold',8,@isnumeric);

p.parse(varargin{:});
param = p.Results;

%set constants:
param.theta = 2*1.454;
param.mu = 1.1215*pi/2;

% parameters for temperature dependence from Malvoisin et al. (2012)
bc = 3640;
cc = 8759;
T0 = 623.6;
Ac = 808.3;
 % temperature dependence from Malvoisin et al. (2012)
k = 1e-13*Ac*exp(-bc/param.T)*(1-exp(-cc*(1/param.T-1/T0)));
% molar volume of serpentine
Vm= 1.073E-04;
% gas constant
R =8.314E+00;

% reaction rate in m/s
Kl = k/(R*param.T)*Vm;

% crystallization pressure

%get thermodynamic parameters from the supcrt file:
A = textread(param.supcrtfile);
P_vec = A(:,1)/10*1e6;    % pressure in megapascal
T_vec = A(:,2)+273.15;    % temperature in K
omega_vec = 10.^(A(:,4)); % supersaturation index
sc_vec = R.*T_vec.*log(omega_vec)./Vm; % crystallisation pressure

%make interpolation table:
T_vec  = unique(T_vec);
P_vec  = unique(P_vec);
sc_tab = reshape(sc_vec,length(T_vec),length(P_vec));

%max crystallisation stress:
sc = interp2(P_vec,T_vec,sc_tab,param.P,param.T);

% minimum length for crack propagation
cmin = param.Gc*pi*param.Eprime/(4*(gam(1)*sc-param.mu*param.spinf)^2);

% characteristic length
l   = pi/4*param.Eprime/param.mu^2*param.Gc/(sc^2);
% characterisitc time
tau = l/sc/Kl;

param.k = k;
param.Vm_liz = Vm;
param.Vm_brc =24.7e-6;

%molar volume of forsterite
param.Vm_for = 44.66e-6; %m^3/mol

param.sc = sc;
param.cmin = cmin;
param.Kl = Kl;
param.l = l;
param.tau = tau;

% normalised parameters:
param.Enorm = param.Eprime/sc;
param.sigmainf = param.spinf/sc;
param.cminnorm = (gam(1)/param.mu-param.sigmainf)^(-2);
param.hsgrainnorm = param.hsgrain/l;
param.hsgrainminnorm = param.hsgrainmin/l;

%normalised initial conditions
param.c0 = param.cminnorm*1.2;
param.a0 = param.c0*param.ac0;
param.w0 = param.c0/param.Enorm *param.theta*param.sigmainf; %sn=0

