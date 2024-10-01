%% solveStationaryRCE_noEps.m
% Solves stationary RCE without separation rate heterogeneity

masterTime = tic;

%% Initialize economic parameters

% Tastes
beta = 0.99337;
numBetas = 3;
Delta = 0.0045;
betas = linspace(beta-Delta,beta+Delta,numBetas);
xi = 15; 
sigma = 1;

% Wages
phi = 0.955274414328896;
rhoP = (0.9695)^(1/12);
sigmaP = sqrt((1-rhoP^2)*(0.0384/(1-(0.9695)^2))); 
numP = 3;
sigmaT = sqrt(0.0522);
numT = 3;

% UI
zeta = 0.5;
dbarbar = 9;
dbar = 6;
rr = 0.5;

% Technology
abar = 1;
k = 0.044589452479389;
deltaBar = 0.034;
eps_delta_a = 0;
eps_delta_beta = 0;
mbar = 0.195494134376014;
eta = 0.7;
lambda = -0.14;

% Income process
omega0 = 0.33;
omega1 = 0.37;
omega2 = 0.50;

[a,aP,numA,Ta,numUI,Taui,Tua,pEpsilon] = genIncomeProcess(abar,...
    rhoP,sigmaP,numP,sigmaT,numT,zeta);

% Assets 
zbar = -0.15; 
znext_g = -8.980553650372542;

% Equity share cutoff (only used to construct moments regarding disposable income)
equitycutoff = 0.35;

% Retailers
epsilon = 6;
tauR = -0.076;

% Max UI
uimax = 0.51;

% Put all of the above in 'params'
varList = whos;
params = struct;
for ii = 1:length(varList)
    params.(varList(ii).name) = eval(varList(ii).name);
end

%% Initialize other program parameters

sf = [0.016148507290193
   0.020844761445035
   0.025391005254216
   0.050523845318802
   0.062819764679651
   0.075538400141736
   0.131018264209080
   0.161104057947118
   0.194867420108101
   0.020451307490297
   0.024984041054341
   0.029437444923104
   0.051918277044137
   0.061735240547047
   0.073372803351997
   0.080921599924647
   0.102839867738469
   0.134460993640353
   0.043478018428846
   0.046832071959979
   0.050863133335317
   0.059716108024636
   0.067048776516675
   0.077875858263745
   0.020510354821502
   0.033892778947078
   0.066896966984023]';

starting = [0.998336106489185,0.053419932760375,0.981797209109838,sf];

gridPoints = 151;
eps_goods = 0.001;
eps_t = 0.00001;
eps_vac = 0.00005;
eps_bargain = 0.00005;
upper = 4000;
eps_v = 0.0001;
finerStatGrid = 0;
eps_phi = 0.0000001;

forceEnd = 18*60*60;
accuracy = struct('gridPoints',gridPoints,...
    'eps',[eps_goods eps_t eps_vac eps_bargain],'upper',upper,'eps_v',eps_v,...
    'finerStatGrid',finerStatGrid,'eps_phi',eps_phi,'forceEnd',forceEnd);

rng(1,'twister');
[m,t,theta,sf] = StationaryRCE(params,accuracy,starting,'output/RCE_noEps',0,0);

disp(['Code completed in ',num2str(toc(masterTime)),' s.']);