%% solveStationaryRCE_durElast.m
% Solves stationary RCE with high calibrated duration elasticity to UI

masterTime = tic;

%% Initialize economic parameters

% Tastes
beta = 0.99335;
numBetas = 3;
Delta = 0.0045;
betas = linspace(beta-Delta,beta+Delta,numBetas);
xi = 2.5; 
sigma = 1;

% Wages
phi = 0.969004268132240;
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
k = 0.044325280956058;
deltaBar = 0.034;
eps_delta_a = -0.011;
eps_delta_beta = -4.55;
mbar = 0.226263867444381;
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
znext_g = -8.483032016374295;

% Equity share cutoff (only used to construct moments regarding disposable income)
equitycutoff = 0.35;

% Retailers
epsilon = 6;
tauR = -0.076;

% Max UI
uimax = 0.4;

% Put all of the above in 'params'
varList = whos;
params = struct;
for ii = 1:length(varList)
    params.(varList(ii).name) = eval(varList(ii).name);
end

%% Initialize other program parameters

sf = [0.015638869229244
   0.019371586136749
   0.022730616760702
   0.051156954732938
   0.061039164965472
   0.070509232323881
   0.147574017255790
   0.172592495056929
   0.197726135059687
   0.019077730990960
   0.022641578256485
   0.025927726894160
   0.055065251090416
   0.063105129883856
   0.071670232053665
   0.120604221572446
   0.137402204338589
   0.159072185858231
   0.039308772913698
   0.041774833101107
   0.044591022475637
   0.072004169146561
   0.076898890765074
   0.084166529164706
   0.053578851441080
   0.068389269683463
   0.090604543644436]';

starting = [0.998336106489185,0.046538948244988,1.159058643101021,sf];

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
[m,t,theta,sf] = StationaryRCE(params,accuracy,starting,'output/RCE_durElast',0,0);

disp(['Code completed in ',num2str(toc(masterTime)),' s.']);