%% solveStationaryRCE_liquid.m
% Solves stationary RCE calibrated to liquid wealth

masterTime = tic;

%% Initialize economic parameters

% Tastes
beta = 0.99425;
numBetas = 3;
Delta = 0.00125;
betas = linspace(beta-Delta,beta+Delta,numBetas);
xi = 15; 
sigma = 1;

% Wages
phi = 0.952240395125342;
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
k = 0.044441655486790;
deltaBar = 0.034;
eps_delta_a = -0.011;
eps_delta_beta = -14;
mbar = 0.189175943547784;
eta = 0.7;
lambda = -0.14;

% Income process
omega0 = 0.33;
omega1 = 0.37;
omega2 = 0.50;

[a,aP,numA,Ta,numUI,Taui,Tua,pEpsilon] = genIncomeProcess(abar,...
    rhoP,sigmaP,numP,sigmaT,numT,zeta);

% Assets 
zbar = -0.5; 
znext_g = -5.100341360459492;

% Equity share cutoff (only used to construct moments regarding disposable income)
equitycutoff = 0.75;

% Retailers
epsilon = 6;
tauR = -1/epsilon;

% Max UI
uimax = 0.44;

% Put all of the above in 'params'
varList = whos;
params = struct;
for ii = 1:length(varList)
    params.(varList(ii).name) = eval(varList(ii).name);
end

%% Initialize other program parameters

sf = [0.027289358749948
   0.032711152353079
   0.038177831434849
   0.066181774937706
   0.077885979253036
   0.091561659908258
   0.119375921663903
   0.146117879495176
   0.183664442348864
   0.030462609629454
   0.036008245634438
   0.041526030006337
   0.068391172793827
   0.079634033034764
   0.093063163208361
   0.086895855028217
   0.112570653412918
   0.150387209061454
   0.034798842273104
   0.040442032848710
   0.045998176051825
   0.070455998625022
   0.080883153703086
   0.093982054975328
   0.044184858357916
   0.070542664076967
   0.109789409663484]';

starting = [0.998336106489185,0.042013957602298,0.960639281601245,sf];

gridPoints = 101;
eps_goods = 0.001; 
eps_t = 0.00001;
eps_vac = 0.00005;
eps_bargain = 0.00005;
upper = 1000;
eps_v = 0.0001; 
finerStatGrid = 0;
eps_phi = 0.0000001;

forceEnd = 18*60*60;
accuracy = struct('gridPoints',gridPoints,...
    'eps',[eps_goods eps_t eps_vac eps_bargain],'upper',upper,'eps_v',eps_v,...
    'finerStatGrid',finerStatGrid,'eps_phi',eps_phi,'forceEnd',forceEnd);

rng(1,'twister');
[m,t,theta,sf] = StationaryRCE(params,accuracy,starting,'output/RCE_liquid',0,0);

disp(['Code completed in ',num2str(toc(masterTime)),' s.']);