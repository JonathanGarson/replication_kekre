%% ident_Delta.m
% Characterizes identification of Delta

%% Preliminaries 

% Get the number of cores allocated to your job
num_cores = str2num(getenv('SLURM_CPUS_PER_TASK'));

% Start a pool of workers
c = parcluster;
poolobj = parpool(c, num_cores);

%% Now my code...

masterTime = tic;

%% Name parameter and range to be varied (only part that needs to change across ident codes)

varyParam1 = 'Delta'; 
varyRange1 = [0.0033 0.0037 0.0041 0.0045 0.0049];
varyParam2 = 'beta'; % set to nan if only one parameter needs to change
varyRange2 = 0.99335 - (varyRange1-0.0045);
numWorkers = length(varyRange1);

%% Initialize economic parameters

% Tastes
beta = 0.99335;
numBetas = 3;
Delta = 0.0045;
xi = 15; 
sigma = 1;

% Wages
phi = 0.956582594041537;
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
k = 0.044375024348784;
deltaBar = 0.034;
eps_delta_a = -0.011;
eps_delta_beta = -4.55;
mbar = 0.188514800460487; 
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
znext_g = -8.311708508271963;

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

sf = [0.017122685989587
   0.021603185317716
   0.025955422085097
   0.051462689972529
   0.062991467099531
   0.075210602430270
   0.138920127483958
   0.167867614015335
   0.200577236609508
   0.022089357506147
   0.026518038946109
   0.030851544593165
   0.054494843780907
   0.064236860936854
   0.075578981932043
   0.092978861699256
   0.114667882201071
   0.145072709854387
   0.046351305728673
   0.049730604428140
   0.053687506404748
   0.063791631890200
   0.070911708836071
   0.081422516403747
   0.020668168438229
   0.041905858473938
   0.073981567227034]';

starting = [0.998336106489185,0.043815722088853,0.962081725622082,sf];

gridPoints = 151;
eps_goods = 0.005;  
eps_t = 0.0001; 
eps_vac = 0.0005;
eps_bargain = 0.0005; % larger values of eps since simply characterizing
 % identification, not stationary RCE precisely
upper = 4000;
eps_v = 0.0001; 
finerStatGrid = 0;
eps_phi = 0.0000001;

forceEnd = 18*60*60;
accuracy = struct('gridPoints',gridPoints,...
    'eps',[eps_goods eps_t eps_vac eps_bargain],'upper',upper,'eps_v',eps_v,...
    'finerStatGrid',finerStatGrid,'eps_phi',eps_phi,'forceEnd',forceEnd);

%% Parallelize solution

multParams = ~isnan(varyParam2);
parfor i=1:numWorkers
    thisParams = params;
    thisVaryRange1 = varyRange1(i);
    thisVaryRange2 = varyRange2(i);
    
    thisParams = setfield(thisParams,varyParam1,thisVaryRange1);
    if multParams
        thisParams = setfield(thisParams,varyParam2,thisVaryRange2);
    end
    
    % Set parameters which might have been updated
    thisBetas = linspace(thisParams.beta-thisParams.Delta,thisParams.beta+thisParams.Delta,numBetas);
    thisParams = setfield(thisParams,'betas',thisBetas);    
    
    try
        rng(1,'twister');
        [m,t,theta,sf] = StationaryRCE(thisParams,accuracy,starting,...
            strcat('output/ident_',varyParam1,'_',num2str(i)),i,1);
        format long;
        disp(['WORKER ',num2str(i),': m=',num2str(m),' t=',num2str(t),' theta=',num2str(theta)]);
        for i_e=1:length(sf)
            disp(['WORKER ',num2str(i),': sf(',num2str(i_e),') = ',num2str(sf(i_e))]);
        end   
    catch e
        warning(strcat('ERROR ON WORKER ',num2str(i),'. SHUTTING THIS WORKER DOWN.'));
        disp(getReport(e));
    end
end

disp(['Parallel code completed in ',num2str(toc(masterTime)),' s.']);

%% Shut down pool of workers
 delete(poolobj)