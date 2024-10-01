%% simGR_seprate_time1.m
% Given separation rate shocks, calibrates one period of discount factor 
%  shock to match unemployment in 05/2008 and construct Jacobian.  
%  Assumes iota = 0.975 for concreteness.

masterTime = tic;

% Load data to be replicated
urdata = importdata('output/ur05081214.csv');

% Load separation rate shocks
deltaBarShocks = importdata('output/seprate05081214.csv');
deltaBarPersist = 0.489132300622970;

% First set znextDiff, goodsDiff, tDiff, equityValDiff, vacpostDiff, and 
%  dwageDiff implied by computational error in steady-state
load(strcat('output/trans_baseline_noshocks_sticky.mat'),'znextTrans','zsupplyTrans',...
    'cTrans','prodTrans','tTrans','t','mTrans','m','equityValTrans','equityVal',...
    'vacpostTrans','sfNewHiresTrans','wTrans','w');
znextDiff = znextTrans - zsupplyTrans;
goodsDiff = cTrans - prodTrans;
tDiff = tTrans - t;
mDiff = mTrans - m;
equityValDiff = equityValTrans - equityVal;
dvacpostDiff = vacpostTrans - sfNewHiresTrans;
dwageDiff = wTrans-w;
cTransSS = cTrans;

clearvars znextTrans zsupplyTrans cTrans prodTrans tTrans t ...
    mTrans m equityValTrans equityVal vacpostTrans sfNewHiresTrans wTrans w;

load('output/RCE_baseline','m','p_e','eStates','uStates',...
    'numBetas','xi','sigma','numA','numP','rhoP','sigmaP','numT','sigmaT',...
    'numUI','Tua','dbarbar',...
    'k','eta','pEpsilon','phi',...
    'betas','Ta','Taui','zeta','dbar','rr','deltaBar',...
    'eps_delta_a','eps_delta_beta','mbar','lambda',...
    'abar','a','aP','zbar','znext_g','omega0','omega1','omega2',...
    'qStat','z','p_e_tilde','phi_e_tilde','phi_u_tilde','gridpointsStat',...
    'dbarbar','theta','t','w','tauR','uimax','RCEName');

load('output/jacobian_baseline_1.mat');
H1 = H;
clearvars H;
load('output/jacobian_baseline_2.mat');
H2 = H;
clearvars H;
load('output/jacobian_baseline_3.mat');
H3 = H;
clearvars H;
load('output/jacobian_baseline_4.mat');
H4 = H;
clearvars H;
H = [H1,H2,H3,H4];
Tbar = 481;
% correct for iota = 0.975 
H(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),:) = ...
    (1-0.975)*H(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),:);
H(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),4*Tbar-2:(4+eStates)*(Tbar-1)+1) = ...
    H(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),4*Tbar-2:(4+eStates)*(Tbar-1)+1) - ... 
    0.975*eye(eStates*(Tbar-1),eStates*(Tbar-1));
Hss = H;

invariantParams = ...
    struct('numBetas',numBetas,'xi',xi,'sigma',sigma,...
    'numA',numA,'numP',numP,'rhoP',rhoP,'sigmaP',sigmaP,'numT',numT,'sigmaT',sigmaT,...
    'numUI',numUI,'Tua',Tua,'dbarbar',dbarbar,...
    'k',k,'eta',eta,'pEpsilon',pEpsilon,'phi',phi,'tauR',tauR);

betasTrans = repmat(betas',1,Tbar);
TaTrans = repmat(Ta,1,1,Tbar);
TauiTrans = repmat(Taui,1,1,Tbar);
zetaTrans = zeta*ones(1,Tbar);
dbarTrans = dbar*ones(1,Tbar);
rrTrans = rr*ones(1,Tbar);
deltaBarTrans = repmat(deltaBar',1,Tbar);
eps_delta_aTrans = repmat(eps_delta_a,1,Tbar);
eps_delta_betaTrans = repmat(eps_delta_beta,1,Tbar);
mbarTrans = repmat(mbar,1,Tbar);
lambdaTrans = lambda*ones(1,Tbar);
abarTrans = repmat(abar,1,Tbar);
aTrans = repmat(a',1,Tbar);
aPTrans = repmat(aP',1,Tbar);
zbarTrans = zbar*ones(1,Tbar);
znext_gTrans = znext_g*ones(1,Tbar);
omega0Trans = repmat(omega0,1,Tbar);
omega1Trans = repmat(omega1,1,Tbar);
omega2Trans = repmat(omega2,1,Tbar);
uimaxTrans = repmat(uimax,1,Tbar);
g = zeros(1,Tbar);

varyingParams = ...
    struct('Tbar',Tbar,'betas',betasTrans,...
    'Ta',TaTrans,'Taui',TauiTrans,...
    'zeta',zetaTrans,'dbar',dbarTrans,'rr',rrTrans,'deltaBar',deltaBarTrans,...
    'eps_delta_a',eps_delta_aTrans,'eps_delta_beta',eps_delta_betaTrans,...
    'mbar',mbarTrans,'lambda',lambdaTrans,...
    'abar',abarTrans,'a',aTrans,'aP',aPTrans,'zbar',zbarTrans,'znext_g',znext_gTrans,...
    'omega0',omega0Trans,'omega1',omega1Trans,'omega2',omega2Trans,'uimax',uimaxTrans,'g',g);  

[equityShares] = ...
    initializeEquity(0.35,...
    (1/m)*qStat,z,[p_e_tilde*phi_e_tilde,(1-p_e_tilde)*phi_u_tilde]);
equityShares_e = repmat(equityShares,1,eStates);
equityShares_u = repmat(equityShares,1,uStates);

% Characterize equilibria
epsSticky = [0.000001 0.000001 0.000001 0.00001 0.000001 0.000001];
epsMacro = 0.001;

transParamsSticky = ...
    struct('psi',360,'adjResource',0,'epsilon',6,'phiInt',1/m-1,'phiPi',1.5,'phiC',1/12,'phiRho',0,'zlb',1,'fixedM',0,'fixedm',0,...
    'iota',0.975,'znext_gT1',0,'rhoZnext_g',nan,'znext_gT2',nan,'fixt',nan,...
    'forceEnd',72000);

% Generate duration of UI in Great Recession
[~,dbarPathMonths,dbarShockMonths] = genUIGreatRecession(Tbar);

try
    GR(Tbar,1,1,8*5.2500e-05,0.95,...
        deltaBarShocks,deltaBarPersist,dbarPathMonths,dbarShockMonths,Hss,[],urdata,0,...
        transParamsSticky,...
        RCEName,'output/GR_seprate_time1',invariantParams,varyingParams,...
        equityShares_e,equityShares_u,epsSticky,epsMacro,znextDiff,goodsDiff,...
        tDiff,mDiff,equityValDiff,dvacpostDiff,(1-0.975)*dwageDiff,cTransSS,...
        0,[],[]);             
catch err
    disp(['Error: ',err.identifier]);
end

disp(['Code completed in ',num2str(toc(masterTime)),' s.']);