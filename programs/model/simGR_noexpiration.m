%% simGR_noexpiration.m
% Simulates discount factor shocks calibrated in simGR_expiration
%  without UI expiration in December 2013

masterTime = tic;

% Load data to be replicated
urdata = importdata('output/ur05081214.csv');

% First set znextDiff, goodsDiff, tDiff, equityValDiff, vacpostDiff, and 
%  dwageDiff implied by computational error in steady-state
load(strcat('output/trans_baseline_dbarbar22_noshocks_sticky.mat'),'znextTrans','zsupplyTrans',...
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

load('output/RCE_baseline_dbarbar22','m','p_e','eStates','uStates',...
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
Hiota0 = [H1,H2,H3,H4];
Tbar = 481;

load('output/jacobian_time1_1.mat');
H1_time1 = H;
clearvars H;
load('output/jacobian_time1_2.mat');
H2_time1 = H;
clearvars H;
load('output/jacobian_time1_3.mat');
H3_time1 = H;
clearvars H;
load('output/jacobian_time1_4.mat');
H4_time1 = H;
clearvars H;
load('output/GR_time1.mat');
p_ePath1 = p_ePath(1);

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
    'iota',0.94,'znext_gT1',0,'rhoZnext_g',nan,'znext_gT2',nan,'fixt',nan,...
    'forceEnd',144000);
    
% correct for iota
Hss = Hiota0;
Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),:) = ...
    (1-0.94)*Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),:);
Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),4*Tbar-2:(4+eStates)*(Tbar-1)+1) = ...
    Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),4*Tbar-2:(4+eStates)*(Tbar-1)+1) - ... 
    0.94*eye(eStates*(Tbar-1),eStates*(Tbar-1));
selectCols = ones(1,(4+eStates)*(Tbar-1)+1);
selectRows = ones(1,(7+eStates)*(Tbar-1)+1);
for i=1:(4+eStates)
    selectCols(1+i*(Tbar-1)) = 0;
end
for i=1:(7+eStates)
    selectRows(1+i*(Tbar-1)) = 0;
end
Hss1 = Hss(:,selectCols==1);
Hss1 = Hss1(selectRows==1,:);
H_time1 = [H1_time1,H2_time1,H3_time1,H4_time1];
H_time1(5*(Tbar-1)-4+1:5*(Tbar-1)-4+eStates*(Tbar-2),:) = ...
    ((1-0.94)/0.025)*H_time1(5*(Tbar-1)-4+1:5*(Tbar-1)-4+eStates*(Tbar-2),:);
H_time1(5*(Tbar-1)-4+1:5*(Tbar-1)-4+eStates*(Tbar-2),4*(Tbar-1)-2:(4+eStates)*(Tbar-2)+1) = ...
    H_time1(5*(Tbar-1)-4+1:5*(Tbar-1)-4+eStates*(Tbar-2),4*(Tbar-1)-2:(4+eStates)*(Tbar-2)+1) - ... 
    ((0.94-0.975)/0.025)*eye(eStates*(Tbar-2),eStates*(Tbar-2));  
dH = (H_time1 - Hss1)/(p_e - p_ePath1);  

shocksStruct = load('output/GR_iota0p94_expiration.mat','shocks');
betaShocks = shocksStruct.shocks;

% set to values starting in period 56
selectCols = ones(1,(4+eStates)*(Tbar-2)+1);
selectRows = ones(1,(7+eStates)*(Tbar-2)+1);
for i=1:(4+eStates)
    selectCols(1+i*(Tbar-2)-56+2:1+i*(Tbar-2)) = 0;
end
for i=1:(7+eStates)
    selectRows(1+i*(Tbar-2)-56+2:1+i*(Tbar-2)) = 0;
end
Hss_56 = Hss1(:,selectCols==1);
Hss_56 = Hss_56(selectRows==1,:);
dH_56 = dH(:,selectCols==1);
dH_56 = dH_56(selectRows==1,:);     

% Generate duration of UI in Great Recession assuming last extension 
%  would last for 24 months but then expires in December 2013
[~,dbarPathMonths,dbarShockMonths] = genUIGreatRecession_noexpiration(Tbar,24);

try
    GR(Tbar,80,0,betaShocks(1),0.95,...
        zeros(1,80),0,dbarPathMonths,dbarShockMonths,Hss_56,dH_56,urdata,0,...
        transParamsSticky,...
        RCEName,'output/GR_iota0p94_noexpiration',invariantParams,varyingParams,...
        equityShares_e,equityShares_u,epsSticky,epsMacro,znextDiff,goodsDiff,...
        tDiff,mDiff,equityValDiff,dvacpostDiff,(1-0.94)*dwageDiff,cTransSS,...
        56,'output/GR_iota0p94_ui_time56',betaShocks);             
catch err
    disp(['Error: ',err.identifier]);
end

disp(['Code completed in ',num2str(toc(masterTime)),' s.']);