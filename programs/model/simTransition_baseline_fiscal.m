%% simTransition_baseline_fiscal.m
% Simulates government spending shocks under alternative assumptions on 
%  monetary policy and financing, in all cases with iota = 0.94 as in
%  calibration of Great Recession

%% Preliminaries 

% Get the number of cores allocated to your job
num_cores = str2num(getenv('SLURM_CPUS_PER_TASK'));

% Start a pool of workers
c = parcluster;
poolobj = parpool(c, num_cores);

%% Now my code...

masterTime = tic;

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

load('output/RCE_baseline');
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
g = [(prodStat/100)*ones(1,12),zeros(1,Tbar-12)];

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

% Adjust H for real wage rigidity
Hss = H;
Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),:) = ...
    (1-0.94)*Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),:);
Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),4*Tbar-2:(4+eStates)*(Tbar-1)+1) = ...
    Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),4*Tbar-2:(4+eStates)*(Tbar-1)+1) - ... 
    0.94*eye(eStates*(Tbar-1),eStates*(Tbar-1));

% Define cells of H and transParams: 
epsSticky = [0.000001 0.000001 0.000001 0.00001 0.000001 0.000001];
HssCell = cell(1,4);
transParamsCell = cell(1,4);
suffixCell = {'','_fixedi','_debt','_debt_fixedi'};

for i=1:4    
    thisHss = Hss;
    
    switch i
        case 1
            % 1: budget-balanced
            HssCell{i} = Hss;
            transParamsCell{i} = struct('psi',360,'adjResource',0,'epsilon',6,'phiInt',1/m-1,'phiPi',1.5,'phiC',1/12,'phiRho',0,'zlb',1,'fixedM',0,'fixedm',0,...
                'iota',0.94,'znext_gT1',0,'rhoZnext_g',nan,'znext_gT2',nan,'fixt',nan,...
                'forceEnd',72000); 
        case 2
            % 2: budget-balanced, fixed M for 18 months
            % adjust for fixed M in first 18 periods
            thisHss = thisHss(1:(end-2*(Tbar-1)),:) - ...            
               [zeros(4*Tbar-3,size(thisHss,2));...
               thisHss(end-2*(Tbar-1)+1:end-2*(Tbar-1)+18,:);...
               zeros(eStates*(Tbar-1)+((Tbar-1)-18),size(thisHss,2))]; 
            HssCell{i} = thisHss;
            transParamsCell{i} = struct('psi',360,'adjResource',0,'epsilon',6,'phiInt',1/m-1,'phiPi',1.5,'phiC',1/12,'phiRho',0,'zlb',1,'fixedM',18,'fixedm',0,...
                'iota',0.94,'znext_gT1',0,'rhoZnext_g',nan,'znext_gT2',nan,'fixt',nan,...
                'forceEnd',72000); 
        case 3
            % 3: deficit-financed
            HssCell{i} = Hss;
            transParamsCell{i} = struct('psi',360,'adjResource',0,'epsilon',6,'phiInt',1/m-1,'phiPi',1.5,'phiC',1/12,'phiRho',0,'zlb',1,'fixedM',0,'fixedm',0,...
                'iota',0.94,'znext_gT1',24,'rhoZnext_g',0.95,'znext_gT2',180,'fixt',t*ones(1,Tbar),...
                'forceEnd',72000); 
        case 4
            % 4: deficit-financed, fixed M for 18 months
            % adjust for fixed M in first 18 periods
            thisHss = thisHss(1:(end-2*(Tbar-1)),:) - ...            
               [zeros(4*Tbar-3,size(thisHss,2));...
               thisHss(end-2*(Tbar-1)+1:end-2*(Tbar-1)+18,:);...
               zeros(eStates*(Tbar-1)+((Tbar-1)-18),size(thisHss,2))];            
            HssCell{i} = thisHss;
            transParamsCell{i} = struct('psi',360,'adjResource',0,'epsilon',6,'phiInt',1/m-1,'phiPi',1.5,'phiC',1/12,'phiRho',0,'zlb',1,'fixedM',18,'fixedm',0,...
                'iota',0.94,'znext_gT1',24,'rhoZnext_g',0.95,'znext_gT2',180,'fixt',t*ones(1,Tbar),...
                'forceEnd',72000); 
    end
end

% Solve for transitional dynamics
parfor i=1:4
    warning('off','all'); % disable warnings when running on cluster
    
    try
        [~,~,~,~,~,~,~,...
            ~,~,~,~,~,~,~,~,~,...
            ~,~,~,~,~,~,~,~,~,~,~] = Transition(transParamsCell{i},...
            RCEName,strcat('output/trans_baseline_fiscal',suffixCell{i}),invariantParams,varyingParams,...
            equityShares_e,equityShares_u,epsSticky,i,znextDiff,goodsDiff,...
            tDiff,mDiff,equityValDiff,dvacpostDiff,(1-0.94)*dwageDiff,cTransSS,[],HssCell{i},[],0,0,0);       
    catch e
        disp(['Error on worker',num2str(i),': ',e.message]);
    end
end

disp(['Parallel code completed in ',num2str(toc(masterTime)),' s.']);

%% Shut down pool of workers
 delete(poolobj)