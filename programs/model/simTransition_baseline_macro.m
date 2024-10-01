%% simTransition_baseline_macro.m
% Simulates various macro shocks under iota = 0.94 and sticky prices

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

epsSticky = [0.000001 0.000001 0.000001 0.00001 0.000001 0.000001];

Hss = H;
% adjust for real wage rigidity
Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),:) = ...
    (1-0.94)*Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),:);
Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),4*Tbar-2:(4+eStates)*(Tbar-1)+1) = ...
    Hss(5*Tbar-4+1:5*Tbar-4+eStates*(Tbar-1),4*Tbar-2:(4+eStates)*(Tbar-1)+1) - ... 
    0.94*eye(eStates*(Tbar-1),eStates*(Tbar-1));    

transParams = struct('psi',360,'adjResource',0,'epsilon',6,'phiInt',1/m-1,'phiPi',1.5,'phiC',1/12,'phiRho',0,'zlb',1,'fixedM',0,'fixedm',0,...
    'iota',0.94,'znext_gT1',0,'rhoZnext_g',nan,'znext_gT2',nan,'fixt',nan,...
    'forceEnd',72000); 
    
varyingParamsCell = cell(1,5);
suffixCell = {'beta','delta','zbar','abar','mbar'};

for i=1:5
    thisVaryingParams = varyingParams;    
    switch i
        case 1 % beta shock
            betaAR1 = zeros(1,Tbar);
            betaAR1(1) = 0.0005;
            for time2=2:180
                betaAR1(time2) = 0.95*betaAR1(time2-1);
            end          
            thisBetasTrans = betasTrans + repmat(betaAR1,numBetas,1);
            thisVaryingParams.betas = thisBetasTrans;
        case 2 % separation rate shock
            deltaAR1 = zeros(1,Tbar);
            deltaAR1(1) = 0.001;
            for time2=2:180
                deltaAR1(time2) = 0.95*deltaAR1(time2-1);
            end
            thisDeltaBarTrans = deltaBarTrans + deltaAR1;
            thisVaryingParams.deltaBar = thisDeltaBarTrans;
        case 3 % zbar shock
            zbarAR1 = zeros(1,Tbar);
            zbarAR1(1) = 0.01;
            for time2=2:180
                zbarAR1(time2) = 0.95*zbarAR1(time2-1);
            end
            thisZbarTrans = zbarTrans + zbarAR1;
            thisVaryingParams.zbar = thisZbarTrans;
        case 4 % productivity shock
            abarAR1 = zeros(1,Tbar);
            abarAR1(1) = 0.0005;            
            for time2=2:180
                abarAR1(time2) = 0.95*abarAR1(time2-1);
            end      
            thisAbarTrans = abarTrans.*(1+abarAR1);
            thisATrans = exp(log(aTrans) + repmat(abarAR1,size(aTrans,1),1));
            thisAPTrans = exp(log(aPTrans) + repmat(abarAR1,size(aPTrans,1),1));
            thisVaryingParams.abar = thisAbarTrans;
            thisVaryingParams.a = thisATrans;
            thisVaryingParams.aP = thisAPTrans;
        case 5 % mbar shock
            mbarAR1 = zeros(1,Tbar);
            mbarAR1(1) = -0.005;
            for time2=2:180
                 mbarAR1(time2) = 0.95*mbarAR1(time2-1);
            end          
            thisMbarTrans = mbarTrans + mbarAR1;
            thisVaryingParams.mbar = thisMbarTrans;
    end

    varyingParamsCell{i} = thisVaryingParams;
end

% Solve for transitional dynamics
parfor i=1:5
    warning('off','all'); % disable warnings when running on cluster
    
    try
        [~,~,~,~,~,~,~,...
            ~,~,~,~,~,~,~,~,~,...
            ~,~,~,~,~,~,~,~,~,~,~] = Transition(transParams,...
            RCEName,strcat('output/trans_baseline_iota0p94_',suffixCell{i}),invariantParams,varyingParamsCell{i},...
            equityShares_e,equityShares_u,epsSticky,i,znextDiff,goodsDiff,...
            tDiff,mDiff,equityValDiff,dvacpostDiff,(1-0.94)*dwageDiff,cTransSS,[],Hss,[],0,0,0);    			
    catch e
        disp(['Error on worker',num2str(i),': ',e.message]);
    end
end

disp(['Parallel code completed in ',num2str(toc(masterTime)),' s.']);

%% Shut down pool of workers
 delete(poolobj)