%% simTransition_liquid_noshocks.m
% Simulates transitional dynamics starting from liquid wealth RCE
%  without any shocks to net out computational error from dynamics with shocks

masterTime = tic;

RCEName = 'output/RCE_liquid';
transName = 'output/trans_liquid_noshocks';
load(strcat(RCEName,'.mat'));

invariantParams = ...
    struct('numBetas',numBetas,'xi',xi,'sigma',sigma,...
    'numA',numA,'numP',numP,'rhoP',rhoP,'sigmaP',sigmaP,'numT',numT,'sigmaT',sigmaT,...
    'numUI',numUI,'Tua',Tua,'dbarbar',dbarbar,...
    'k',k,'eta',eta,'pEpsilon',pEpsilon,'phi',phi,'tauR',tauR);

Tbar = 481;

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
    initializeEquity(0.75,...
    (1/m)*qStat,z,[p_e_tilde*phi_e_tilde,(1-p_e_tilde)*phi_u_tilde]);
equityShares_e = repmat(equityShares,1,eStates);
equityShares_u = repmat(equityShares,1,uStates);

epsSticky = [0.1 0.00005 0.00005 0.0001 0.00005 0.00005];
transParamsSticky = ...
    struct('psi',360,'adjResource',0,'epsilon',6,'phiInt',1/m-1,'phiPi',1.5,...
    'phiC',0,'phiRho',0,'zlb',1,'fixedM',0,'fixedm',0,...
    'iota',0,'znext_gT1',0,'rhoZnext_g',nan,'znext_gT2',nan,'fixt',nan,...
    'forceEnd',72000);
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = Transition(transParamsSticky,...
    RCEName,transName,invariantParams,varyingParams,...
    equityShares_e,equityShares_u,epsSticky,0,zeros(1,Tbar),zeros(1,Tbar),...
    zeros(1,Tbar),zeros(1,Tbar),0,zeros(1,Tbar),zeros(eStates,Tbar),ones(1,Tbar),[],[],[],0,0,0);

disp(['Code completed in ',num2str(toc(masterTime)),' s.']);