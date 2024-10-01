%% simJacobian_baseline.m
% Characterizes Jacobian starting from the baseline stationary RCE
%  columns (inputs): equityVal
%                    m(1:Tbar-1)
%                    theta(1:Tbar-1)
%                    t(1:Tbar-1)
%                    mu(1:Tbar-1)
%                    w(1,1:Tbar-1)
%                    ...
%                    w(eStates,1:Tbar-1)
%  rows (outputs): equityValTrans - equityVal - equityValDiffSS
%                  znextTrans(1:Tbar-1) - zsupplyTrans(1:Tbar-1) - znextDiffSS(1:Tbar-1)
%                  cTrans(1:Tbar-1) - prodTrans(1:Tbar-1) - goodsDiffSS(1:Tbar-1)
%                  vacpostTrans(1:Tbar-1) - sfNewHiresTrans(1:Tbar-1) - dvacpostDiffSS(1:Tbar-1)
%                  tTrans(1:Tbar-1)-t(1:Tbar-1)-tDiffSS(1:Tbar-1)
%                  mTrans(1:Tbar-1)-m(1:Tbar-1)-mDiffSS(1:Tbar-1)
%                  wTrans(1,1:Tbar-1)-w(1,1:Tbar-1)-dwageDiffSS(1,1:Tbar-1)
%                  ...
%                  wTrans(eStates,1:Tbar-1)-w(eStates,1:Tbar-1)-dwageDiffSS(eStates,1:Tbar-1)
%                  MTrans(1:Tbar-1) [used to consider fixed M subperiods]
% 				   cTrans(1:Tbar-1) [used to assess sensitivity to Taylor rule responding to output gap]
% Run in 4 components, so that convergence can be achieved in ~4 hours on
%  cluster.  Only need to copy this code 4 times, changing `varGroup':
% varGroup = 1: fills out columns wrt non wage variables
% varGroup = 2: fills out columns wrt wage variables 1,...,(eStates/3)
% varGroup = 3: fills out columns wrt wage variables (eStates/3)+1,...,2*(eStates/3)
% varGroup = 4: fills out columns wrt wage variables 2*(eStates/3)+1,...,3*(eStates/3)

%% Preliminaries 

% Get the number of cores allocated to your job
num_cores = str2num(getenv('SLURM_CPUS_PER_TASK'));

% Start a pool of workers
c = parcluster;
poolobj = parpool(c, num_cores);

%% Now my code...

masterTime = tic;
varGroup = 2;
if varGroup==1
    numVars = 4;
else
    numVars = 9;
end
Tbar = 481;
numWorkers = poolobj.NumWorkers;
numWorkersPerVar = numWorkers/numVars;
numPeriodsPerWorker = ceil((Tbar-1)/(numWorkers/numVars));

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

clearvars znextTrans zsupplyTrans prodTrans tTrans t ...
    mTrans m equityValTrans equityVal vacpostTrans sfNewHiresTrans wTrans w;

load('output/RCE_baseline');

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

% Initialize epsilon to be extremely loose, so that first derivatives can
%  be computed
epsSticky = [1 1 1 1 1 1];
transParamsSticky = ...
    struct('psi',360,'adjResource',0,'epsilon',6,'phiInt',1/m-1,'phiPi',1.5,'phiC',1/12,'phiRho',0,'zlb',1,'fixedM',0,'fixedm',0,...
    'iota',0,'znext_gT1',0,'rhoZnext_g',nan,'znext_gT2',nan,'fixt',nan,...
    'forceEnd',72000);

% Initialize Jacobian to be computed (see definitions)
if varGroup == 1
    H = nan(1+(7+eStates)*(Tbar-1),1+numVars*(Tbar-1));
else
    H = nan(1+(7+eStates)*(Tbar-1),numVars*(Tbar-1));
end
starting = struct('m',m*ones(1,Tbar),...
    'theta',theta*ones(1,Tbar),...
    't',t*ones(1,Tbar),...
    'equityVal',(1/m)*qStat,...
    'w',repmat(w',1,Tbar),...
    'mu',ones(1,Tbar)*(epsilon/(epsilon-1))*(1+tauR));

warning('off','all'); % disable warnings when running on cluster

% If varGroup = 1, compute first derivatives with respect to equityVal
if varGroup == 1
    thisStarting = starting;
    thisStarting = setfield(thisStarting,'equityVal',starting.equityVal+0.0001);
    [thisM,thisT,~,thisEquityVal,thisW,~,~,...
        thisZnextTrans,thisZsupplyTrans,thisCTrans,thisProdTrans,thisTTrans,...
        thismTrans,thisEquityValTrans,thisVacpostTrans,thisSfNewHiresTrans,thisWTrans,~,~,...
        ~,~,thisMTrans,~,~,~,~,~] = Transition(transParamsSticky,...
        RCEName,[],invariantParams,varyingParams,...
        equityShares_e,equityShares_u,epsSticky,0,znextDiff,goodsDiff,...
        tDiff,mDiff,equityValDiff,dvacpostDiff,dwageDiff,cTrans,thisStarting,[],[],0,0,0); 
    H(1,1) = (thisEquityValTrans-thisEquityVal-equityValDiff)/0.0001;
    H(2:Tbar,1) = (thisZnextTrans(1:Tbar-1)-thisZsupplyTrans(1:Tbar-1)-znextDiff(1:Tbar-1))'/0.0001;
    H(Tbar+1:2*Tbar-1,1) = (thisCTrans(1:Tbar-1)-thisProdTrans(1:Tbar-1)-goodsDiff(1:Tbar-1))'/0.0001;
    H(2*Tbar:3*Tbar-2,1) = (thisVacpostTrans(1:Tbar-1)-thisSfNewHiresTrans(1:Tbar-1)-dvacpostDiff(1:Tbar-1))'/0.0001;
    H(3*Tbar-1:4*Tbar-3,1) = (thisTTrans(1:Tbar-1)-thisT(1:Tbar-1)-tDiff(1:Tbar-1))'/0.0001;
    H(4*Tbar-2:5*Tbar-4,1) = (thismTrans(1:Tbar-1)-thisM(1:Tbar-1)-mDiff(1:Tbar-1))'/0.0001;
    for eInd=1:eStates
        H(5*Tbar-4+(eInd-1)*(Tbar-1)+1:5*Tbar-4+eInd*(Tbar-1),1) = (thisWTrans(eInd,1:Tbar-1)-thisW(eInd,1:Tbar-1)-dwageDiff(eInd,1:Tbar-1))'/0.0001;
    end
    H(5*Tbar-4+eStates*(Tbar-1)+1:6*Tbar-5+eStates*(Tbar-1),1) = (thisMTrans(1:Tbar-1)-m)'/0.0001;
	H(6*Tbar-5+eStates*(Tbar-1)+1:7*Tbar-6+eStates*(Tbar-1),1) = (thisCTrans(1:Tbar-1)-cTrans(1:Tbar-1))'/0.0001;
end

% Loop over other inputs and compute first derivatives
% Workers ordered by input and then the segment of time they hit
% If varGroup = 1, define
%  Hpar = reshape(H(:,2:end),(7*Tbar-6+eStates*(Tbar-1))*numPeriodsPerWorker,numWorkers)
% If varGroup > 1, define
%  Hpar = reshape(H,(7*Tbar-6+eStates*(Tbar-1))*numPeriodsPerWorker,numWorkers)
Hpar = nan((7*Tbar-6+eStates*(Tbar-1))*numPeriodsPerWorker,numWorkers);
parfor i=1:poolobj.NumWorkers 
    warning('off','all'); % disable warnings when running on cluster
    
    timeStart = mod(i-1,ceil((Tbar-1)/numPeriodsPerWorker))*numPeriodsPerWorker+1;
    timeEnd = min((mod(i-1,ceil((Tbar-1)/numPeriodsPerWorker))+1)*numPeriodsPerWorker,Tbar-1);
    if varGroup == 1
        input = ceil(i/numWorkersPerVar);
    else
        input = (varGroup-2)*numVars + ceil(i/numWorkersPerVar);
    end
    thisHpar = nan((7*Tbar-6+eStates*(Tbar-1))*numPeriodsPerWorker,1);
    for time=timeStart:timeEnd
        thisStarting = starting;
        thisVaryingParams = varyingParams;
        if varGroup == 1
            switch input
                case 1
                    newM = starting.m;
                    newM(time) = newM(time)+0.0001;
                    thisStarting = setfield(thisStarting,'m',newM);
                case 2
                    newTheta = starting.theta;
                    newTheta(time) = newTheta(time)+0.0001;
                    thisStarting = setfield(thisStarting,'theta',newTheta);
                case 3
                    newT = starting.t;
                    newT(time) = newT(time)+0.0001;
                    thisStarting = setfield(thisStarting,'t',newT);  
                case 4
                    newMu = starting.mu;
                    newMu(time) = newMu(time)+0.0001;
                    thisStarting = setfield(thisStarting,'mu',newMu);                              
            end
        else
            newW = starting.w;
            newW(input,time) = newW(input,time)+0.0001;
            thisStarting = setfield(thisStarting,'w',newW);
        end
        
        [thisM,thisT,~,thisEquityVal,thisW,~,~,...
            thisZnextTrans,thisZsupplyTrans,thisCTrans,thisProdTrans,thisTTrans,...
            thismTrans,thisEquityValTrans,thisVacpostTrans,thisSfNewHiresTrans,thisWTrans,~,~,...
            ~,~,thisMTrans,~,~,~,~,~] = Transition(transParamsSticky,...
            RCEName,[],invariantParams,thisVaryingParams,...
            equityShares_e,equityShares_u,epsSticky,i,znextDiff,goodsDiff,...
            tDiff,mDiff,equityValDiff,dvacpostDiff,dwageDiff,cTrans,thisStarting,[],[],0,0,0);        
        thisHpar((time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+1,1) = ...
            (thisEquityValTrans-thisEquityVal-equityValDiff)/0.0001;
        thisHpar((time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+2:(time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+Tbar,1) = ...
            (thisZnextTrans(1:Tbar-1)-thisZsupplyTrans(1:Tbar-1)-znextDiff(1:Tbar-1))'/0.0001;
        thisHpar((time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+Tbar+1:(time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+2*Tbar-1,1) = ...
            (thisCTrans(1:Tbar-1)-thisProdTrans(1:Tbar-1)-goodsDiff(1:Tbar-1))'/0.0001;
        thisHpar((time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+2*Tbar:(time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+3*Tbar-2,1) = ...
            (thisVacpostTrans(1:Tbar-1)-thisSfNewHiresTrans(1:Tbar-1)-dvacpostDiff(1:Tbar-1))'/0.0001;
        thisHpar((time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+3*Tbar-1:(time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+4*Tbar-3,1) = ...
            (thisTTrans(1:Tbar-1)-thisT(1:Tbar-1)-tDiff(1:Tbar-1))'/0.0001;
        thisHpar((time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+4*Tbar-2:(time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+5*Tbar-4,1) = ...
            (thismTrans(1:Tbar-1)-thisM(1:Tbar-1)-mDiff(1:Tbar-1))'/0.0001;
        eInd=1;
        while eInd <= eStates 
            thisHpar((time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+5*Tbar-4+(eInd-1)*(Tbar-1)+1:(time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+5*Tbar-4+eInd*(Tbar-1),1) = (thisWTrans(eInd,1:Tbar-1)-thisW(eInd,1:Tbar-1)-dwageDiff(eInd,1:Tbar-1))'/0.0001;
            eInd = eInd + 1;
        end
        thisHpar((time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+5*Tbar-3+eStates*(Tbar-1):(time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+6*Tbar-5+eStates*(Tbar-1),1) = ...
            (thisMTrans(1:Tbar-1)-m)'/0.0001;
        thisHpar((time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+6*Tbar-4+eStates*(Tbar-1):(time-timeStart)*(7*Tbar-6+eStates*(Tbar-1))+7*Tbar-6+eStates*(Tbar-1),1) = ...
            (thisCTrans(1:Tbar-1)-cTrans(1:Tbar-1))'/0.0001;			
    end
    Hpar(:,i) = thisHpar;
end
Hpar_reshape = reshape(Hpar,(7*Tbar-6+eStates*(Tbar-1)),numPeriodsPerWorker*numWorkers);
% populate H
colsPerInput = numPeriodsPerWorker*numWorkersPerVar;
colsToSelect = nan(1,numVars*(Tbar-1));
for varInd=1:numVars
    colsToSelect((varInd-1)*(Tbar-1)+1:varInd*(Tbar-1)) = ...
        [(varInd-1)*colsPerInput+1:(varInd-1)*colsPerInput+Tbar-1];
end
if varGroup == 1
    H(:,2:end) = Hpar_reshape(:,colsToSelect);
else
    H = Hpar_reshape(:,colsToSelect);
end

save(['output/jacobian_baseline_',num2str(varGroup)],'H');
disp(['Parallel code completed in ',num2str(toc(masterTime)),' s.']);

%% Shut down pool of workers
delete(poolobj)