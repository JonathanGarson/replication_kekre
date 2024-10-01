function [m,t,theta,equityVal,w,mu,znext_g,...
    znextTrans,zsupplyTrans,cTrans,prodTrans,tTrans,...
    mTrans,equityValTrans,vacpostTrans,sfNewHiresTrans,wTrans,swTrans,sfTrans,...
    yTrans,PiTrans,MTrans,p_e,p_u,sbarTrans,avgW,dist2] = ...
    Transition(transParams,...
    RCEName,toName,invariantParams,varyingParams,equityShares_eSS,equityShares_uSS,...
    eps,numTabs,znextDiffSS,goodsDiffSS,tDiffSS,...
    mDiffSS,equityValDiffSS,dvacpostDiffSS,dwageDiffSS,cTransSS,starting,H,dist1,dampen,pe,saveWorkspace)
%% Transition.m
% Solves for transitional dynamics
%   transParams = parameters specific to transitional dynamics
%     psi = Rotemberg adj costs in [0,inf]
%     adjResource = 1 if adjustment costs are resource costs, 0 otherwise
%     epsilon = elasticity of substitution across varieties
%     phiInt = Taylor rule intercept 
%     phiPi = Taylor rule coefficient on inflation
%     phiC = Taylor rule coefficient on output
%     phiRho = Taylor rule smoothing parameter
%     zlb = indicator for whether Taylor Rule should be subject to ZLB (1=yes)
%     fixedM = number of periods that nominal interest rate is fixed
%     fixedm = number of periods that real interest rate is fixed 
%     iota = weight on steady-state real wage schedule
%     znext_gT1 = period of simulation through which debt financing occurs 
%      (if = 0, then taxes balance budget each period and znext_gT2 
%      and fixt below are irrelevant)
%     znext_gT2 = period of simulation at which debt is assumed to 
%      return to SS level
%     fixt = fixed path of t used for first znext_gT1 periods
%   RCEName = filename from which stationary RCE should be called
%   toName = (optional) filename to which transitional dynamics should be saved
%   invariantParams = parameters assumed to be invariant over transition
%   varyingParams = parameters which may vary over transition
%   equityShares_{i}SS = assumed fraction of portfolio held in firm equity 
%    in the stationary RCE for i \in {e,u} agents.  A gridpointsStat x iStates matrix.
%   eps = convergence criteria, where multiple rows imply dynamics will be
%    saved at multiple points.  Each row i interpreted as:
%     [eps_znext(i),eps_t(i),eps_q0(i),eps_vac(i),eps_wage(i)] in flex price case
%     [eps_goods(i),eps_t(i),eps_m(i),eps_q0(i),eps_vac(i),eps_wage(i)] in sticky price case
%      where eps_m(i) only relevant if transType.psi < inf 
%   numTabs = number of tabs in output
%   znextDiffSS = znextTrans - zsupplyTrans used as benchmark for converg
%   goodsDiffSS = cTrans - prodTrans used as benchmark for converg
%   tDiffSS = tTrans - t used as benchmark for converg
%   mDiffSS = mTrans - m used as benchmark for converg
%   equityValDiffSS = equityValTrans - equityVal used as benchmark for converg
%   dvacpostDiffSS = vacpostTrans - sfNewHiresTrans used as a benchmark for converg
%   dwageDiffSS = wTrans - w used as a benchmark for converg
%   cTransSS = cTrans used as a benchmark for cTrans in Taylor rule 
%    (only relevant if phiC > 0)
%   starting = (optional) starting macro aggregates
%   H = (optional) jacobian to use in quasi-Newton algorithm: 
%    columns: equityVal
%             m(1:Tbar-1)
%             theta(1:Tbar-1)
%             t(1:Tbar-1)
%             mu(1:Tbar-1)
%             w(1,1:Tbar-1)
%             ...
%             w(eStates,1:Tbar-1)
%    rows: equityValTrans - equityVal - equityValDiffSS
%          znextTrans(1:Tbar-1) - zsupplyTrans(1:Tbar-1) - znextDiffSS(1:Tbar-1)
%          cTrans(1:Tbar-1) - prodTrans(1:Tbar-1) - goodsDiffSS(1:Tbar-1)
%           [use one of the above two but not both]
%          vacpostTrans(1:Tbar-1) - sfNewHiresTrans(1:Tbar-1) - dvacpostDiffSS(1:Tbar-1)
%          tTrans(1:Tbar-1)-t(1:Tbar-1)-tDiffSS(1:Tbar-1)
%          mTrans(1:Tbar-1)-m(1:Tbar-1)-mDiffSS(1:Tbar-1)
%          wTrans(1,1:Tbar-1)-w(1,1:Tbar-1)-dwageDiffSS(1,1:Tbar-1)
%          ...
%          wTrans(eStates,1:Tbar-1)-w(eStates,1:Tbar-1)-dwageDiffSS(eStates,1:Tbar-1)
%          MTrans(1:Tbar-1)
%   dist1 = (optional) struct holding initial distribution of agents along state space
%   dampen = parameter used to dampen updating of (m,theta,...) 
%   pe = 1 if PE dynamics are being characterized, 0 otherwise
%    note pe = 1 is assumed to be accompanied by psi > 0 (sticky prices) and H provided
%   saveWorkspace = (optional) indicator for saving entire workspace 

overallRuntime = tic;
tabs = '.';
for i=1:numTabs
    tabs = sprintf('   %s',tabs);
end

%% Set parameters describing transition

global psi adjResource epsilon phiInt phiPi phiC phiRho zlb fixedM fixedm iota znext_gT1 rhoZnext_g znext_gT2 fixt forceEnd;

psi = transParams.psi;
if psi == 0
    prices = 'flex';
    saveFilename = strcat(toName,'_flex');
else
    prices = 'sticky';
    saveFilename = strcat(toName,'_sticky');
end
adjResource = transParams.adjResource;

epsilon = transParams.epsilon;
phiInt = transParams.phiInt;
phiPi = transParams.phiPi;
phiC = transParams.phiC;
phiRho = transParams.phiRho;
zlb = transParams.zlb;
fixedM = transParams.fixedM;
fixedm = transParams.fixedm;
iota = transParams.iota;
znext_gT1 = transParams.znext_gT1;
rhoZnext_g = transParams.rhoZnext_g;
znext_gT2 = transParams.znext_gT2;
fixt = transParams.fixt;
forceEnd = transParams.forceEnd;

%% Set economic parameters invariant over transition

global numBetas xi sigma ...
    numA numP rhoP sigmaP numT sigmaT ...
    numUI Tua dbarbar ...
    k eta ...
    pEpsilon phi ...
    tauR;

numBetas = invariantParams.numBetas;
xi = invariantParams.xi;
sigma = invariantParams.sigma;
numA = invariantParams.numA;
numP = invariantParams.numP;
rhoP = invariantParams.rhoP;
sigmaP = invariantParams.sigmaP;
numT = invariantParams.numT;
sigmaT = invariantParams.sigmaT;
numUI = invariantParams.numUI;
Tua = invariantParams.Tua;
dbarbar = invariantParams.dbarbar;
k = invariantParams.k;
eta = invariantParams.eta;
pEpsilon = invariantParams.pEpsilon;
phi = invariantParams.phi;
tauR = invariantParams.tauR;

%% Set economic parameters which may vary over transition

global Tbar betas ...
    Ta Taui...
    zeta dbar rr ...
    deltaBar eps_delta_a eps_delta_beta mbar lambda ...
    abar a aP...
    zbar znext_g ...
    omega0 omega1 omega2 uimax g;

Tbar = varyingParams.Tbar;
betas = varyingParams.betas(:,1:Tbar);
Ta = varyingParams.Ta(:,:,1:Tbar);
Taui = varyingParams.Taui(:,:,1:Tbar);
zeta = varyingParams.zeta(:,1:Tbar);
dbar = varyingParams.dbar(:,1:Tbar);
rr = varyingParams.rr(:,1:Tbar);
deltaBar = varyingParams.deltaBar(:,1:Tbar);
eps_delta_a = varyingParams.eps_delta_a(:,1:Tbar);
eps_delta_beta = varyingParams.eps_delta_beta(:,1:Tbar);
mbar = varyingParams.mbar(:,1:Tbar);
lambda = varyingParams.lambda(:,1:Tbar);
abar = varyingParams.abar(:,1:Tbar);
a = varyingParams.a(:,1:Tbar);
aP = varyingParams.aP(:,1:Tbar);
zbar = varyingParams.zbar(:,1:Tbar);
znext_g = varyingParams.znext_g(:,1:Tbar);
omega0 = varyingParams.omega0(:,1:Tbar);
omega1 = varyingParams.omega1(:,1:Tbar);
omega2 = varyingParams.omega2(:,1:Tbar);
uimax = varyingParams.uimax(:,1:Tbar);
g = varyingParams.g(:,1:Tbar);

%% Define initial and final steady-states

SS = load(strcat(RCEName,'.mat'));
mSS = SS.m;
tSS = SS.t;
thetaSS = SS.theta;
zbarSS = SS.zbar;
zSS = SS.zStat;
v_eSS = SS.v_eStat;
v_uSS = SS.v_uStat;
sSS = SS.sStat;
c_eSS = SS.c_eStat;
znext_eSS = SS.znext_eStat;
c_uSS = SS.c_uStat;
znext_uSS = SS.znext_uStat;
p_e_tildeSS = SS.p_e_tilde;
phi_e_tildeSS = SS.phi_e_tilde;
phi_u_tildeSS = SS.phi_u_tilde;
znext_gSS = SS.znext_g;
qSS = SS.qStat;
dbarSS = SS.dbar;
phi_eSS = SS.phi_e;
p_eSS = SS.p_e;
sbarSS = SS.sbar;
abarSS = SS.abar;
aSS = SS.a;
wSS = SS.w;
bSS = SS.b;
zetaSS = SS.zeta;
sfSS = SS.sf;
deltaBarSS = SS.deltaBar;
clearvars SS;

%% Set tuning parameters

% Define the number of employed and unemployed states
eStates = numA*numBetas;
uStates = (dbarbar+1)*numT*numUI*numBetas;

numEps = size(eps,1);
epsInd = 1;
if strcmp(prices,'flex')
    eps_znext = eps(epsInd,1);
    eps_t = eps(epsInd,2);
    eps_equityVal = eps(epsInd,3);
    eps_vac = eps(epsInd,4);
    eps_wage = eps(epsInd,5);
else
    eps_znext = eps(epsInd,1);
    eps_t = eps(epsInd,2);
    eps_m = eps(epsInd,3);
    eps_equityVal = eps(epsInd,4);
    eps_vac = eps(epsInd,5);
    eps_wage = eps(epsInd,6);
end

if isempty(H)
    delta_m = 0.002;
    delta_theta = 3;
    delta_w = 0.1;
    delta_mu = 0.1;
    delta_t = 1;
    delta_equityVal = 1; 
else 
    if strcmp(prices,'flex')
        if iota < 1
            Hinv = H([1:Tbar,2*Tbar:(4*Tbar-3),(5*Tbar-3):(5*Tbar-4+eStates*(Tbar-1))],...
                [1:(3*Tbar-2),(4*Tbar-2):(4*Tbar-3+eStates*(Tbar-1))])^(-1);
        else
            Hinv = H([1:Tbar,2*Tbar:(4*Tbar-3)],1:(3*Tbar-2))^(-1);            
        end
    else
        if iota < 1
            if ~pe
                Hinv = H([1:Tbar,2*Tbar:4*Tbar-3,4*Tbar-3+1+fixedm:(5*Tbar-4+eStates*(Tbar-1))],...
                    [1,2+fixedm:(4*Tbar-3+eStates*(Tbar-1))])^(-1);
            else
                Hinv = H([3*(Tbar-1)+2:3*(Tbar-1)+Tbar,5*(Tbar-1)+2:5*(Tbar-1)+1+eStates*(Tbar-1)],...
                    [2*(Tbar-1)+2:2*(Tbar-1)+Tbar,4*(Tbar-1)+2:4*(Tbar-1)+1+eStates*(Tbar-1)])^(-1);
            end
        else
            if ~pe
                Hinv = H([1:Tbar,2*Tbar:4*Tbar-3,4*Tbar-3+1+fixedm:(5*Tbar-4)],...
                    [1,2+fixedm:(4*Tbar-3)])^(-1);            
            else
                Hinv = H(3*(Tbar-1)+2:3*(Tbar-1)+Tbar,...
                    2*(Tbar-1)+2:2*(Tbar-1)+Tbar)^(-1);                
            end
        end            
    end
end

weights = exp(-0.015*(1:Tbar));

tol = 1e-12; % meant to weed out computational error

%% Other preliminaries before starting algorithm

if isempty(dist1)
    % Define the fixed z grid used throughout the algorithm (=zSS)
    z = zSS;
    gridPoints = length(z);
    p_e_tilde1 = p_e_tildeSS;
    phi_e_tilde1 = phi_e_tildeSS;
    phi_u_tilde1 = phi_u_tildeSS;
    lastPosDens = find(max(phi_e_tilde1,[],2) > 0,1,'last');
    
    equityShares_e = interp1(zSS,equityShares_eSS,z);
    equityShares_u = interp1(zSS,equityShares_uSS,z);
    
    zbar0 = zbarSS;
else
    z = dist1.z;
    gridPoints = length(z);
    p_e_tilde1 = dist1.p_e_tilde1;
    phi_e_tilde1 = dist1.phi_e_tilde1;
    phi_u_tilde1 = dist1.phi_u_tilde1;
    lastPosDens = find(max(phi_e_tilde1,[],2) > 0,1,'last');
     
    equityShares_e = equityShares_eSS;
    equityShares_u = equityShares_uSS;
    
    zbar0 = dist1.zbar0;
end

% Define zFine, which includes 30 gridpoints in between the last two 
%  gridpoints of z (used for backward policy iteration *only*)
zFine = [z(1:gridPoints-2);linspace(z(gridPoints-1),z(gridPoints),30)'];
gridPointsFine = length(zFine);
origGridInds = [1:gridPoints-1,gridPointsFine]';

% Define marginal utility of real income for employed
uPrime = @(c) c.^(-sigma);

% Define functions returning job-finding and vacancy-filling probabilities
lambdaAdj = nan(dbarbar+1,Tbar);
lambdaAdj(1:(dbarSS+1),:) = repmat(lambda,dbarSS+1,1).*repmat((0:dbarSS)',1,Tbar);
lambdaAdj((dbarSS+2):(dbarbar+1),:) = repmat(lambda*dbarSS,dbarbar-dbarSS,1);
p = @(t,time) repmat(repmat((mbar(time).^(1-eta)).*(t.^eta),(dbarbar+1),1).*...
    exp(lambdaAdj(:,time)),numT*numUI*numBetas,1); % will be of size [dbarbar+1]*numT*numUI*numBetas x length(time)
p0 = @(t,time) (mbar(time).^(1-eta)).*(t.^eta);
q = @(t,time) (mbar(time).^(1-eta)).*(t.^(eta-1));

% And now define matrices defining (in general) date-specific separation 
%  rates and {eState,uState}->{eState,uState} transition probabilities

% - first by computing number of beta and wage for each state when employed
eBetas = ceil((1:numA*numBetas)/(numA));
eW = repmat(1:numA,1,numBetas);
% - and then by defining transition functions for next state 
ee = @(thisState) cumsum([floor((thisState-1)/numA)*numA+1;...
   ones(numA-1,length(thisState))],1);
ee_prob = @(thisState,time) Ta(eW(thisState)',:,time)';

eu = @(thisState) cumsum([(eBetas(thisState)-1).*((dbarbar+1)*numT*numUI)+1;...
        ones(numT*numUI-1,length(thisState))*(dbarbar+1)],1);
eu_prob = @(thisState,time) kron(Taui(eW(thisState),:,time)',pEpsilon');
Tee = cell(1,Tbar);
Teu = cell(1,Tbar);
for time=1:Tbar
    Tee{time} = sparse(reshape(repmat(1:eStates,size(ee(1:eStates),1),1),eStates*size(ee(1:eStates),1),1),...
        reshape(ee(1:eStates),eStates*size(ee(1:eStates),1),1),...
        reshape(ee_prob(1:eStates,time),eStates*size(ee(1:eStates),1),1),...
        eStates,eStates);
    Teu{time} = sparse(reshape(repmat(1:eStates,size(eu(1:eStates),1),1),eStates*size(eu(1:eStates),1),1),...
        reshape(eu(1:eStates),eStates*size(eu(1:eStates),1),1),...
        reshape(eu_prob(1:eStates,time),eStates*size(eu(1:eStates),1),1),...
        eStates,uStates);    
end

% - analogously, by computing number of beta, UI schedule, 
%   endowment trans component (if applic), and duration for each state when unemployed
uBetas = ceil((1:(dbarbar+1)*numT*numUI*numBetas)/...
    ((dbarbar+1)*numT*numUI));
uUI = repmat(ceil((1:(dbarbar+1)*numT*numUI)/((dbarbar+1)*numT)),...
    1,numBetas);
uEndow = repmat(ceil((1:(dbarbar+1)*numT)/(dbarbar+1)),...
    1,numUI*numBetas);
uDur = repmat(1:(dbarbar+1),1,numT*numUI*numBetas);
% - and then by defining transition functions for next state
ue = @(thisState) cumsum([(uBetas(thisState)-1).*(numA) + 1;...
    ones(numA-1,length(thisState))],1);
ue_prob = @(thisState) Tua((uUI(thisState)-1)*numT+uEndow(thisState),:)';
uu = @(thisState) cumsum([(uBetas(thisState)-1)*((dbarbar+1)*numT*numUI) + ...
    (uUI(thisState)-1)*((dbarbar+1)*numT) + ...
    uDur(thisState) + (uDur(thisState)<(dbarbar+1));...
    (dbarbar+1)*ones(numT-1,length(thisState))],1);
uu_prob = @(thisState) repmat(pEpsilon',1,length(thisState));
Tue = sparse(reshape(repmat(1:uStates,size(ue(1:uStates),1),1),uStates*size(ue(1:uStates),1),1),...
    reshape(ue(1:uStates),uStates*size(ue(1:uStates),1),1),...
    reshape(ue_prob(1:uStates),uStates*size(ue(1:uStates),1),1),...
    uStates,eStates);
Tuu = sparse(reshape(repmat(1:uStates,size(uu(1:uStates),1),1),uStates*size(uu(1:uStates),1),1),...
    reshape(uu(1:uStates),uStates*size(uu(1:uStates),1),1),...
    reshape(uu_prob(1:uStates),uStates*size(uu(1:uStates),1),1),...
    uStates,uStates);
 % ({u}->{e,u} transition matrices are time-invariant)
 
% - and finally computing a matrix of separation rates at each (eState,time)
deltaMat = repmat(deltaBar,eStates,1).*...
    (repmat(ones(numA,Tbar) + ...
    repmat(eps_delta_a/deltaBarSS,numA,1).*(log(aP)-repmat(log(abar),numA,1)),numBetas,1) + ...
    repmat(eps_delta_beta/deltaBarSS,numA*numBetas,1).*...
    (betas(eBetas,:)-repmat(mean(betas,1),numA*numBetas,1)));

% Finally, initialize final period values and consumption policies
v_eTbar = interp1q(zSS,v_eSS,zFine);
c_eTbar = interp1q(zSS,c_eSS,zFine);
znext_eTbar = interp1q(zSS,znext_eSS,zFine);
v_uTbar = interp1q(zSS,v_uSS,zFine);
c_uTbar = interp1q(zSS,c_uSS,zFine);
znext_uTbar = interp1q(zSS,znext_uSS,zFine);
sTbar = interp1q(zSS,sSS,zFine);

%% Run shooting algorithm to characterize transition

% Keep track of convergence
converged = 0;
numIter = 0;
d = @(v1,v2) max(max(abs((v1 - v2)),[],1),[],2);

% Step 0) Guess initial paths for macro aggregates
if isempty(starting)
    m = ones(1,Tbar)*mSS;
    theta = ones(1,Tbar)*thetaSS;
    t = ones(1,Tbar)*tSS;
    equityVal = (1/mSS)*qSS;
    w = repmat(wSS',1,Tbar);
    mu = ones(1,Tbar)*(epsilon/(epsilon-1))*(1+tauR);
else
    if strcmp(prices,'flex')
        m = starting.m;
        theta = starting.theta;
        t = starting.t;
        equityVal = starting.equityVal;
        w = starting.w;
        mu = ones(1,Tbar)*(epsilon/(epsilon-1))*(1+tauR);
    else
        m = starting.m;
        theta = starting.theta;
        t = starting.t;
        equityVal = starting.equityVal;
        w = starting.w;
        mu = starting.mu;
    end
end
lastAtZLB = 0;

while converged == 0
    singleIter = tic;
    numIter = numIter + 1;
    
    disp([tabs,'Iteration ',num2str(numIter)]);
    
    % Solve for UI consistent with conjectured wages
    b = nan(uStates,Tbar);
    for time=1:(Tbar-1)
        wP = pEpsilon*reshape(w(:,time)',numT,numP*numBetas);
        [thisB,~] = genTransfers(numBetas,numP,numT,dbarbar,dbar(time),rr(time)*(1-omega0(time)),...
            w(:,time)',wP,omega1(time),omega2(time),zeta(time),uimax(time));  
        b(:,time) = thisB';
    end
    b(:,Tbar) = bSS';
        
    % Define income of employed and unemployed workers
    y_e = w - repmat(t,eStates,1); % a numA*numBetas x Tbar matrix
    y_u = b; % a [dbarbar+1]*numT*numUI*numBetas x Tbar matrix 

    % EGM-augmented VFI backwards
    [v_e,v_u,s,c_e,znext_e,c_u,znext_u,lastPolDens] = ...
        iteratePoliciesBackward(zFine,v_eTbar,c_eTbar,znext_eTbar,...
        v_uTbar,c_uTbar,znext_uTbar,sTbar,...
        m,p(theta,1:Tbar),deltaMat,Tee,Teu,Tue,Tuu,...
        y_e,y_u,[zbar0,zbar(1:(Tbar-1))],zbar,...
        betas(eBetas,:),betas(uBetas,:),origGridInds(lastPosDens),origGridInds,'on',tabs);
    display([tabs,'Worker policies solved backwards.']);
    
    % Iterate policies forward to characterize evolution of phi
    %  - first define (beginning of period) value of equity that would
    %  have prevailed prior to any shock
    if isempty(dist1)
        initEquityVal = (1/mSS)*qSS;
    else
        initEquityVal = dist1.initEquityVal;
    end
    %  - then account for equity revaluation among equityholders, given 
    %  value of equity that would have prevailed prior to any shock 
    phi_tilde1 = [p_e_tilde1*phi_e_tilde1,...
        (1-p_e_tilde1)*phi_u_tilde1];
    phi_1_reval = revalueEquity(...
        initEquityVal,equityVal,[equityShares_e,equityShares_u],...
        z,phi_tilde1,lastPosDens);     
    phi_e_tilde1_reval = phi_1_reval(:,1:eStates)/...
        sum(sum(phi_1_reval(:,1:eStates)));
    phi_u_tilde1_reval = phi_1_reval(:,(eStates+1):(eStates+uStates))/...
        sum(sum(phi_1_reval(:,(eStates+1):(eStates+uStates))));     
    %  - ignore the distribution following lastPol if aggregate wealth
    %     above it is below 1e-20; else, throw an error
    lastPolDens1_reval = sum(phi_1_reval(lastPolDens+1:gridPoints,:),2)'*...
        z(lastPolDens+1:gridPoints);
    if lastPolDens1_reval < 1e-20
        phi_1_reval(lastPolDens+1:gridPoints,:) = 0;
        phi_1_reval(1:lastPolDens,:) = phi_1_reval(1:lastPolDens,:)*(1/(1-lastPolDens1_reval));
        phi_e_tilde1_reval(lastPolDens+1:gridPoints,:) = 0;
        phi_e_tilde1_reval(1:lastPolDens,:) = phi_e_tilde1_reval(1:lastPolDens,:)*(1/(1-lastPolDens1_reval));        
        phi_u_tilde1_reval(lastPolDens+1:gridPoints,:) = 0;
        phi_u_tilde1_reval(1:lastPolDens,:) = phi_u_tilde1_reval(1:lastPolDens,:)*(1/(1-lastPolDens1_reval));
    else
        display([tabs,'Wealth held by the tail without policies > 1e-20.']);
    end
    %  - then iterate forward from _reval distributions      
    [p_e_tilde,phi_e_tilde,phi_u_tilde,p_e,phi_e,phi_u] = ...
       iterateWorkersForward(z,s,znext_e,znext_u,...
       p(theta,1:Tbar),deltaMat,Tee,Teu,Tue,Tuu,...
       p_e_tilde1,phi_e_tilde1_reval,phi_u_tilde1_reval,'on',tabs);
   
    display([tabs,'Worker distribution iterated forwards.']);    

    % Check that no agents are at artifical upper bound of z along
    % transition 
    densUpperBound = max(max(squeeze(phi_e(gridPoints,:,:))));
    if densUpperBound > 0
        display([tabs,'Some agents are at upper bound of z during transition!']);
        %pause;
    end
    
    % Iterate backwards to compute firm surplus 
    sfTrans = solveFirmBackwards(Tee,deltaMat,m,w,mu,sfSS);
    display([tabs,'Given m, w, and mu, firm surplus solved backwards.']);
    
    % Define unemployment rates by duration and average wage
    p_u = nan(dbarbar+1,Tbar);
    for dur=1:(dbarbar+1)
        p_u(dur,:) = (1-p_e).*sum(squeeze(sum(phi_u(:,uDur==dur,:),2)),1); 
    end   
    avgW = sum(w.*squeeze(sum(phi_e,1)),1);    
    
    % Compute macro aggregates, assess convergence and update 
        
    converged = 1; 

    % Lump-sum tax (and government debt, if znext_gT1 > 0)
    if znext_gT1 > 0
        znext_gLagged = nan(1,znext_gT1+1);
        znext_gLagged(1) = znext_gSS;
        for tInd = 1:znext_gT1
            znext_gLagged(tInd+1) = (1/m(tInd))*(znext_gLagged(tInd) + ...
                p_e(tInd)*fixt(tInd) - g(tInd) - ...
                (1-p_e(tInd)).*nansum(reshape(sum(phi_u(:,:,tInd),1),uStates,1).*b(:,tInd),1));          
        end
        znext_g(1:znext_gT1) = znext_gLagged(2:(znext_gT1+1));
        for tInd = (znext_gT1+1):(znext_gT2-1)
            znext_g(tInd) = znext_gSS + (znext_g(tInd-1)-znext_gSS)*rhoZnext_g;
        end
        znext_g(znext_gT2:Tbar) = znext_gSS;
    end
    tTrans = (1./p_e).*(m.*znext_g - [znext_gSS,znext_g(1:(Tbar-1))] + g + ...
        (1-p_e).*nansum(reshape(sum(phi_u,1),uStates,Tbar).*b,1));        
    tDist = d((tTrans-t).*weights,tDiffSS.*weights); 
    tDiff = (tTrans-t)-tDiffSS;    
    
    % Aggregate borrowing and supply of assets
    znextTrans = nansum(squeeze(nansum(znext_e.*phi_e.*repmat(reshape(p_e,1,1,Tbar),gridPoints,eStates,1),2)) + ...
        squeeze(nansum(znext_u.*phi_u.*repmat(reshape((1-p_e),1,1,Tbar),gridPoints,uStates,1),2)),1);
    sbarTrans = nansum(squeeze(nansum(...
        repmat(reshape(p(ones(1,Tbar),1:Tbar)./...
        repmat(mbar.^(1-eta),uStates,1),1,uStates,Tbar),gridPoints,1,1).*...
        s.*phi_u_tilde.*repmat(reshape(1-p_e_tilde,1,1,Tbar),gridPoints,uStates,1),2)),1);
    yTrans = sum(a(eW,:).*reshape(sum(phi_e,1),eStates,Tbar).*repmat(p_e,eStates,1),1) - ...
        abar.*k.*theta.*sbarTrans;
    if strcmp(prices,'sticky')
        PiTrans = solvePiBackwards(mu,m,yTrans);
        if psi < inf
            acTrans = (psi/2)*(PiTrans.^2).*yTrans;
        else
            acTrans = zeros(1,Tbar);
        end
    else
        PiTrans = zeros(1,Tbar);
        acTrans = zeros(1,Tbar);
    end
    qTrans = valueEquity(m,...
        [sum((repmat(a,numBetas,1)-w).*reshape(sum(phi_e,1),eStates,Tbar).*repmat(p_e,eStates,1),1) - ...
        abar.*k.*theta.*sbarTrans - ...
        adjResource*acTrans,((1-mSS)/mSS)*qSS],qSS);
    zsupplyTrans = -znext_g + (1./m).*qTrans;
    znextDist = d((znextTrans-zsupplyTrans).*weights,znextDiffSS.*weights);
    znextDiff = (znextTrans-zsupplyTrans)-znextDiffSS;    

    % Aggregate consumption and production
    cTrans = nansum(squeeze(nansum(c_e.*phi_e.*repmat(reshape(p_e,1,1,Tbar),gridPoints,eStates,1),2)) + ...
        squeeze(nansum(c_u.*phi_u.*repmat(reshape(1-p_e,1,1,Tbar),gridPoints,uStates,1),2)),1);
    prodTrans = yTrans - adjResource*acTrans;
    goodsDist = d((cTrans+g-prodTrans).*weights,goodsDiffSS.*weights); 
    goodsDiff = (cTrans+g-prodTrans)-goodsDiffSS; 

    % Revalued equity at the beginning of date 1
    equityValTrans = (repmat(a(:,1),numBetas,1)-w(:,1))'*nansum(phi_e(:,:,1),1)'*p_e(1) - ...
        abar(1)*k*theta(1)*sbarTrans(1) - ...
        adjResource*acTrans(1) + ...
        qTrans(1);
    equityValDist = equityValTrans - equityVal - equityValDiffSS;
    
    % Real interest rate (only under sticky prices)
    if strcmp(prices,'sticky')
        MTrans = nan(1,Tbar);
        if fixedM == 0
            if fixedm == 0
                if isempty(dist1)
                    MTrans0 = mSS;
                else
                    MTrans0 = dist1.MTrans0;
                end
                if zlb 
                    MTrans(1) = (1+max((1-phiRho)*(phiInt + phiPi*PiTrans(1) + phiC*(cTrans(1)+g(1)-cTransSS(1))/cTransSS(1))+phiRho*(1/MTrans0-1),0)).^(-1);
                else
                    MTrans(1) = (1+(1-phiRho)*(phiInt + phiPi*PiTrans(1) + phiC*(cTrans(1)+g(1)-cTransSS(1))/cTransSS(1))+phiRho*(1/MTrans0-1)).^(-1);
                end
            else
                MTrans(1) = mSS/(1+PiTrans(2));
            end
        else
            MTrans(1) = mSS;
        end
        for time=2:Tbar
            if time > fixedM+1
                if time > fixedm+1
                    if zlb
                        MTrans(time) = (1+max((1-phiRho)*(phiInt + phiPi*PiTrans(time) + phiC*(cTrans(time)+g(time)-cTransSS(time))/cTransSS(time))+phiRho*(1/MTrans(time-1)-1),0)).^(-1);
                    else
                        MTrans(time) = (1+(1-phiRho)*(phiInt + phiPi*PiTrans(time) + phiC*(cTrans(time)+g(time)-cTransSS(time))/cTransSS(time))+phiRho*(1/MTrans(time-1)-1)).^(-1);
                    end
                else
                    MTrans(time) = mSS/(1+PiTrans(time+1));
                end
            else
                MTrans(time) = mSS;
            end
        end
        mTrans = MTrans.*[1+PiTrans(2:Tbar),1];
        mDist = d((mTrans-m).*weights,mDiffSS.*weights);
        mDiff = (mTrans-m)-mDiffSS;
        if zlb
            atZLB = sum(MTrans==1);
        else
            atZLB = 0;
        end
    else
        mTrans = m;
        MTrans = m; % just set to m as placeholder
    end
    
    % Vacancy posting
    newHireWeightingTrans = nan(gridPoints,eStates,Tbar);
    v_u_weightingTrans = nan(gridPoints,eStates,Tbar);
    for time=1:Tbar
        newHireWeightingTrans(:,:,time) = ((repmat((p(1,time)./...
            mbar(time).^(1-eta))',gridPoints,1).*...
            s(:,:,time).*phi_u_tilde(:,:,time))*Tue)/...
            nansum(squeeze(nansum(...
            repmat(reshape(p(1,time)./...
            repmat(mbar(time).^(1-eta),uStates,1),1,uStates),gridPoints,1).*...
            s(:,:,time).*phi_u_tilde(:,:,time),2)),1);
        v_u_weightingTrans(:,:,time) = (v_u(:,:,time).*...
            (repmat((p(1,time)./...
            mbar(time).^(1-eta))',gridPoints,1).*...
            s(:,:,time).*phi_u_tilde(:,:,time))/...
            nansum(squeeze(nansum(...
            repmat(reshape(p(1,time)./...
            repmat(mbar(time).^(1-eta),uStates,1),1,uStates),gridPoints,1).*...
            s(:,:,time).*phi_u_tilde(:,:,time),2)),1))*Tue;
    end    
    newHireWeightingTrans(isnan(newHireWeightingTrans)) = 0; % set to 0 if s = nan and phi = 0
    v_u_weightingTrans(isnan(v_u_weightingTrans)) = 0; % set to 0 if v_u,s = nan and phi = 0
    
    vacpostTrans = (mu.^(-1)).*abar.*k./q(theta,1:Tbar);
    sfNewHiresTrans = sum(sfTrans.*squeeze(nansum(newHireWeightingTrans,1)),1);
    dvacpostDist = d((vacpostTrans-sfNewHiresTrans).*weights,dvacpostDiffSS.*weights);
    dvacpostDiff = vacpostTrans - sfNewHiresTrans - dvacpostDiffSS;
    
    % Nash bargaining
    uPrime_eVecTrans = squeeze(nansum(uPrime(c_e).*newHireWeightingTrans,1));
    dvVecTrans = squeeze(nansum(v_e.*newHireWeightingTrans - v_u_weightingTrans,1));        
    swTrans = dvVecTrans./uPrime_eVecTrans;
    wTrans = iota*repmat(wSS',1,Tbar)+(1-iota)*(sfTrans + w - ((1-phi)/phi)*swTrans);
    dwageDist = d((wTrans-w).*repmat(weights,eStates,1),...
        dwageDiffSS.*repmat(weights,eStates,1)); 
    dwageDiff = wTrans-w-dwageDiffSS;
    
    % Report
    disp([tabs,'d(znextTrans,zsupplyTrans) = ',num2str(znextDist)]);
    disp([tabs,'d(cTrans+g,prodTrans) = ',num2str(goodsDist)]);
    disp([tabs,'d(tTrans,guess) = ',num2str(tDist)]);
    if strcmp(prices,'sticky')
        disp([tabs,'d(mTrans,guess) = ',num2str(mDist)]);
    end
    disp([tabs,'equityVal - guess = ',num2str(abs(equityValDist))]);
    disp([tabs,'dvacpost = ',num2str(dvacpostDist)]);
    disp([tabs,'dwage = ',num2str(dwageDist)]);
    if strcmp(prices,'sticky')
        disp([tabs,'theta(1:5) = ',num2str(theta(1:5))]); % to monitor progress
    end    
    disp(' ');
        
    % Assess convergence and update
    if isempty(H) % no Jacobian provided, so slowly update based on equilibrium errors
        if znextDist > eps_znext
            m = m + delta_m*znextDiff.*weights; 
            converged = 0;
        end 
        if tDist > eps_t
            t = t + delta_t*tDiff.*weights;
            converged = 0;
        end
        if strcmp(prices,'sticky')
            if mDist > eps_m
                mu = mu - delta_mu*mDiff.*weights;
                converged = 0;
            end
        end
        if abs(equityValDist) > eps_equityVal
            equityVal = equityVal + delta_equityVal*equityValDist; 
            converged = 0;
        end
        if dvacpostDist > eps_vac
            theta = theta - delta_theta*dvacpostDiff.*weights;
            converged = 0;
        end
        if dwageDist > eps_wage
            w = w + delta_w*dwageDiff.*repmat(weights,eStates,1);
            converged = 0;
        end
    else % use Jacobian to update using Newton method
        if strcmp(prices,'flex')
            if abs(equityValDist) > eps_equityVal || znextDist > eps_znext || dvacpostDist > eps_vac || tDist > eps_t || dwageDist > eps_wage
                if iota < 1
                    endog = [equityVal;m(1:Tbar-1)';theta(1:Tbar-1)';t(1:Tbar-1)';reshape(w(:,1:Tbar-1)',(Tbar-1)*eStates,1)];
                    endognew = endog - Hinv*...
                        [equityValTrans - equityVal - equityValDiffSS;...
                        znextTrans(1:Tbar-1)' - zsupplyTrans(1:Tbar-1)' - znextDiffSS(1:Tbar-1)';...
                        vacpostTrans(1:Tbar-1)' - sfNewHiresTrans(1:Tbar-1)' - dvacpostDiffSS(1:Tbar-1)';...
                        tTrans(1:Tbar-1)'-t(1:Tbar-1)'-tDiffSS(1:Tbar-1)';...
                        reshape(wTrans(:,1:Tbar-1)'-w(:,1:Tbar-1)'-dwageDiffSS(:,1:Tbar-1)',(Tbar-1)*eStates,1)];
                    endognew = dampen*endog + (1-dampen)*endognew;
                    w(:,1:Tbar-1) = reshape(endognew(3*Tbar-1:end),Tbar-1,eStates)';
                else
                    endog = [equityVal;m(1:Tbar-1)';theta(1:Tbar-1)';t(1:Tbar-1)'];
                    endognew = endog - Hinv*...
                        [equityValTrans - equityVal - equityValDiffSS;...
                        znextTrans(1:Tbar-1)' - zsupplyTrans(1:Tbar-1)' - znextDiffSS(1:Tbar-1)';...
                        vacpostTrans(1:Tbar-1)' - sfNewHiresTrans(1:Tbar-1)' - dvacpostDiffSS(1:Tbar-1)';...
                        tTrans(1:Tbar-1)'-t(1:Tbar-1)'-tDiffSS(1:Tbar-1)']; 
                    endognew = dampen*endog + (1-dampen)*endognew;
                end
                equityVal = endognew(1);
                m(1:Tbar-1) = endognew(2:Tbar)';
                theta(1:Tbar-1) = endognew(Tbar+1:2*Tbar-1)';
                t(1:Tbar-1) = endognew(2*Tbar:3*Tbar-2)';
                converged = 0;
            end
        else
            if ~pe
                if abs(equityValDist) > eps_equityVal || znextDist > eps_znext || dvacpostDist > eps_vac || tDist > eps_t || mDist > eps_m || dwageDist > eps_wage
                    if atZLB ~= lastAtZLB % adjust Jacobian to account for ZLB over atZLB periods
                        if atZLB >= 1   
                            Hzlb = H - ...
                                [zeros(4*Tbar-3,size(H,2));...
                                H((5+eStates)*Tbar-4-eStates+1:(5+eStates)*Tbar-4-eStates+atZLB,:);...
                                zeros((3+eStates)*(Tbar-1)-atZLB,size(H,2))];
                        else
                            Hzlb = H;
                        end
                        if iota < 1
                            Hinv = Hzlb([1:Tbar,2*Tbar:(5*Tbar-4+eStates*(Tbar-1))],...
                                1:(4*Tbar-3+eStates*(Tbar-1)))^(-1);
                        else
                            Hinv = Hzlb([1:Tbar,2*Tbar:(5*Tbar-4)],...
                                1:(4*Tbar-3))^(-1);            
                        end
                        lastAtZLB = atZLB;
                    end
                    if iota < 1
                        endog = [equityVal;m(fixedm+1:Tbar-1)';theta(1:Tbar-1)';t(1:Tbar-1)';mu(1:Tbar-1)';reshape(w(:,1:Tbar-1)',(Tbar-1)*eStates,1)];
                        endognew = endog - Hinv*...
                            [equityValTrans - equityVal - equityValDiffSS;...
                            znextTrans(1:Tbar-1)' - zsupplyTrans(1:Tbar-1)' - znextDiffSS(1:Tbar-1)';...
                            vacpostTrans(1:Tbar-1)' - sfNewHiresTrans(1:Tbar-1)' - dvacpostDiffSS(1:Tbar-1)';...
                            tTrans(1:Tbar-1)'-t(1:Tbar-1)'-tDiffSS(1:Tbar-1)';...
                            mTrans(fixedm+1:Tbar-1)'-m(fixedm+1:Tbar-1)'-mDiffSS(fixedm+1:Tbar-1)';...
                            reshape(wTrans(:,1:Tbar-1)'-w(:,1:Tbar-1)'-dwageDiffSS(:,1:Tbar-1)',(Tbar-1)*eStates,1)];
                        endognew = dampen*endog + (1-dampen)*endognew;                        
                        w(:,1:Tbar-1) = reshape(endognew(4*Tbar-2-fixedm:end),Tbar-1,eStates)';
                    else
                        endog = [equityVal;m(fixedm+1:Tbar-1)';theta(1:Tbar-1)';t(1:Tbar-1)';mu(1:Tbar-1)'];
                        endognew = endog - Hinv*...
                            [equityValTrans - equityVal - equityValDiffSS;...
                            znextTrans(1:Tbar-1)' - zsupplyTrans(1:Tbar-1)' - znextDiffSS(1:Tbar-1)';...
                            vacpostTrans(1:Tbar-1)' - sfNewHiresTrans(1:Tbar-1)' - dvacpostDiffSS(1:Tbar-1)';...
                            tTrans(1:Tbar-1)'-t(1:Tbar-1)'-tDiffSS(1:Tbar-1)';...
                            mTrans(fixedm+1:Tbar-1)'-m(fixedm+1:Tbar-1)'-mDiffSS(fixedm+1:Tbar-1)'];
                        endognew = dampen*endog + (1-dampen)*endognew;                                                
                    end
                    equityVal = endognew(1);
                    m(fixedm+1:Tbar-1) = endognew(2:Tbar-fixedm)';
                    theta(1:Tbar-1) = endognew(Tbar+1-fixedm:2*Tbar-1-fixedm)';
                    t(1:Tbar-1) = endognew(2*Tbar-fixedm:3*Tbar-2-fixedm)';
                    mu(1:Tbar-1) = endognew(3*Tbar-1-fixedm:4*Tbar-3-fixedm)';
                    converged = 0;
                end        
            else
                if tDist > eps_t || dwageDist > eps_wage
                    if iota < 1
                        endog = [t(1:Tbar-1)';reshape(w(:,1:Tbar-1)',(Tbar-1)*eStates,1)];
                        endognew = endog - Hinv*...
                            [tTrans(1:Tbar-1)'-t(1:Tbar-1)'-tDiffSS(1:Tbar-1)';...
                            reshape(wTrans(:,1:Tbar-1)'-w(:,1:Tbar-1)'-dwageDiffSS(:,1:Tbar-1)',(Tbar-1)*eStates,1)];
                        endognew = dampen*endog + (1-dampen)*endognew;                                                
                        w(:,1:Tbar-1) = reshape(endognew(Tbar:end),Tbar-1,eStates)';
                    else
                        endog = [t(1:Tbar-1)'];
                        endognew = endog - Hinv*...
                            [tTrans(1:Tbar-1)'-t(1:Tbar-1)'-tDiffSS(1:Tbar-1)'];              
                        endognew = dampen*endog + (1-dampen)*endognew;                          
                    end
                    t(1:Tbar-1) = endognew(1:Tbar-1)';
                    converged = 0;
                end             
            end
        end
    end
    
    % Save if converged
    if converged == 1
        if epsInd < numEps
            thisFilename = strcat(saveFilename,'_eps',num2str(epsInd),'.mat');
            converged = 0;
            epsInd = epsInd + 1;          
            if strcmp(prices,'flex')
                eps_znext = eps(epsInd,1);
                eps_t = eps(epsInd,2);
                eps_equityVal = eps(epsInd,3);
                eps_vac = eps(epsInd,4);
                eps_wage = eps(epsInd,5);
            else
                eps_znext = eps(epsInd,1);
                eps_t = eps(epsInd,2);
                eps_m = eps(epsInd,3);
                eps_equityVal = eps(epsInd,4);
                eps_vac = eps(epsInd,5);
                eps_wage = eps(epsInd,6);
            end
        else
            thisFilename = strcat(saveFilename,'.mat');
        end   
        phi_u_uStates = reshape(sum(phi_u,1),uStates,Tbar);
        % update varyingParams if anything changed
        varyingParams.znext_g = znext_g;
        
        % set dist2 
        dist2 = struct('z',z,'p_e_tilde1',p_e_tilde(2),...
            'phi_e_tilde1',phi_e_tilde(:,:,2),'phi_u_tilde1',phi_u_tilde(:,:,2),...
            'MTrans0',MTrans(1),...
            'initEquityVal',(repmat(a(:,2),numBetas,1)-w(:,2))'*nansum(phi_e(:,:,2),1)'*p_e(2) - ...
            abar(2)*k*theta(2)*sbarTrans(2) - adjResource*acTrans(2) + qTrans(2),...
            'zbar0',zbar(1),'p_e_tilde',p_e_tilde);
        if ~isempty(toName)
            save(thisFilename,'transParams','invariantParams','varyingParams',...
                'm','t','theta','equityVal','w','mu','b','znext_g',...
                'sbarTrans','znextTrans','zsupplyTrans','cTrans','prodTrans','yTrans',...
                'tTrans','mTrans','equityValTrans','vacpostTrans','sfNewHiresTrans',...
                'wTrans','swTrans','sfTrans',...
                'PiTrans','mu','MTrans','acTrans',...
                'p_e_tilde','p_e','p_u',...
                'avgW','phi_u_uStates');
        end
        
        if converged == 1 && exist('saveWorkspace') && saveWorkspace == 1
            save(thisFilename);
        end            
    end
    
    % Force end after `forceEnd' hours
    if toc(overallRuntime) > forceEnd    
        converged = 1;
        thisFilename = strcat(saveFilename,'_forceEnd.mat');
        phi_u_uStates = reshape(sum(phi_u,1),uStates,Tbar);
        % update varyingParams if anything changed
        varyingParams.znext_g = znext_g;
        if ~isempty(toName)
            save(thisFilename,'transParams','invariantParams','varyingParams',...
                'm','t','theta','equityVal','w','mu','b','znext_g',...
                'sbarTrans','znextTrans','zsupplyTrans','cTrans','prodTrans','yTrans',...
                'tTrans','mTrans','equityValTrans','vacpostTrans','sfNewHiresTrans',...
                'wTrans','swTrans','sfTrans',...
                'PiTrans','mu','MTrans','acTrans',...
                'p_e_tilde','p_e','p_u',...
                'avgW','phi_u_uStates'); 
        end
    end      
    
    disp([tabs,'Time elapsed: = ',num2str(toc(singleIter))]);
    disp(' ');
end

toc(overallRuntime);