function [m,t,theta,sf] = StationaryRCE(params,accuracy,starting,RCEName,numTabs,ident)
%% StationaryRCE.m
% Characterizes stationary RCE given
%  params: a struct holding all economic parameters
%  accuracy: a struct holding gridPoints, eps, upper, eps_v, finerStatGrid, 
%   eps_phi, and forceEnd
%  starting: a vector holding starting [m,t,theta,sf]
%  RCEName: name to save to ([] means nothing is saved)
%  numTabs: number of tabs to indent when printing to screen
%  ident: indicator variable = 1 if code is being run to characterize
%   identification of parameters, rather than precisely solving for RCE

overallRuntime = tic;

tabs = '.';
for i=1:numTabs
    tabs = sprintf('   %s',tabs);
end

%% Set economic parameters

global numBetas betas xi sigma ...
    numP sigmaP rhoP numT sigmaT numA w Ta ...
    numUI Taui Tua dbarbar dbar rr b zeta phi ...
    k deltaBar eps_delta_a eps_delta_beta mbar eta lambda ...
    abar a aP ...
    zbar znext_g ...
    omega0 omega1 omega2 pEpsilon equitycutoff...
    epsilon tauR uimax;

% Tastes 
numBetas = params.numBetas;
betas = params.betas;
xi = params.xi;
sigma = params.sigma;

% Wages
numP = params.numP;
sigmaP = params.sigmaP;
rhoP = params.rhoP;
numT = params.numT;
sigmaT = params.sigmaT;
numA = params.numA;
phi = params.phi;
Ta = params.Ta;

% UI 
numUI = params.numUI;
Taui = params.Taui;
Tua = params.Tua;
dbarbar = params.dbarbar;
dbar = params.dbar;
rr = params.rr;
zeta = params.zeta;
uimax = params.uimax;

% Technology
abar = params.abar;
a = params.a;
aP = params.aP;
k = params.k;
deltaBar = params.deltaBar;
eps_delta_a = params.eps_delta_a;
eps_delta_beta = params.eps_delta_beta;
mbar = params.mbar;
eta = params.eta;
lambda = params.lambda;

% Assets 
zbar = params.zbar;
znext_g = params.znext_g;

% Non-labor endowment
omega0 = params.omega0;
omega1 = params.omega1;
omega2 = params.omega2;
pEpsilon = params.pEpsilon;

% Equity share cutoff (only used to compute moments regarding disposable income)
equitycutoff = params.equitycutoff;

% Retailers
epsilon = params.epsilon;
tauR = params.tauR;

%% Set tuning parameters

eps = accuracy.eps;
eps_goods = eps(1);
eps_t = eps(2);
eps_vac = eps(3);
eps_bargain = eps(4);

if ~ident
    delta_m = 0.0005;
    delta_t = 0.25; 
    delta_theta = 0.25;
else
    delta_m = 0.001;
    delta_t = 0.5; 
    delta_theta = 2;    
end

gridPoints = accuracy.gridPoints;
upper = accuracy.upper;

tol = 1e-12; % meant to weed out computational error

%% Other preliminaries before starting algorithm
% Speed up algorithm despite more heterogeneity by keeping
% policies and value functions in columns indexed by agent's state
% variables (except for assets, which index the rows).  Throughout, will
% keep `order of heterogeneity' consistent:
%  For the employed:
%   Beta (#=numBetas)
%   Productivity (#=numA)
%  For the unemployed:
%   Beta (#=numBetas)
%   UI schedule (#=numUI)
%   Transitory component of prod (#=numT)
%   Duration of unemployment (#=dbarbar+1)

% Define marginal utility of real income for employed
uPrime = @(c) c.^(-sigma);

% Define functions returning job-finding and vacancy-filling probabilities
lambdaAdj = nan(1,dbarbar+1);
lambdaAdj(1:(dbar+1)) = lambda*(0:dbar);
lambdaAdj((dbar+2):(dbarbar+1)) = ones(1,dbarbar-dbar)*lambda*dbar;
p = @(t) repmat((mbar^(1-eta))*exp(lambdaAdj)*(t)^(eta),1,numT*numUI*numBetas); % will be of size 1 x [dbarbar+1]*numT*numUI*numBetas
p0 = @(t) (mbar^(1-eta))*(t)^(eta); 
q = @(t) (mbar^(1-eta))*t^(eta-1); 

% Now define matrices {eState,uState}->{eState,uState} 
%  transition probabilities

% - first by defining the number of employed and unemployed
eStates = numA*numBetas;
uStates = (dbarbar+1)*numT*numUI*numBetas;

% - then by computing number of beta and wage for each state when employed
eBetas = ceil((1:numA*numBetas)/(numA));
eA = repmat(1:numA,1,numBetas);
% - and then by defining transition functions for next state
ee = @(thisState) cumsum([floor((thisState-1)/numA)*numA+1;...
   ones(numA-1,length(thisState))],1);
ee_prob = @(thisState) Ta(eA(thisState)',:)';
eu = @(thisState) cumsum([(eBetas(thisState)-1).*((dbarbar+1)*numT*numUI)+1;...
    ones(numT*numUI-1,length(thisState))*(dbarbar+1)],1);
eu_prob = @(thisState) kron(Taui(eA(thisState),:)',pEpsilon');
Tee = sparse(reshape(repmat(1:eStates,size(ee(1:eStates),1),1),eStates*size(ee(1:eStates),1),1),...
    reshape(ee(1:eStates),eStates*size(ee(1:eStates),1),1),...
    reshape(ee_prob(1:eStates),eStates*size(ee(1:eStates),1),1),...
    eStates,eStates);
Teu = sparse(reshape(repmat(1:eStates,size(eu(1:eStates),1),1),eStates*size(eu(1:eStates),1),1),...
    reshape(eu(1:eStates),eStates*size(eu(1:eStates),1),1),...
    reshape(eu_prob(1:eStates),eStates*size(eu(1:eStates),1),1),...
    eStates,uStates);

% - analogously, by computing number of beta, UI schedule, 
%   endowment trans component, and duration for each state when unemployed
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

% Now define separation rate for each state of employment
deltaVec = nan(1,eStates);
for iBeta=1:numBetas
    deltaVec((iBeta-1)*numA+1:iBeta*numA) = ...
        max(deltaBar + eps_delta_a*(log(aP)-log(abar)) + eps_delta_beta*(betas(iBeta)-mean(betas)),1e-6);
end

% Finally initialize asset space z
z = zbar + (logspace(0,3,gridPoints)-1)'/999*(upper-zbar); 

%% Run algorithm to characterize stationary RCE

% Keep track of convergence
converged = 0;

% Step 0) Choose initial aggregates
m = starting(1); 
t = starting(2);
theta = starting(3);
sf = starting(4:eStates+3);

% Solve for steady-state gross mark-up
mu = (epsilon/(epsilon-1))*(1+tauR); 

while(converged == 0)
    toDisp = [tabs,'***NEW m = ',num2str(m),' t = ', num2str(t),' theta = ', num2str(theta)];
    disp(toDisp);
    for i_e=1:eStates
        disp([tabs,'*   sf(',num2str(i_e),') = ',num2str(sf(i_e))]);
    end    

    % Solve for wages consistent with conjectured surplus
    w = (1/mu)*repmat(a,1,numBetas) - sf + (1-deltaVec).*m.*(Tee*sf')';
    wP = pEpsilon*reshape(w,numT,numP*numBetas);
     
    % Solve for UI consistent with those wages
    [b,~] = genTransfers(numBetas,numP,numT,dbarbar,dbar,rr*(1-omega0),...
        w,wP,omega1,omega2,zeta,uimax);
    
    % Define y_e and y_u accordingly
    y_e = w - t; % a 1 x numT*numP*numBetas vector
    y_u = b; % a 1 x (dbarbar+1)*numT*numUI*numBetas vector
       
    % Characterize optimal policies and value functions until
    % znext_e doesn't hit the artifical upper bound
    upper = accuracy.upper;
    upperBoundHit = 1;
    while upperBoundHit == 1
        [~,v_e,v_u,s,~,znext_e,~,znext_u] = ...
            statPoliciesEGM(m,p(theta),deltaVec,Tee,Teu,Tue,Tuu,...
            y_e,y_u,betas(eBetas),betas(uBetas),accuracy.eps_v,~ident,z);

        eSaver = numBetas*numA;
        % theoretically, this will hold the index of the employed
        %  agent who saves the most at least at high levels of wealth (highest
        %  beta and highest w)
        if z(gridPoints) - znext_e(gridPoints,eSaver) > tol
            upperBoundHit = 0;
        else
            slope = (znext_e(gridPoints,eSaver) - znext_e(gridPoints-1,eSaver))/...
                (z(gridPoints)-z(gridPoints-1));
            interpolatedUpper = (znext_e(gridPoints,eSaver) - slope*z(gridPoints))/...
                (1-slope); % so that by linear interpolation, 
            % znext_e(interpolatedUpper,eSaver) should = interpolatedUpper
            
            % Then extend z to 1.1*interpolatedUpper
            z = zbar + (logspace(0,3,gridPoints)-1)'/999*(1.1*interpolatedUpper-zbar); 
            
            disp([tabs,'Iterative algorithm upper bound hit: expanding it to ',...
                num2str(1.1*interpolatedUpper)]);           
        end
    end
    
    % Alert if p(theta)s = 1 for some unemployed agents
    if sum(sum(repmat(p(theta),gridPoints,1).*s >= 1)) > 0
        disp([tabs,'Some agents at max effort!']);
        pause;
    end
    
    % To eliminate miniscule computational error, enforce:
    znext_e(znext_e <= zbar) = zbar;
    znext_u(znext_u <= zbar) = zbar;
    
    % Construct zStat
    if accuracy.finerStatGrid == 1
        % Define gridLength to be used in constructing a (weakly) finer grid
        gridLength = .9*min(z(gridPoints) - ...
            znext_e(gridPoints,eSaver),...
            min(z(2:gridPoints)-z(1:(gridPoints-1)))); 
        % the above requires that znext_e(upper bound,i) < upper bound for all i, 
        %  which is necessary for an invariant distribution to exist
        
        % Construct grid
        zStat = linspace(z(1),z(gridPoints),ceil((z(gridPoints)-z(1))/gridLength))';
        
        % An alternative which introduces denser gridpoints at end of space
        %  to ease transitional algorithm down the line
        gridLength = .9*min(z(gridPoints) - ...
            znext_e(gridPoints,eSaver),...
            z(gridPoints)-z(gridPoints-1)); 
        zStat = [z(1:(gridPoints-2));...
            linspace(z(gridPoints-1),z(gridPoints),ceil((z(gridPoints)-z(gridPoints-1))/gridLength))'];
    else
        zStat = z;
    end
    
    % Characterize implied invariant distribution
    [p_e_tilde,phi_e_tilde,phi_u_tilde,p_e,phi_e,phi_u,...
        znext_eStat,znext_uStat,sStat] = ...
        obtainInvariantSparse(z,s,znext_e,znext_u,p(theta),deltaVec,...
        Tee,Teu,Tue,Tuu,zStat,accuracy.eps_phi,~ident);
    gridpointsStat = length(zStat);
    
    % Interpolate along zStat, which may be weakly finer than z
    v_eStat = interp1q(z,v_e,zStat); 
    c_eStat = repmat(y_e,gridpointsStat,1) + repmat(zStat,1,eStates) - m*znext_eStat;
    
    v_uStat = interp1q(z,v_u,zStat);
    c_uStat = repmat(y_u,gridpointsStat,1) + repmat(zStat,1,uStates) - m*znext_uStat;
        
    % Define worker surplus relevant for Nash bargaining
    uPrime_eVec = sum(uPrime(c_eStat).*((repmat(p(1)/(mbar^(1-eta)),gridpointsStat,1).*...
        sStat.*phi_u_tilde)*Tue),1)/...
        sum(sum(repmat(p(1)/(mbar^(1-eta)),gridpointsStat,1).*...
        sStat.*phi_u_tilde));
    dvVec = sum(v_eStat.*((repmat(p(1)/(mbar^(1-eta)),gridpointsStat,1).*...
        sStat.*phi_u_tilde)*Tue) - ...
        (v_uStat.*repmat(p(1)/(mbar^(1-eta)),gridpointsStat,1).*...
        sStat.*phi_u_tilde)*Tue,1)/...
        sum(sum(repmat(p(1)/(mbar^(1-eta)),gridpointsStat,1).*...
        sStat.*phi_u_tilde));
    sw = dvVec./uPrime_eVec;   
    
    % Step 3) Compute macro aggregates, assess convergence and update 
    
    converged = 1;
    
    % Aggregate borrowing and supply of assets 
    znextStat = sum(sum(znext_eStat.*phi_e*p_e)) + ...
        sum(sum(znext_uStat.*phi_u*(1-p_e)));
    sbar = sum(sum(repmat(p(1)/(mbar^(1-eta)),gridpointsStat,1).*...
        sStat.*phi_u_tilde*(1-p_e_tilde)));
    qStat = (m/(1-m))*((repmat(a,1,numBetas)-w)*sum(phi_e,1)'*p_e - abar*k*theta*sbar);
    zsupplyStat = -znext_g + (1/m)*qStat;
    % Lump-sum tax
    tStat = (1/(p_e))*((m-1)*znext_g + ...
        (1-p_e)*sum(sum(phi_u,1).*b));            
    % Aggregate consumption and production 
    cStat = sum(sum(c_eStat.*phi_e*p_e)) + ...
        sum(sum(c_uStat.*phi_u*(1-p_e)));  
    prodStat = a(eA)*sum(phi_e,1)'*p_e - abar*k*theta*sbar;
    % Optimal vacancy posting
    sfNewHiresStat = sum(sf.*sum(((repmat(p(1)/(mbar^(1-eta)),gridpointsStat,1).*...
        sStat.*phi_u_tilde)*Tue),1)/...
        sum(sum(repmat(p(1)/(mbar^(1-eta)),gridpointsStat,1).*...
        sStat.*phi_u_tilde)));
    dvacpost = (1/mu)*abar*k/q(theta)-sfNewHiresStat; 
    % Nash bargaining
    dbargain = sw-(phi/(1-phi))*sf;    
    
    % Report
    disp([tabs,'znext - zsupply = ',num2str(znextStat-zsupplyStat)]);
    disp([tabs,'cons - prod = ',num2str(cStat-prodStat)]);
    disp([tabs,'t - guess = ', num2str(tStat-t)]);
    disp([tabs,'d vac post = ', num2str(dvacpost)]);
    disp([tabs,'worker surp - firm surp = ']);
    for i_e=1:eStates
        disp([tabs,num2str(dbargain(i_e))]);
    end      
    disp(' ');
    
    % Assess convergence and update
    if abs(cStat-prodStat) > eps_goods
        m = m - delta_m*(cStat-prodStat);
        converged = 0;
    end
    
    if abs(tStat - t) > eps_t
        t = t + delta_t*(tStat - t);
        converged = 0;
    end
    
    if abs(dvacpost) > eps_vac
        theta = theta - delta_theta*dvacpost;
        converged = 0;
    end
    
    if max(abs(dbargain)) > eps_bargain
        sf = sf + (1-phi)/phi*dbargain;
        converged = 0;
    end    
    
    % Force end after `forceEnd' time, if specified
    if isfield(accuracy,'forceEnd') && toc(overallRuntime) > accuracy.forceEnd    
        converged = 1;
    end
end

%% Compute additional moments of interest

% Simulate agents for use in computing several moments below
[zInd,eTildeStatus,eStateTildeInd,uTildeStatus,uStateTildeInd,...
    eStatus,eStateInd,uStatus,uStateInd] = ...
    simStationaryAgents(p_e_tilde,phi_e_tilde,phi_u_tilde,zStat,...
    znext_eStat,znext_uStat,sStat,p(theta),deltaVec,...
    Tee,Teu,Tue,Tuu,2000000,73,0,0); % higher sample necessary for EU coefficients
[zIndSep,eTildeStatusSep,eStateTildeIndSep,uTildeStatusSep,uStateTildeIndSep,...
    eStatusSep,eStateIndSep,uStatusSep,uStateIndSep] = ...
    simStationaryAgents(p_e_tilde,phi_e_tilde,phi_u_tilde,zStat,...
    znext_eStat,znext_uStat,sStat,p(theta),deltaVec,...
    Tee,Teu,Tue,Tuu,200000,73,0,1);          

% (Annual) interest rate

r = (1/m-1)*12;

% Patterns of wealth

pdf_z = sum(p_e*phi_e,2) + ...
    sum((1-p_e)*phi_u,2);
cdf_z = cumsum(pdf_z);

ind_50 = find(cdf_z >= 0.5,1);   

wealth_to_moincome = sum(pdf_z.*zStat)/prodStat; 

dwealth_ue = (sum(zStat.*sum(phi_u,2)) - sum(zStat.*sum(phi_e,2))) / ...
    sum(zStat.*pdf_z);    
dwealth_ue_to_moincome = (sum(zStat.*sum(phi_u,2)) - sum(zStat.*sum(phi_e,2))) / ...
    prodStat; 

z0 = find(zStat >= 0,1);
share_wealth_lt0 = cdf_z(z0-1);

% Patterns of household income

if zeta < 1
    frac_unemp_UI = sum(sum(phi_u(:,uUI>numP&uDur<=dbar)));
else
    frac_unemp_UI = sum(sum(phi_u(:,uDur<=dbar))); 
end

% first define index of persistent wage component for each eState
eWP = reshape(repmat(ceil(numT/2):numT:...
    (numBetas-1)*numA+(numP-1)*numT+ceil(numT/2),numT,1),...
    1,eStates);

% - income of unemp HH receiving benefits / pre-job-loss (persist inc, among UI exhaustees only)
if zeta < 1
    unempExhaust = (eStatusSep(:,1)==1)&(uStatusSep(:,dbar+2)==1)&(uDur(uStateIndSep(:,dbar+2))==(dbar+1))'&(uUI(uStateIndSep(:,dbar+2))>numP)';
else
    unempExhaust = (eStatusSep(:,1)==1)&(uStatusSep(:,dbar+2)==1)&(uDur(uStateIndSep(:,dbar+2))==(dbar+1))';
end 
hhinc_ue_UI = mean(y_u(uStateIndSep(unempExhaust,2)) ./ y_e(eWP(eStateIndSep(unempExhaust,1))));

% - income of unemp HH not receiving benefits / pre-job-loss 
hhinc_ue_NoUI = mean(y_u(uStateIndSep(unempExhaust,dbar+2)) ./ y_e(eWP(eStateIndSep(unempExhaust,1))));

% Unemployment rate, overall and by duration

ur = 1-p_e;

share_ur_lt = sum(sum(phi_u(:,uDur>6))); 

% Duration elasticity to potential benefit duration, where I focus on
%  the responses *among UI recipients* only (and count exits to 
%  retirement as a completed UI spell).  This reflects my understanding 
%  that the literature estimating these elasticities measures unemployment 
%  using weeks of UI claimed (e.g., Meyer (1990) or Card and Levine (2000)), 
%  and thus reflects responses conditional on UI receipt.  Note also
%  that the code below assumes dbar > 0.
%  Note furthermore that the duration elasticity accounts for any change in
%    endowment in response to the change in UI policy    

znext_uStat(znext_uStat <= zbar) = zbar;
dur = simDurationAlt(phi_u(:,uDur==1)/sum(sum(phi_u(:,uDur==1))),...
    uDur==1,zStat,znext_uStat,sStat,p(theta),Tuu,0);

bHiUIDur = reshape(b,dbarbar+1,numT*numUI*numBetas);
bHiUIDur(dbar+1,:) = bHiUIDur(dbar,:);
bHiUIDur = reshape(bHiUIDur,1,(dbarbar+1)*numT*numUI*numBetas);
if zeta == 1
    unempOfInterest = (uDur==1);
else
    unempOfInterest = (uDur==1) & (uUI > numP);
end
y_uHiUIDur = bHiUIDur; % a 1 x [dbarbar+1]*numT*numUI*numBetas vector
y_e = w - t; % a 1 x numA*numBetas or 1 x numA*numBetas vector 

[durHiUIDur,~,~,~,~,~] = computeDurAltUI(zStat,...
    phi_u(:,unempOfInterest)/sum(sum(phi_u(:,unempOfInterest))),...
    unempOfInterest,m,p(theta),deltaVec,Tee,Teu,Tue,Tuu,y_e,y_uHiUIDur,...
    betas(eBetas),betas(uBetas),accuracy.eps_v,0,1e-12);
elast_dur = ((durHiUIDur - dur)/dur)/(1/dbar);

conv_tightness = theta*(sbar/((1-p_e_tilde)));
mowage_per_vac = (k/q(theta))/(1-omega0);

gridpointsStat = length(zStat);

% Ganong-Noel consumption decline through exhaustion

if zeta < 1       
    gnSample = (eStatusSep(:,1)==1)&(uStatusSep(:,dbar+1)==1)&(uDur(uStateIndSep(:,dbar+1))==dbar)'&(uUI(uStateIndSep(:,2))>numP)';        
else
    gnSample = (eStatusSep(:,1)==1)&(uStatusSep(:,dbar+1)==1)&(uDur(uStateIndSep(:,dbar+1))==dbar)';  
end
 % includes agents who were employed at date 1, unemployed and receiving UI
 %  at date 2, and remain unemployed in the month UI benefits expire 
 %  (but could have found a job thereafter!)
gnCe = c_eStat((eStateIndSep(gnSample,1)-1)*(gridpointsStat) + zIndSep(gnSample,1));  
gnCuReceipt = mean(c_uStat((uStateIndSep(gnSample,2:(dbar+1))-1)*(gridpointsStat) + zIndSep(gnSample,2:(dbar+1))),2);
gnCuPostExhaust = ...
    c_uStat((uStateIndSep(gnSample,dbar+2)-1)*(gridpointsStat) + zIndSep(gnSample,dbar+2)).*...
    uStatusSep(gnSample,dbar+2) + ...
    c_eStat((eStateIndSep(gnSample,dbar+2)-1)*(gridpointsStat) + zIndSep(gnSample,dbar+2)).*...
    eStatusSep(gnSample,dbar+2);
gn17receipt = mean(gnCuReceipt./gnCe);
gn17exhaust = mean(gnCuPostExhaust./gnCe);

% GN change in spending / change in income at exhaustion (conditional 
%  only on those remaining unemployed after exhaustion, and only focused on
%  one-month changes)
if zeta < 1
    gnSample2 = (eStatusSep(:,1)==1)&(uStatusSep(:,dbar+1)==1)&(uStatusSep(:,dbar+2)==1)&(uDur(uStateIndSep(:,dbar+1))==dbar)'&(uUI(uStateIndSep(:,2))>numP)';        
else
    gnSample2 = (eStatusSep(:,1)==1)&(uStatusSep(:,dbar+1)==1)&(uStatusSep(:,dbar+2)==1)&(uDur(uStateIndSep(:,dbar+1))==dbar)';         
end

gnCuPreExhaust2 = c_uStat((uStateIndSep(gnSample2,dbar+1)-1)*(gridpointsStat) + zIndSep(gnSample2,dbar+1));
gnCuPostExhaust2 = c_uStat((uStateIndSep(gnSample2,dbar+2)-1)*(gridpointsStat) + zIndSep(gnSample2,dbar+2));
gndCu2 = mean(gnCuPostExhaust2 - gnCuPreExhaust2);
gnYuPreExhaust2 = y_u(uStateIndSep(gnSample2,dbar+1))';
gnYuPostExhaust2 = y_u(uStateIndSep(gnSample2,dbar+2))';
gndYu2 = mean(gnYuPostExhaust2 - gnYuPreExhaust2);
gn17dCdY2 = gndCu2/gndYu2;     

% MPCs out of ~$500 rebate 
% - first compute expected consumption next month as a function of today's states
[c_eStat1,c_uStat1] = computeExpectedCons(zStat,sStat,znext_eStat,znext_uStat,p(theta),deltaVec,...
        Tee,Teu,Tue,Tuu,c_eStat,c_uStat);

% - then compute expected consumption in the month following next as a function of
%  today's states
[c_eStat2,c_uStat2] = computeExpectedCons(zStat,sStat,znext_eStat,znext_uStat,p(theta),deltaVec,...
        Tee,Teu,Tue,Tuu,c_eStat1,c_uStat1);

% - finally aggregate into expected consumption over next 3 months        
c_eQ = c_eStat + c_eStat1 + c_eStat2;
c_uQ = c_uStat + c_uStat1 + c_uStat2;

dz = (500/6761)*prodStat; % using 6761 as estimate of avg HH income in 2004 dollars among emp and unemp HH

zAfterRebate = zStat(1:(gridpointsStat-1)) + dz;

% - by employment status and duration
c_eQAfterRebate = interp1q(zStat,c_eQ,zAfterRebate);
dc_eQ = c_eQAfterRebate - c_eQ(1:(gridpointsStat-1),:);
mpc_eQ = dc_eQ/dz;
c_uQAfterRebate = interp1q(zStat,c_uQ,zAfterRebate);
dc_uQ = c_uQAfterRebate - c_uQ(1:(gridpointsStat-1),:);
mpc_uQ = dc_uQ/dz;

mpc_eQAvg = sum(sum(mpc_eQ .* phi_e(1:(gridpointsStat-1),:)))/...
    sum(sum(phi_e(1:(gridpointsStat-1),:)));
mpc_uQAvg = sum(sum(mpc_uQ .* phi_u(1:(gridpointsStat-1),:)))/...
    sum(sum(phi_u(1:(gridpointsStat-1),:)));
mpc_uQAvgByDur = nan(1,dbarbar+1);
phi_uByDur = nan(1,dbarbar+1);
for d=1:(dbarbar+1)
   mpc_uQAvgByDur(d) = sum(sum(mpc_uQ(:,uDur==d).*phi_u(1:(gridpointsStat-1),uDur==d)))/...
       sum(sum(phi_u(1:(gridpointsStat-1),uDur==d))); 
   phi_uByDur(d) = sum(sum(phi_u(1:(gridpointsStat-1),uDur==d)));
end

mpc_QAvg = p_e*mpc_eQAvg + (1-p_e)*mpc_uQAvg;
mpc_uSTQAvg = sum(mpc_uQAvgByDur(1:3).*phi_uByDur(1:3))/sum(phi_uByDur(1:3));
mpc_uMTQAvg = sum(mpc_uQAvgByDur(4:6).*phi_uByDur(4:6))/sum(phi_uByDur(4:6));
mpc_uLTQAvg = sum(mpc_uQAvgByDur(7:(dbarbar+1)).*phi_uByDur(7:(dbarbar+1)))/sum(phi_uByDur(7:(dbarbar+1)));

% - by pre-tax income
[equityShares] = ...
    initializeEquity(equitycutoff,...
    (1/m)*qStat,zStat,[p_e_tilde*phi_e_tilde,(1-p_e_tilde)*phi_u_tilde]);
ydis_e = repmat(y_e,gridpointsStat,1) + repmat(equityShares.*zStat.*(1-m),1,eStates);
ydis_u = repmat(y_u,gridpointsStat,1) + repmat(equityShares.*zStat.*(1-m),1,uStates);
    
pretaxy_e = t + ydis_e;
pretaxy_u = ydis_u;

[pretaxy,~,indpretaxy] = unique([pretaxy_e(:);pretaxy_u(:)]);
temp_phi = [p_e*phi_e(:);(1-p_e)*phi_u(:)];
pdf_pretaxy = zeros(length(pretaxy),1);
for temp_ind = 1:length(pretaxy)
    pdf_pretaxy(temp_ind) = sum(temp_phi(indpretaxy==temp_ind));
end
cdf_pretaxy = cumsum(pdf_pretaxy);

ind_e_lt75k = (pretaxy_e(1:end-1,:) <= 0.92*prodStat);  
ind_u_lt75k = (pretaxy_u(1:end-1,:) <= 0.92*prodStat);

mpc_lt75k = (p_e*sum(sum(mpc_eQ(ind_e_lt75k).*phi_e([ind_e_lt75k;false(1,eStates)]))) + ...
    (1-p_e)*sum(sum(mpc_uQ(ind_u_lt75k).*phi_u([ind_u_lt75k;false(1,uStates)]))))/...
    (p_e*sum(sum(phi_e([ind_e_lt75k;false(1,eStates)]))) + ...
    (1-p_e)*sum(sum(phi_u([ind_u_lt75k;false(1,uStates)]))));   
mpc_gt75k = (p_e*sum(sum(mpc_eQ(~ind_e_lt75k).*phi_e([~ind_e_lt75k;false(1,eStates)]))) + ...
    (1-p_e)*sum(sum(mpc_uQ(~ind_u_lt75k).*phi_u([~ind_u_lt75k;false(1,uStates)]))))/...
    (p_e*sum(sum(phi_e([~ind_e_lt75k;false(1,eStates)]))) + ...
    (1-p_e)*sum(sum(phi_u([~ind_u_lt75k;false(1,uStates)])))); 

% - by wealth quintile
ind_20 = find(cdf_z >= 0.2,1);
ind_40 = find(cdf_z >= 0.4,1);
ind_60 = find(cdf_z >= 0.6,1);
ind_80 = find(cdf_z >= 0.8,1);

ind_ZQ1 = [ones(ind_20-1,1);...
    (0.2-cdf_z(ind_20-1))/(cdf_z(ind_20)-cdf_z(ind_20-1));...
    zeros(gridpointsStat-ind_20,1)];
ind_ZQ2 = [zeros(ind_20-1,1);...
    (cdf_z(ind_20)-0.2)/(cdf_z(ind_20)-cdf_z(ind_20-1));...
    ones(ind_40-1-ind_20,1);...
    (0.4-cdf_z(ind_40-1))/(cdf_z(ind_40)-cdf_z(ind_40-1));...
    zeros(gridpointsStat-ind_40,1)];
ind_ZQ3 = [zeros(ind_40-1,1);...
    (cdf_z(ind_40)-0.4)/(cdf_z(ind_40)-cdf_z(ind_40-1));...
    ones(ind_60-1-ind_40,1);...
    (0.6-cdf_z(ind_60-1))/(cdf_z(ind_60)-cdf_z(ind_60-1));...
    zeros(gridpointsStat-ind_60,1)]; 
ind_ZQ4 = [zeros(ind_60-1,1);...
    (cdf_z(ind_60)-0.6)/(cdf_z(ind_60)-cdf_z(ind_60-1));...
    ones(ind_80-1-ind_60,1);...
    (0.8-cdf_z(ind_80-1))/(cdf_z(ind_80)-cdf_z(ind_80-1));...
    zeros(gridpointsStat-ind_80,1)];
ind_ZQ5 = [zeros(ind_80-1,1);...
    (cdf_z(ind_80)-0.8)/(cdf_z(ind_80)-cdf_z(ind_80-1));...
    ones(gridpointsStat-ind_80,1)];
    
mpc_ZQ1 = (p_e*sum(sum(repmat(ind_ZQ1(1:end-1),1,eStates).*mpc_eQ.*phi_e(1:end-1,:))) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ1(1:end-1),1,uStates).*mpc_uQ.*phi_u(1:end-1,:))))/...
    (p_e*sum(sum(repmat(ind_ZQ1(1:end-1),1,eStates).*phi_e(1:end-1,:))) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ1(1:end-1),1,uStates).*phi_u(1:end-1,:))));
mpc_ZQ2 = (p_e*sum(sum(repmat(ind_ZQ2(1:end-1),1,eStates).*mpc_eQ.*phi_e(1:end-1,:))) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ2(1:end-1),1,uStates).*mpc_uQ.*phi_u(1:end-1,:))))/...
    (p_e*sum(sum(repmat(ind_ZQ2(1:end-1),1,eStates).*phi_e(1:end-1,:))) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ2(1:end-1),1,uStates).*phi_u(1:end-1,:))));
mpc_ZQ3 = (p_e*sum(sum(repmat(ind_ZQ3(1:end-1),1,eStates).*mpc_eQ.*phi_e(1:end-1,:))) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ3(1:end-1),1,uStates).*mpc_uQ.*phi_u(1:end-1,:))))/...
    (p_e*sum(sum(repmat(ind_ZQ3(1:end-1),1,eStates).*phi_e(1:end-1,:))) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ3(1:end-1),1,uStates).*phi_u(1:end-1,:))));
mpc_ZQ4 = (p_e*sum(sum(repmat(ind_ZQ4(1:end-1),1,eStates).*mpc_eQ.*phi_e(1:end-1,:))) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ4(1:end-1),1,uStates).*mpc_uQ.*phi_u(1:end-1,:))))/...
    (p_e*sum(sum(repmat(ind_ZQ4(1:end-1),1,eStates).*phi_e(1:end-1,:))) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ4(1:end-1),1,uStates).*phi_u(1:end-1,:))));   
mpc_ZQ5 = (p_e*sum(sum(repmat(ind_ZQ5(1:end-1),1,eStates).*mpc_eQ.*phi_e(1:end-1,:))) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ5(1:end-1),1,uStates).*mpc_uQ.*phi_u(1:end-1,:))))/...
    (p_e*sum(sum(repmat(ind_ZQ5(1:end-1),1,eStates).*phi_e(1:end-1,:))) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ5(1:end-1),1,uStates).*phi_u(1:end-1,:))));

% Annual MPC out of one month income

c_eA = c_eQ;
c_uA = c_uQ;  
perAhead = 3;
c_eStatToAvg = c_eStat2;
c_uStatToAvg = c_uStat2;         
while perAhead <= 11
    [c_eStatNext,c_uStatNext] = computeExpectedCons(zStat,sStat,znext_eStat,znext_uStat,p(theta),deltaVec,...
            Tee,Teu,Tue,Tuu,c_eStatToAvg,c_uStatToAvg);
    c_eA = c_eA + c_eStatNext;
    c_uA = c_uA + c_uStatNext;
    c_eStatToAvg = c_eStatNext;
    c_uStatToAvg = c_uStatNext;             
    perAhead = perAhead + 1;
end

dzA = prodStat;

zAfterRebateA = zStat(1:(gridpointsStat-1)) + dzA;

c_eAAfterRebate = interp1q(zStat,c_eA,zAfterRebateA);
dc_eA = c_eAAfterRebate - c_eA(1:(gridpointsStat-1),:);
mpc_eA = dc_eA/dzA;
c_uAAfterRebate = interp1q(zStat,c_uA,zAfterRebateA);
dc_uA = c_uAAfterRebate - c_uA(1:(gridpointsStat-1),:);
mpc_uA = dc_uA/dzA;

mpc_eAAvg = sum(sum(mpc_eA .* phi_e(1:(gridpointsStat-1),:)))/...
    sum(sum(phi_e(1:(gridpointsStat-1),:)));
mpc_uAAvg = sum(sum(mpc_uA .* phi_u(1:(gridpointsStat-1),:)))/...
    sum(sum(phi_u(1:(gridpointsStat-1),:)));

% Wage-EU and wealth-EU relationships

euSample = (eStatus(:,1)==1)&(eStatus(:,13)==1|uStatus(:,13)==1);
[weuCoeff,weuInt] = regress(uStatus(euSample,13),...
    [ones(sum(euSample,1),1),log(w((eBetas(eStateInd(euSample,1))-1)*numA+eA(eStateInd(euSample,1))))']);          
[wealtheuCoeff,wealtheuInt] = regress(uStatus(euSample,13),...
    [ones(sum(euSample,1),1),zStat(zInd(euSample,1))/prodStat]);     

% Consumption shares by wealth quintile
aggc_ZQ1 = p_e*sum(sum(repmat(ind_ZQ1,1,eStates).*phi_e.*c_eStat)) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ1,1,uStates).*phi_u.*c_uStat));
aggc_ZQ2 = p_e*sum(sum(repmat(ind_ZQ2,1,eStates).*phi_e.*c_eStat)) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ2,1,uStates).*phi_u.*c_uStat));  
aggc_ZQ3 = p_e*sum(sum(repmat(ind_ZQ3,1,eStates).*phi_e.*c_eStat)) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ3,1,uStates).*phi_u.*c_uStat));  
aggc_ZQ4 = p_e*sum(sum(repmat(ind_ZQ4,1,eStates).*phi_e.*c_eStat)) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ4,1,uStates).*phi_u.*c_uStat));
aggc_ZQ5 = p_e*sum(sum(repmat(ind_ZQ5,1,eStates).*phi_e.*c_eStat)) + ...
    (1-p_e)*sum(sum(repmat(ind_ZQ5,1,uStates).*phi_u.*c_uStat));  
aggc = p_e*sum(sum(phi_e.*c_eStat)) + ...
    (1-p_e)*sum(sum(phi_u.*c_uStat));      

% Max UI relative to avg wage of UI recipients
wP_of_recipients = (1-omega0)*reshape(repmat(reshape(repmat(reshape(reshape(wP,numP,numBetas),...
    1,numP*numBetas),numT,1),1,numT*numP*numBetas),dbar,1),1,dbar*numT*numP*numBetas);
if zeta == 1
    max_rr = uimax/(sum(phi_u(:,uDur <= dbar),1)*wP_of_recipients'/sum(sum(phi_u(uDur <= dbar))));
else
    max_rr = uimax/(sum(phi_u(:,uDur <= dbar & uUI > numP),1)*wP_of_recipients'/sum(sum(phi_u(:,uDur <= dbar & uUI > numP))));
end 

%% Save results (if RCEName isn't empty)

if ~isempty(RCEName)
    if ~ident
        save(strcat(RCEName,'.mat'));
    else
        save(strcat(RCEName,'.mat'),'m','t','theta','sf','r','wealth_to_moincome',...
            'share_wealth_lt0','dwealth_ue_to_moincome',...
            'frac_unemp_UI','hhinc_ue_UI','hhinc_ue_NoUI',...
            'ur','share_ur_lt',...
            'elast_dur','conv_tightness','mowage_per_vac',...
            'gn17receipt','gn17exhaust','gn17dCdY2',...
            'mpc_QAvg','mpc_eQAvg','mpc_uQAvg','mpc_uSTQAvg',...
            'mpc_uMTQAvg','mpc_uLTQAvg','mpc_eAAvg','mpc_uAAvg',...
            'mpc_lt75k','mpc_gt75k',...
            'mpc_ZQ1','mpc_ZQ2','mpc_ZQ3','mpc_ZQ4','mpc_ZQ5',...
            'weuCoeff','weuInt','wealtheuCoeff','wealtheuInt',...
            'aggc_ZQ1','aggc_ZQ2','aggc_ZQ3','aggc_ZQ4','aggc_ZQ5',...
            'max_rr');
    end
end

toc(overallRuntime);