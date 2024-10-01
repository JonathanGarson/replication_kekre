function [zInd,eTildeStatus,eStateTildeInd,uTildeStatus,uStateTildeInd,...
    eStatus,eStateInd,uStatus,uStateInd] = ...
    simStationaryAgents(p_e_tilde,phi_e_tilde,phi_u_tilde,z,...
    znext_e,znext_u,s,pVec,deltaVec,Tee,Teu,Tue,Tuu,...
    numAgents,T,toDisplay,sepOnly)  
%% simStationaryAgents.m
% Simulates numAgents sampled from the stationary RCE for T periods

timeTaken = tic;

%% Initialize objects holding returned simulation paths of interest
zInd = nan(numAgents,T);
eTildeStatus = logical(zeros(numAgents,T));
eStateTildeInd = ones(numAgents,T); % if agent isn't e at a given date,
% this will contain dummy val of 1.  But this means that 
% indexing with this matrix should also be accompanied by eTildeStatus==1
% so that dummy values aren't counted!
uTildeStatus = logical(zeros(numAgents,T));
uStateTildeInd = ones(numAgents,T); % if agent isn't u at a given date,
% this will contain dummy val of 1.  But this means that 
% indexing with this matrix should also be accompanied by uTildeStatus==1
% so that dummy values aren't counted!
eStatus = logical(zeros(numAgents,T));
eStateInd = ones(numAgents,T); % if agent isn't e at a given date,
% this will contain dummy val of 1.  But this means that 
% indexing with this matrix should also be accompanied by eStatus==1
% so that dummy values aren't counted!
uStatus = logical(zeros(numAgents,T));
uStateInd = ones(numAgents,T); % if agent isn't u at a given date,
% this will contain dummy val of 1.  But this means that 
% indexing with this matrix should also be accompanied by uStatus==1
% so that dummy values aren't counted!

%% Initialize agents' state variables at beginning of date 1
%  (where agents are ordered (WLOG) e,u at date 1)

initEtilde = round(numAgents*p_e_tilde);
initUtilde = numAgents - initEtilde;

phi_e_tildeCDF = [0;cumsum(phi_e_tilde(:))]; 
phi_u_tildeCDF = [0;cumsum(phi_u_tilde(:))];
% need to add zeros so that I correctly assign the right fraction of agents
%  at the borrowing constraint below

eTildeStatus(1:initEtilde,1) = 1;
uTildeStatus((initEtilde+1):numAgents,1) = 1;

randNums = rand(numAgents,1);
gridPoints = length(z);

[~,~,initEStateTilde] = histcounts(randNums(1:initEtilde),...
    phi_e_tildeCDF);
zInd(1:initEtilde,1) = mod(initEStateTilde-1,gridPoints)+1;
eStateTildeInd(1:initEtilde,1) = ceil(initEStateTilde/gridPoints);

[~,~,initUStateTilde] = histcounts(randNums((initEtilde+1):numAgents),...
    phi_u_tildeCDF);
zInd((initEtilde+1):numAgents,1) = mod(initUStateTilde-1,gridPoints)+1;
uStateTildeInd((initEtilde+1):numAgents,1) = ceil(initUStateTilde/gridPoints);

%% Create matrices summarizing z transitions, similar to algorithm in 
%   genSparseTransition invoked in obtainInvariantSparse.m

eStates = size(phi_e_tilde,2);
uStates = size(phi_u_tilde,2);

znext_e_bin = nan(gridPoints,eStates);
znext_e_IndLo = nan(gridPoints,eStates);
znext_e_IndHi = nan(gridPoints,eStates);
prob_e_Lo = nan(gridPoints,eStates);
for eState=1:eStates
    [~,~,znext_e_bin(:,eState)] = histcounts(znext_e(:,eState),z);
    znext_e_IndLo(:,eState) = znext_e_bin(:,eState);
    znext_e_IndHi(:,eState) = znext_e_bin(:,eState)+1;
    prob_e_Lo(:,eState) = (z(znext_e_IndHi(:,eState)) - znext_e(:,eState))./...
        (z(znext_e_IndHi(:,eState))-z(znext_e_IndLo(:,eState)));
end

znext_u_bin = nan(gridPoints,uStates);
znext_u_IndLo = nan(gridPoints,uStates);
znext_u_IndHi = nan(gridPoints,uStates);
prob_u_Lo = nan(gridPoints,uStates);
for uState=1:uStates
    [~,~,znext_u_bin(:,uState)] = histcounts(znext_u(:,uState),z);
    znext_u_IndLo(:,uState) = znext_u_bin(:,uState);
    znext_u_IndHi(:,uState) = znext_u_bin(:,uState)+1;
    prob_u_Lo(:,uState) = (z(znext_u_IndHi(:,uState)) - znext_u(:,uState))./...
        (z(znext_u_IndHi(:,uState))-z(znext_u_IndLo(:,uState)));
end

%% Create matrices summarizing {eState,uState}->{eState,uState} transitions

Temp = repmat(pVec,gridPoints,1).*s;
if sepOnly == 1
    Tsep1 = repmat(deltaVec/max(deltaVec),gridPoints,1);
    % holds the probability of separation at the end of period 1 if
    %  sepOnly == 1.  implies that relative probability of separation 
    %  across types of employed workers is unchanged from deltaVec, but 
    %  overall level of separation has been scaled up
    TsepBase = repmat(deltaVec,gridPoints,1);
else
    Tsep1 = repmat(deltaVec,gridPoints,1);
    TsepBase = repmat(deltaVec,gridPoints,1);
end

%% Simulate

% Generate random numbers used in simulation
randEmp = rand(numAgents,T);
randAssets = rand(numAgents,T);
randSep = rand(numAgents,T);

randEE = rand(numAgents,T);
randEU = rand(numAgents,T);

randUE = rand(numAgents,T);
randUU = rand(numAgents,T);

getNextState = @(probs,rands) sum(repmat(rands,1,size(probs,2)) > cumsum(probs,2),2)+1;
 % for a matrix probs where each row constains the PDF of transitioning 
 %  from that rowState to each columnState, as well as a vector of random number
 %  draws for each rowState, returns the randomly chosen columnState for
 %  each rowState

% Simulate
for t=1:T
    if t==1
        Tsep = Tsep1;
    else
        Tsep = TsepBase;
    end
    
    % Simulate employment status in the middle of the period  
    %  - U --> E
    ue = uTildeStatus(:,t)==1 & ...
        randEmp(:,t)<Temp((uStateTildeInd(:,t)-1)*gridPoints + zInd(:,t));
    eStatus(ue,t) = 1;
    eStateInd(ue,t) = getNextState(Tue(uStateTildeInd(ue,t),:),randUE(ue,t));
       
    %  - U --> U (simple)
    uu = uTildeStatus(:,t)==1 & ...
        randEmp(:,t)>=Temp((uStateTildeInd(:,t)-1)*gridPoints + zInd(:,t));
    uStatus(uu,t) = 1;
    uStateInd(uu,t) = uStateTildeInd(uu,t);
    
    %  - E --> E (simple)  
    ee = eTildeStatus(:,t)==1;
    eStatus(ee,t) = 1;
    eStateInd(ee,t) = eStateTildeInd(ee,t);    
    
    if t < T
        % Simulate next period's assets based on this period's
        %  consumption-savings problem
        zInd(eStatus(:,t)==1,t+1) = ...
            znext_e_IndLo((eStateInd(eStatus(:,t)==1,t)-1)*gridPoints + zInd(eStatus(:,t)==1,t)).*...
            (randAssets(eStatus(:,t)==1,t) < ...
            prob_e_Lo((eStateInd(eStatus(:,t)==1,t)-1)*gridPoints + zInd(eStatus(:,t)==1,t))) + ...
            znext_e_IndHi((eStateInd(eStatus(:,t)==1,t)-1)*gridPoints + zInd(eStatus(:,t)==1,t)).*...
            (randAssets(eStatus(:,t)==1,t) >= ...
            prob_e_Lo((eStateInd(eStatus(:,t)==1,t)-1)*gridPoints + zInd(eStatus(:,t)==1,t)));
        zInd(uStatus(:,t)==1,t+1) = ...
            znext_u_IndLo((uStateInd(uStatus(:,t)==1,t)-1)*gridPoints + zInd(uStatus(:,t)==1,t)).*...
            (randAssets(uStatus(:,t)==1,t) < ...
            prob_u_Lo((uStateInd(uStatus(:,t)==1,t)-1)*gridPoints + zInd(uStatus(:,t)==1,t))) + ...
            znext_u_IndHi((uStateInd(uStatus(:,t)==1,t)-1)*gridPoints + zInd(uStatus(:,t)==1,t)).*...
            (randAssets(uStatus(:,t)==1,t) >= ...
            prob_u_Lo((uStateInd(uStatus(:,t)==1,t)-1)*gridPoints + zInd(uStatus(:,t)==1,t))); 

        % Simulate flows into employment at the start of next period
        % - E --> E
        ee = eStatus(:,t)==1 & ...
            randSep(:,t)>Tsep((eStateInd(:,t)-1)*gridPoints + zInd(:,t));
        eTildeStatus(ee,t+1) = 1;
        eStateTildeInd(ee,t+1) = getNextState(Tee(eStateInd(ee,t),:),randEE(ee,t));

        % Simulate flows into unemployment at the start of next period
        % - E --> U
        eu = eStatus(:,t)==1 & ...
            randSep(:,t)<=Tsep((eStateInd(:,t)-1)*gridPoints + zInd(:,t));
        uTildeStatus(eu,t+1) = 1;
        uStateTildeInd(eu,t+1) = getNextState(Teu(eStateInd(eu,t),:),randEU(eu,t));

        % - U --> U    
        uu = uStatus(:,t)==1;
        uTildeStatus(uu,t+1) = 1;
        uStateTildeInd(uu,t+1) = getNextState(Tuu(uStateInd(uu,t),:),randUU(uu,t));
    end
end

%% Plot (if toDisplay == 1)

if toDisplay == 1
    display(['Consumption simulation complete: ', ...
        num2str(t), ' periods, ', num2str(toc(timeTaken)), 's, ', ...
        num2str(timeTaken/t), 's/period.']);
    
    simLambdaE = sum(eStatus,1)./(sum(eStatus,1)+sum(uStatus,1));
    display('Unemployment rate:');
    disp([1-simLambdaE(1),1-simLambdaE(T)]);
    % obviously shouldn't expect this to be stationary when sepOnly = 1
    
    figure;
    h1 = subplot(1,2,1);
    histogram(h1,z(zInd(eStatus(:,1)==1,1)),20,'Normalization','probability');
    title('Simulated distribution of employed agents by assets, t=1');
    h2 = subplot(1,2,2);
    histogram(h2,z(zInd(eStatus(:,T)==1,T)),20,'Normalization','probability');
    title('Simulated distribution of employed agents by assets, t=T');

    figure;
    h3 = subplot(1,2,1);
    histogram(h3,z(zInd(uStatus(:,1)==1,1)),20,'Normalization','probability');
    title('Simulated distribution of unemployed agents by assets, t=1');
    h4 = subplot(1,2,2);
    histogram(h4,z(zInd(uStatus(:,T)==1,T)),20,'Normalization','probability');
    title('Simulated distribution of unemployed agents by assets, t=T');    
end