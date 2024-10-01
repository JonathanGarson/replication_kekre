function [p_e_tilde,phi_e_tilde,phi_u_tilde,p_e,phi_e,phi_u,...
    znext_eStat,znext_uStat,sStat] = ...
    obtainInvariantSparse(z,s,znext_e,znext_u,pVec,deltaVec,...
    Tee,Teu,Tue,Tuu,zStat,eps_phi,toDisplay)
%% obtainInvariantSparse.m
% Approximates invariant distribution given policies s, znext_e, znext_u
%  over z, job-finding probabilities pVec = p(theta), separation rates in 
%  deltaVec, {eState,uState}->{eState,uState} transition 
%  matrices Tee,Teu,Tue,Tuu, and zStat for use in constructing (weakly) finer grid

tic;

% Step 1) Characterize the potentially new grid
gridpointsStat = length(zStat);

% Step 2) Linearly interpolate s, znext_e, znext_u, znext_r along zStat
znext_eStat = interp1q(z,znext_e,zStat);
znext_uStat = interp1q(z,znext_u,zStat);
sStat = interp1q(z,s,zStat);

% To avoid computational error, enforce borrowing constraint:
znext_eStat(znext_eStat < zStat(1)) = zStat(1);
znext_uStat(znext_uStat < zStat(1)) = zStat(1);

% Step 3) define job-finding probabilities and separation rates (given the 
%  potentially new grid zStat) and transition matrices for assets (given
%  the policy functions in znext_eStat,znext_uStat)

Temp = repmat(pVec,gridpointsStat,1).*sStat;
Tsep = repmat(deltaVec,gridpointsStat,1);

eStates = size(znext_e,2);
uStates = size(znext_u,2);
Tznext_e = cell(1,eStates);
Tznext_u = cell(1,uStates);

for eState=1:eStates
    Tznext_e{1,eState} = genSparseTransition(znext_eStat(:,eState),zStat);
end
for uState=1:uStates
    Tznext_u{1,uState} = genSparseTransition(znext_uStat(:,uState),zStat);
end

% Step 4) Initialize distributions (for now, uniform prior across 
%  assets, 90%/10% employed/unemployed split, and uniform prior across 
%  beta, wage, UI receipt, and duration heterogeneity)
p_e_tilde = 0.9;
phi_e_tilde = ones(gridpointsStat,eStates)/(gridpointsStat*eStates);
phi_u_tilde = ones(gridpointsStat,uStates)/(gridpointsStat*uStates);

% Define preliminaries for iteration
d = @(dist1,dist2) sqrt(nansum(nansum((dist1 - dist2).^2)));

converged = 0;
numIter = 1;

while(converged == 0 && numIter <= 50000)
    % Step 5) Update distribution
    [~,~,~,p_e_tildeNext,phi_e_tildeNext,phi_u_tildeNext] = ...
        simOnePeriod(p_e_tilde,phi_e_tilde,phi_u_tilde,...
        Temp,Tznext_e,Tznext_u,Tsep,Tee,Teu,Tue,Tuu);
    
    % Step 6) Assess convergence and update
    dist = d([p_e_tilde*phi_e_tilde,...
        (1-p_e_tilde)*phi_u_tilde],...
        [p_e_tildeNext*phi_e_tildeNext,...
        (1-p_e_tildeNext)*phi_u_tildeNext]);
    
    if dist > eps_phi
        p_e_tilde = p_e_tildeNext;
        phi_e_tilde = phi_e_tildeNext;
        phi_u_tilde = phi_u_tildeNext;
        numIter = numIter + 1;
        if toDisplay == 1
            if mod(numIter,100) == 0
                display(['Iter = ',num2str(numIter)]);
            end
        end
    else
        converged = 1;
        p_e_tilde = p_e_tildeNext;
        phi_e_tilde = phi_e_tildeNext;
        phi_u_tilde = phi_u_tildeNext;
        
        % compute final distributions in the middle of the period      
        p_e = sum(sum((p_e_tilde*phi_e_tilde + ...
            ((1-p_e_tilde)*phi_u_tilde.*Temp)*Tue)));
        phi_e = (1/p_e)*(p_e_tilde*phi_e_tilde + ...
            ((1-p_e_tilde)*phi_u_tilde.*Temp)*Tue);
        p_u = sum(sum((1-p_e_tilde)*phi_u_tilde.*(1-Temp))); 
         % = 1-p_e, but use this next to avoid miniscule computational error
        phi_u = (1/p_u)*((1-p_e_tilde)*phi_u_tilde.*(1-Temp));
    end
end     

timeTaken = toc;

% Display performance
if converged == 1
    if toDisplay == 1
        disp(['Stationary dist obtained: ',num2str(timeTaken),...
            's, ',num2str(numIter), ' iterations, ',...
            num2str(timeTaken/numIter), ' s/iteration.  Distance = ',...
            num2str(dist)]);
    end
else
    disp('Stationary dist unable to converge in 50000 iterations');
    pause;
end