function [phiFin] = revalueEquity(valInit,valFin,equityShares,z,phiInit,zLastPosDens)
%% revalueEquity.m
% Revalues asset positions given an unexpected shock to firm equity from
%  valInit to valFin, overall positions distributed by phiInit over 
%  assets z, and equity portfolio shares distributed by equityShares over 
%  assets.  Returns the distribution over revalued asset positions.

states = size(equityShares,2); % number of other state variables besides assets

% Compute implied initial aggregate wealth in equity
equityAggInit = sum(sum(repmat(z,1,states).*equityShares.*phiInit));

% Scale equity shares (uniformly) to be consistent with actual initial
%  aggregate wealth in equity in my stationary RCE
scaleFactor = valInit/equityAggInit;
scaledEquityShares = scaleFactor*equityShares;

% Compute implied wealth in bonds vs. equity (which implicitly makes most
%  sense if equityShares>0 <=> z>0)
bondWealth = repmat(z,1,states).*(1-scaledEquityShares);
equityWealthInit = repmat(z,1,states).*scaledEquityShares;

% Compute revalued equity wealth
equityWealthFin = (valFin/valInit)*equityWealthInit;

% Compute revalued overall assets
revaluedZ = bondWealth + equityWealthFin;

% Check to make sure no agents (with positive density) have revalued wealth 
%  below the assumed lower bound of z or above the assumed upper bound of z
if min(min(revaluedZ(1:zLastPosDens,:))) < z(1) || ...
        max(max(revaluedZ(1:zLastPosDens,:))) > z(length(z))
    error('Revaluation pushes density of agents outside assumed bounds!');
else
    outOfBounds = revaluedZ(:,:) < z(1) | revaluedZ(:,:) > z(length(z));
    revaluedZ(outOfBounds) = nan; % ensures no issues with constructing sparse transition below 
end

% Compute distribution over revalued assets
phiFin = zeros(length(z),states);
for state=1:states
    Trevalue = genSparseTransition(revaluedZ(:,state),z);
    phiFin(:,state) = (phiInit(:,state)'*Trevalue)';
end