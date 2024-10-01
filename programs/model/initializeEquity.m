function [equityShares] = initializeEquity(startQuant,valInit,z,phiInit)
%% initializeEquity.m
% Initializes equity shares using piecewise log-linear specification,
%  assuming that wealth is the only state variable relevant

gridPoints = length(z);
rawEquityShares = nan(gridPoints,1);
cdf_z = cumsum(nansum(phiInit,2));

% Zero shares below startQuant
indStart = find(cdf_z > startQuant,1); 
rawEquityShares(1:(indStart-1)) = 0;

% Log-linear shares after startQuant (placeholder slope of 0.05, normalized below anyway)
rawEquityShares(indStart:gridPoints) = 0.05*(log(z(indStart:gridPoints)-log(z(indStart))));

% Zero share in final 1% of gridpoints
indAdjust = floor(0.99*gridPoints);
rawEquityShares(indAdjust:gridPoints) = 0;

% Scale shares to be consistent with aggregate value of equity 
equityAggInit = sum(z.*rawEquityShares.*nansum(phiInit,2));
scaleFactor = valInit/equityAggInit;
equityShares = scaleFactor*rawEquityShares;