function Pi = solvePiBackwards(mu,m,y)
%% solvePiBackwards.m
% Given mu, m, and y, returns inflation consistent with retailer
%  optimization under Rotemberg pricing

global epsilon psi Tbar tauR;

Pi = nan(1,Tbar);
Pi(Tbar) = 0; % assume zero inflation SS
time = Tbar-1;
while time > 0
    Pi(time) = (-1/2) + ...
        (1/2)*sqrt(1 - 4*...
        (((epsilon-1)/psi)*(1-(epsilon/(epsilon-1))*(1+tauR)*(mu(time)^(-1)))-m(time)*Pi(time+1)*(1+Pi(time+1))*y(time+1)/y(time)));    
    time = time - 1;
end
