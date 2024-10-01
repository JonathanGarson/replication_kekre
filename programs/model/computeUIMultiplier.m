function [avgDUnemp,avgDR,UImult,UImult12,UImult_pe] = ...
    computeUIMultiplier(bSS,initName,peName,finName)
%% computeUIMultiplier.m
% Computes change in unemployment, real interest rate, and UI multiplier
% (both in GE and PE) starting from initName to finName

init = load(initName);
if ~isempty(peName)
    pe = load(peName);
end
fin = load(finName);

% Define mechanical change in taxes (not accounting for behavioral response)
numBetas = init.invariantParams.numBetas;
numA = init.invariantParams.numA;
numP = init.invariantParams.numP;
numT = init.invariantParams.numT;
numUI = init.invariantParams.numUI;
dbarbar = init.invariantParams.dbarbar;
pEpsilon = init.invariantParams.pEpsilon;
Tbar = init.varyingParams.Tbar;
zeta = init.varyingParams.zeta(:,1:Tbar);
omega0 = init.varyingParams.omega0(:,1:Tbar);
omega1 = init.varyingParams.omega1(:,1:Tbar);
omega2 = init.varyingParams.omega2(:,1:Tbar);
uStates = (dbarbar+1)*numT*numUI*numBetas;

dbar = fin.varyingParams.dbar(:,1:Tbar);
rr = fin.varyingParams.rr(:,1:Tbar);
uimax = fin.varyingParams.uimax(:,1:Tbar);

bFin = nan(uStates,Tbar);
for time=1:(Tbar-1)
    wP = pEpsilon*reshape(init.w(:,time)',numT,numP*numBetas);
    [thisB,~] = genTransfers(numBetas,numP,numT,dbarbar,dbar(time),rr(time)*(1-omega0(time)),...
        init.w(:,time)',wP,omega1(time),omega2(time),zeta(time),uimax(time));  
    bFin(:,time) = thisB';
end
bFin(:,Tbar) = bSS';
bDiff = (1-init.p_e).*nansum(init.phi_u_uStates.*(bFin-init.b),1);

% Find last month where UI is changed
endUI = find(abs(bDiff)>1e-10,1,'last'); 
% So for shock in month 1 only, returns 1,  
%        shock in month 10 only, returns 10,  
%        shock in months 1-12 only, returns 12, etc.

% Compute multipliers and other statistics of interest
initUnemp = 1-init.p_e;
finUnemp = 1-fin.p_e;
avgDUnemp = mean(finUnemp(1:endUI)-initUnemp(1:endUI));
avgDR = mean((1./fin.m(1:endUI)-1)*12 - (1./init.m(1:endUI)-1)*12);
initC = init.cTrans;
finC = fin.cTrans;
UImult = sum(finC(1:endUI) - initC(1:endUI))/sum(bDiff(1:endUI));
UImult12 = sum(finC(1:12) - initC(1:12))/sum(bDiff(1:12));
if ~isempty(peName)
    peC = pe.cTrans;
    UImult_pe = sum(peC(1:endUI) - initC(1:endUI))/sum(bDiff(1:endUI));
else
    UImult_pe = nan;
end
