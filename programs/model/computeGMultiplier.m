function [avgDUnemp,avgDR,Gmult,Gmult12,Gmult_pe] = ...
    computeGMultiplier(initName,peName,finName)
%% computeGMultiplier.m
% Computes change in unemployment, real interest rate, and govt spending 
%  multiplier (both in GE and PE) starting from initName to finName

init = load(initName);
if ~isempty(peName)
    pe = load(peName);
end
fin = load(finName);

gDiff = fin.varyingParams.g - init.varyingParams.g;

% Find last month where g is changed
endg = find(abs(gDiff)>1e-10,1,'last'); 

% Compute multipliers and other statistics of interest
initUnemp = 1-init.p_e;
finUnemp = 1-fin.p_e;
avgDUnemp = mean(finUnemp(1:endg)-initUnemp(1:endg));
avgDR = mean((1./fin.m(1:endg)-1)*12 - (1./init.m(1:endg)-1)*12);
initC = init.cTrans;
finC = fin.cTrans;
Gmult = sum(finC(1:endg)+gDiff(1:endg) - initC(1:endg))/sum(gDiff(1:endg));
Gmult12 = sum(finC(1:12)+gDiff(1:12) - initC(1:12))/sum(gDiff(1:12));
if ~isempty(peName)
    peC = pe.cTrans;
    Gmult_pe = sum(peC(1:endg)+gDiff(1:endg) - initC(1:endg))/sum(gDiff(1:endg));
else
    Gmult_pe = nan;
end