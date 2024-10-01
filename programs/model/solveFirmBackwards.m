function [sfTrans,mu] = solveFirmBackwards(Tee,deltaMat,m,w,mu,sfSS)
%% solveFirmBackwards.m
% Solves firm surplus from each worker type backwards

global Tbar numA a;

eStates = length(sfSS);
sfTrans = nan(length(sfSS),Tbar);
sfTrans(:,Tbar) = sfSS';
for time=Tbar-1:-1:1
    sfTrans(:,time) = (mu(time)^(-1))*repmat(a(:,time),eStates/numA,1) - w(:,time) + (1-deltaMat(:,time)).*m(time).*(Tee{time}*sfTrans(:,time+1));
end