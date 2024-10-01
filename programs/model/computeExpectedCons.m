function [expectc_e,expectc_u] = ...
    computeExpectedCons(z,s,znext_e,znext_u,pVec,deltaVec,...
    Tee,Teu,Tue,Tuu,c_e,c_u)
%% computeExpectedCons.m
% Computes the expected values of consumption (policies c_e, c_u)
%  given state variables today

gridpoints = length(z);

% Define job-finding probabilities and separation rates and 
%  transition matrices for assets (given the policy functions in 
%  znext_e,znext_u)

Temp = repmat(pVec,gridpoints,1).*s;
Tsep = repmat(deltaVec,gridpoints,1);

eStates = size(znext_e,2);
uStates = size(znext_u,2);
Tznext_e = cell(1,eStates);
Tznext_u = cell(1,uStates);

for eState=1:eStates
    Tznext_e{1,eState} = genSparseTransition(znext_e(:,eState),z);
end
for uState=1:uStates
    Tznext_u{1,uState} = genSparseTransition(znext_u(:,uState),z);
end

% Expected next period consumption for agents as a function of next 
%  period's assets but this period's other state variables
expectc_e_nextz = (1-Tsep).*(c_e*Tee') + ...
    Tsep.*((Temp.*(c_e*Tue') + (1-Temp).*c_u)*Teu');
expectc_u_nextz = (Temp.*(c_e*Tue') + (1-Temp).*c_u)*Tuu'; 

% Expected next period consumption as a function of this period's assets
expectc_e = nan(gridpoints,eStates);
for eState=1:eStates
    expectc_e(:,eState) = Tznext_e{1,eState}*expectc_e_nextz(:,eState);
end
expectc_u = nan(gridpoints,uStates);
for uState=1:uStates
    expectc_u(:,uState) = Tznext_u{1,uState}*expectc_u_nextz(:,uState);
end