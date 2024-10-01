function [p_e,phi_e,phi_u,p_e_tildeNext,phi_e_tildeNext,phi_u_tildeNext] = ...
    simOnePeriod(p_e_tilde,phi_e_tilde,phi_u_tilde,...
    Temp,Tznext_e,Tznext_u,Tsep,Tee,Teu,Tue,Tuu)
%% simOnePeriod.m
% Calculates the distributions of agents in the middle of the period, and
%  the beginning of next period, given job-finding
%  probabilities in Temp, separation rates in Tsep, and transition matrices 
%  Tznext_e,Tznext_u,Tee,Teu,Tue,Tuu
% Note that Temp is gridpoints x uStates
%           Tznext_e is gridpoints x gridpoints (sparse)
%           Tznext_u is gridpoints x gridpoints (sparse)
%           Tsep is gridpoints x eStates
%           Tee is eStates x eStates (sparse)
%           Teu is eStates x uStates (sparse)
%           Tue is uStates x eStates (sparse)
%           Tuu is uStates x uStates (sparse)

gridpoints = size(phi_e_tilde,1);
eStates = size(phi_e_tilde,2);
uStates = size(phi_u_tilde,2);

% Distribution of agents in the middle of the period
p_e = sum(sum(p_e_tilde*phi_e_tilde + ...
    (((1-p_e_tilde)*phi_u_tilde).*Temp)*Tue));
phi_e = (1/p_e)*(p_e_tilde*phi_e_tilde + ...
    (((1-p_e_tilde)*phi_u_tilde).*Temp)*Tue);
p_u = sum(sum(((1-p_e_tilde)*phi_u_tilde).*(1-Temp))); 
 % = 1-p_e, but use this next to avoid miniscule computational error
phi_u = (1/p_u)*(((1-p_e_tilde)*phi_u_tilde).*(1-Temp));

% Looking ahead, distributions of agents by *next* period's debt but this
%  period's other states
phi_e_znext = zeros(gridpoints,eStates);
phi_u_znext = zeros(gridpoints,uStates);
for eState=1:eStates
    phi_e_znext(:,eState) = (phi_e(:,eState)'*Tznext_e{eState})';
end
for uState=1:uStates
    phi_u_znext(:,uState) = (phi_u(:,uState)'*Tznext_u{uState})';
end

% Next period's distributions of agents 
p_e_tildeNext = sum(sum((p_e*phi_e_znext.*(1-Tsep))*Tee));
phi_e_tildeNext = (1/p_e_tildeNext)*((p_e*phi_e_znext.*(1-Tsep))*Tee);

p_u_tildeNext = sum(sum((p_e*phi_e_znext.*Tsep)*Teu + ...
    (1-p_e)*phi_u_znext*Tuu));
 % = 1 - p_e_tildeNext, but use this next to avoid miniscule computational error
phi_u_tildeNext = (1/p_u_tildeNext)*...
    ((p_e*phi_e_znext.*Tsep)*Teu + ...
    (1-p_e)*phi_u_znext*Tuu);