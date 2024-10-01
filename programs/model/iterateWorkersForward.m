function [p_e_tilde,phi_e_tilde,phi_u_tilde,p_e,phi_e,phi_u] = ...
    iterateWorkersForward(z,s,znext_e,znext_u,pMat,...
    deltaMat,Tee,Teu,Tue,Tuu,p_e_tilde1,phi_e_tilde1,phi_u_tilde1,toDisplay,tabs)
%% iterateWorkersForward.m
% Iterates workers forward starting from
%  {p_e_tilde1,phi_e_tilde1,phi_u_tilde1} given policies in 
%  s, znext_e, znext_u over inherited z gridpoints in z, and with
%  job-finding rates p(theta), separation rates Tsep, 
%  and {e,u}->{e,u} transition matrices Tee,Teu,Tue,Tuu

tic;
global Tbar; 

% Useful constants
gridPoints = length(z);
eStates = size(znext_e,2);
uStates = size(znext_u,2);

% Define job-finding probabilities and separation rates
Temp = repmat(reshape(pMat,[1,uStates,Tbar]),gridPoints,1,1).*s;
Temp(isnan(Temp)) = 0; 
Tsep = repmat(reshape(deltaMat,[1,eStates,Tbar]),gridPoints,1,1);

% Initialize distributions and other preliminaries for iteration
p_e_tilde = nan(1,Tbar);
phi_e_tilde = zeros(gridPoints,eStates,Tbar);
phi_u_tilde = zeros(gridPoints,uStates,Tbar);
p_e = nan(1,Tbar);
phi_e = zeros(gridPoints,eStates,Tbar);
phi_u = zeros(gridPoints,uStates,Tbar);

p_e_tilde(1) = p_e_tilde1;
phi_e_tilde(:,:,1) = phi_e_tilde1;
phi_u_tilde(:,:,1) = phi_u_tilde1;

Tznext_e = cell(1,eStates);
Tznext_u = cell(1,uStates);

t=1;
while t <= Tbar
    % Transform period t's asset policy functions into transition matrices
    for eState=1:eStates
        Tznext_e{1,eState} = genSparseTransition(znext_e(:,eState,t),z);
    end
    for uState=1:uStates
        Tznext_u{1,uState} = genSparseTransition(znext_u(:,uState,t),z);
    end
    
    % Update distribution
    [p_eThis,phi_eThis,phi_uThis,p_e_tildeNext,phi_e_tildeNext,phi_u_tildeNext] = ...
        simOnePeriod(p_e_tilde(t),phi_e_tilde(:,:,t),phi_u_tilde(:,:,t),...
        Temp(:,:,t),Tznext_e,Tznext_u,Tsep(:,:,t),Tee{t},Teu{t},Tue,Tuu); 
    
    p_e(t) = p_eThis;
    phi_e(:,:,t) = phi_eThis;
    phi_u(:,:,t) = phi_uThis;
    if t < Tbar
        p_e_tilde(t+1) = p_e_tildeNext;
        phi_e_tilde(:,:,t+1) = phi_e_tildeNext;
        phi_u_tilde(:,:,t+1) = phi_u_tildeNext;
    end
    t = t+1;
end

timeTaken = toc;
d = @(dist1,dist2) sqrt(nansum(nansum((dist1 - dist2).^2)));
% Display performance if toDisplay = 'on'
if strcmp(toDisplay,'on')
    disp([tabs,'Forward worker iteration complete: ',num2str(timeTaken),...
        's, ',num2str(Tbar), ' iterations, ',...
        num2str(timeTaken/Tbar), ' s/iteration.  d(phi1,phiTbar) = ',...
        num2str(d([p_e_tilde(1)*phi_e_tilde(:,:,1),...
        (1-p_e_tilde(1))*phi_u_tilde(:,:,1)],...
        [p_e_tilde(Tbar)*phi_e_tilde(:,:,Tbar),...
        (1-p_e_tilde(Tbar))*phi_u_tilde(:,:,Tbar)]))]);
end