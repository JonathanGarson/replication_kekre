function [z,v_e,v_u,s,c_e,znext_e,c_u,znext_u] = ...
    statPoliciesEGM(m,pVec,deltaVec,Tee,Teu,Tue,Tuu,y_e,y_u,...
    betas_e,betas_u,eps_v,toDisplay,znext,v_e,v_u)
%% statPoliciesEGM.m
% Given m, p = p(theta), separation rates, {eState,uState}->
%  {eState,uState} transition matrices, y_e, y_u, discount factors 
%  for e, u agents, gridpoints znext, as well as OPTIONAL initial 
%  v_e, v_u, uses endogenous-gridpoint-method-augmented value function
%  iteration to compute optimal policies and value functions over the grid

tic;
global xi sigma zbar;

gridPoints = length(znext);

% Define whether initial v's exist
withStart = exist('v_e','var');

% Define useful functions
largeNeg = -999999999;
if sigma ~= 1
    u = @(c) (c>=0).*((c.^(1-sigma)-1)./(1-sigma))+(c<0).*(largeNeg);
 % setting utility to negative large negative number for negative 
 % consumption ensures that in iteration below, a policy inducing this
 % outcome is never chosen
else
    u = @(c) (c>=0).*(log(c))+(c<0).*(largeNeg);
end
psi = @(eval) (eval).^(xi+1);
psiPrimeInv = @(val) (val./(xi+1)).^(1/xi);
d = @(v1,v2) max(max(abs(v1-v2)));

% Define size of employed, unemployed value+policy matrices
eStates = length(y_e);
uStates = length(y_u);

% Step 1) guess initial value functions if not already provided.  For now 
%  I assume that value functions are identical across worker types.
if withStart == 0
    v_e = repmat(interp1q(znext(1)+(logspace(0,1,gridPoints)-1)'/9*(znext(gridPoints)-znext(1)),...
        log(linspace(10,100,gridPoints))',znext),1,eStates);
    v_u = zeros(gridPoints,uStates);
end

converged = 0;
numIter = 1;
while(converged == 0 && numIter <= 50000)
    % Step 3) construct updated value and policy functions   
    if numIter == 1
        [~,z_e_endog,c_e_endog,vupd_e_endog,gothicv_e,...
            z_u_endog,c_u_endog,vupd_u_endog,gothicv_u] = ...
            oneWorkerEGM(znext,v_e,v_u,-99,-99,...
            m,pVec,deltaVec,Tee,Teu,Tue,Tuu,y_e,y_u,betas_e,betas_u);
    else
        [~,z_e_endog,c_e_endog,vupd_e_endog,gothicv_e,...
            z_u_endog,c_u_endog,vupd_u_endog,gothicv_u] = ...
            oneWorkerEGM(znext,v_e,v_u,c_e,c_u,...
            m,pVec,deltaVec,Tee,Teu,Tue,Tuu,y_e,y_u,betas_e,betas_u);        
    end
    
    % Interpolate new value and consumption policy functions along z = znext
    %  If artificial upper bound on assets is optimally chosen at 
    %  z < max(z), so that for values of z in between, an even higher 
    %  amount of assets is chosen -- that is set to nan using interp1q
    %  (though won't happen in stationary RCE)
  
    c_e = nan(gridPoints,eStates);
    vupd_e = nan(gridPoints,eStates);
    c_u = nan(gridPoints,uStates);
    vupd_u = nan(gridPoints,uStates);
    
    for eState=1:eStates      
        zToAdd = znext(z_e_endog(1,eState) > znext);
        % if non-empty, this means that the lower bound on assets is 
        %  optimally chosen at z > zbar, so that for all values of z 
        %  in between, I know that zbar is chosen            
        z_e_endog_tointerp = [zToAdd;z_e_endog(:,eState)];
        c_e_endog_tointerp = [y_e(eState) - m*zbar + zToAdd;...
            c_e_endog(:,eState)];
        vupd_e_endog_tointerp = [u(y_e(eState) - m*zbar + zToAdd) + ...
            betas_e(eState)*gothicv_e(1,eState);vupd_e_endog(:,eState)];
        % now vupd_e_endog_tointerp will preserve non-linearity in region 
        %  where agents are constrained
        
        c_e(:,eState) = interp1q(z_e_endog_tointerp,...
            c_e_endog_tointerp,znext);
        vupd_e(:,eState) = interp1q(z_e_endog_tointerp,...
            vupd_e_endog_tointerp,znext);
    end

    for uState=1:uStates
        zToAdd = znext(z_u_endog(1,uState) > znext);
        % if non-empty, this means that the lower bound on assets is 
        %  optimally chosen at z > zbar, so that for all values of z 
        %  in between, I know that zbar is chosen 
        z_u_endog_tointerp = [zToAdd;z_u_endog(:,uState)];
        c_u_endog_tointerp = [y_u(uState) - m*zbar + zToAdd;...
            c_u_endog(:,uState)];
        vupd_u_endog_tointerp = [u(y_u(uState) - m*zbar + zToAdd) + ...
            betas_u(uState)*gothicv_u(1,uState);vupd_u_endog(:,uState)];
        % now vupd_u_endog_tointerp will preserve non-linearity in region 
        %  where agents are constrained   

        c_u(:,uState) = interp1q(z_u_endog_tointerp,...
            c_u_endog_tointerp,znext);
        vupd_u(:,uState) = interp1q(z_u_endog_tointerp,...
            vupd_u_endog_tointerp,znext);
    end
    
    % Step 4) compute distance and iterate
    dist = d([v_e,v_u],[vupd_e,vupd_u]);
    if dist > eps_v
        v_e = vupd_e;
        v_u = vupd_u;
        numIter = numIter + 1;
        if toDisplay == 1
            if mod(numIter,100) == 0
                display(['Iter = ',num2str(numIter),' ',num2str(dist)]);
            end
        end
    else
        converged = 1;
        % Set final value functions defined over an INHERITED grid z
        % For simplicity, I take z = znext
        z = znext;
        v_e = vupd_e;
        v_u = vupd_u;
    end
end

% Compute ex-post policy functions
% - for next period's assets, znext = z makes this easy
znext_e = (1/m)*(repmat(y_e,gridPoints,1) + repmat(znext,1,eStates) - c_e);
znext_u = (1/m)*(repmat(y_u,gridPoints,1) + repmat(znext,1,uStates) - c_u);
% - and for consumption, c_e, c_u are already defined over znext = z

% And lastly compute ex-ante policy function for labor supply
ve_bigger = (v_e*Tue'-v_u) >= 0;
sTemp = psiPrimeInv(repmat(pVec,gridPoints,1).*(v_e*Tue' - v_u)).*ve_bigger;
s = min(repmat(1./pVec,gridPoints,1),sTemp); 

timeTaken = toc;

% Display performance if toDisplay == 1
if converged == 1
    if toDisplay == 1
        disp(['Iterative algorithm complete: ',num2str(timeTaken),...
            's, ',num2str(numIter), ' iterations, ',...
            num2str(timeTaken/numIter), ' s/iteration.  Distance = ',...
            num2str(dist)]);
    end
else
    disp('Iterative algorithm unable to converge in 50000 iterations');
    pause;
end