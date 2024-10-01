function [v_e,v_u,s,c_e,znext_e,c_u,znext_u,zLastPol] = ...
    iteratePoliciesBackward(z,...
    v_eTbar,c_eTbar,znext_eTbar,v_uTbar,c_uTbar,znext_uTbar,sTbar,...
    m,pMat,deltaMat,Tee,Teu,Tue,Tuu,y_e,y_u,...
    zbarPrev,zbar,betas_e,betas_u,zLastPosDens,origGridInds,toDisplay,tabs)
%% iteratePoliciesBackward.m
% Iterates policies and value functions backward Tbar periods

tic;
global xi sigma Tbar;

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

% Initialize key variables
gridPoints = length(z);
eStates = size(v_eTbar,2);
uStates = size(v_uTbar,2);

v_e = nan(gridPoints,eStates,Tbar);
v_u = nan(gridPoints,uStates,Tbar);
s = nan(gridPoints,uStates,Tbar);
c_e = nan(gridPoints,eStates,Tbar);
znext_e = nan(gridPoints,eStates,Tbar);
c_u = nan(gridPoints,uStates,Tbar);
znext_u = nan(gridPoints,uStates,Tbar);

v_e(:,:,Tbar) = v_eTbar;
v_u(:,:,Tbar) = v_uTbar;
s(:,:,Tbar) = sTbar;
c_e(:,:,Tbar) = c_eTbar;
znext_e(:,:,Tbar) = znext_eTbar;
c_u(:,:,Tbar) = c_uTbar;
znext_u(:,:,Tbar) = znext_uTbar;

% EGM-augmented VFI backwards
zLastIter = z;
v_eLastIter = v_eTbar;
v_uLastIter = v_uTbar;
c_eLastIter = c_eTbar;
c_uLastIter = c_uTbar;

time = Tbar-1;
while time > 0   
    [~,z_e_endog,c_e_endog,v_e_endog,gothicv_e,...
        z_u_endog,c_u_endog,v_u_endog,gothicv_u] = ...
        oneWorkerEGM(zLastIter,v_eLastIter,v_uLastIter,...
        c_eLastIter,c_uLastIter,m(time),...
        pMat(:,time+1)',deltaMat(:,time)',Tee{time},Teu{time},Tue,Tuu,...
        y_e(:,time)',y_u(:,time)',betas_e(:,time)',betas_u(:,time)');     
   
    % Define relevant portion of z grid at date 'time' (including
    % all gridpoints up to zmax)
    zmax = min([max(z_e_endog,[],1),max(z_u_endog,[],1)]);
    zThisIterInd = (zmax >= z); 
    zThisIter = z(zThisIterInd); 
    
    % Add zbarPrev(time) if not equal to zmin (assume that it is
    %  impossible that zbarPrev(time) will be exactly equal to one of the gridpoints)
    zbarPrevInd = find(zbarPrev(time) <= z,1);
    if zbarPrevInd > 1
        zThisIter = [zThisIter(1:zbarPrevInd-1);zbarPrev(time);zThisIter(zbarPrevInd:end)];
        zThisIterIndMapping = logical([zThisIterInd(1:zbarPrevInd-1);0;zThisIterInd(zbarPrevInd:end)]);
    else
        zThisIterIndMapping = zThisIterInd;
    end
        
    % Interpolate new value and consumption policy functions along zThisIter
    
    this_c_e = nan(length(zThisIter),eStates);
    this_v_e = nan(length(zThisIter),eStates);
    for eState=1:eStates
        zToAdd = z(z_e_endog(1,eState) > z);
        % if non-empty, this means that the lower bound on assets is 
        %  optimally chosen at z > min(z), so that for all values of z 
        %  in between, I know that zbar(time) is chosen
        % note that values of z between [min(z),zbarPrev(time))
        %  technically should not be reached in simulation, but because 
        %  simulation doesn't feature all values of zbar exactly on grid,
        %  they can be
        z_e_endog_tointerp = [zToAdd;z_e_endog(:,eState)];
        c_e_endog_tointerp = [y_e(eState,time) - m(time)*zbar(time) + zToAdd;...
            c_e_endog(:,eState)];
        v_e_endog_tointerp = [u(y_e(eState,time) - m(time)*zbar(time) + zToAdd) + ...
            betas_e(eState,time)*gothicv_e(1,eState);v_e_endog(:,eState)];
        % now v_e_endog_tointerp will preserve non-linearity in region 
        %  where agents are constrained
        
        this_c_e(:,eState) = interp1q(z_e_endog_tointerp,c_e_endog_tointerp,zThisIter);
        this_v_e(:,eState) = interp1q(z_e_endog_tointerp,v_e_endog_tointerp,zThisIter);
        c_e(zThisIterInd,eState,time) = this_c_e(zThisIterIndMapping,eState);
        v_e(zThisIterInd,eState,time) = this_v_e(zThisIterIndMapping,eState);     
    end

    this_c_u = nan(length(zThisIter),uStates);
    this_v_u = nan(length(zThisIter),uStates);    
    for uState=1:uStates
        zToAdd = z(z_u_endog(1,uState) > z);
        % if non-empty, this means that the lower bound on assets is 
        %  optimally chosen at z > zbarPrev(time), so that for all values of z 
        %  in between, I know that zbar(time) is chosen
        
        z_u_endog_tointerp = [zToAdd;z_u_endog(:,uState)];
        c_u_endog_tointerp = [y_u(uState,time) - m(time)*zbar(time) + zToAdd;...
            c_u_endog(:,uState)];
        v_u_endog_tointerp = [u(y_u(uState,time) - m(time)*zbar(time) + zToAdd) + ...
            betas_u(uState,time)*gothicv_u(1,uState);v_u_endog(:,uState)];        
        % now v_u_endog_tointerp will preserve non-linearity in region 
        %  where agents are constrained
        
        this_c_u(:,uState) = interp1q(z_u_endog_tointerp,c_u_endog_tointerp,zThisIter);
        this_v_u(:,uState) = interp1q(z_u_endog_tointerp,v_u_endog_tointerp,zThisIter);
        c_u(zThisIterInd,uState,time) = this_c_u(zThisIterIndMapping,uState);
        v_u(zThisIterInd,uState,time) = this_v_u(zThisIterIndMapping,uState);          
    end
    
    % Finally, compute asset policy functions given inherited z, where I
    %  enforce >= zbar(time) to eliminate miniscule computational error when
    %  iterating distributions forward later
    znext_e(zThisIterInd,:,time) = max((1/m(time))*...
        (repmat(y_e(:,time)',sum(zThisIterInd),1) + ...
        repmat(zThisIter(zThisIterIndMapping),1,eStates) - c_e(zThisIterInd,:,time)),zbar(time));
    znext_u(zThisIterInd,:,time) = max((1/m(time))*...
        (repmat(y_u(:,time)',sum(zThisIterInd),1) + ...
        repmat(zThisIter(zThisIterIndMapping),1,uStates) - c_u(zThisIterInd,:,time)),zbar(time));
   
    % Solve for search
    ve_bigger = (this_v_e*Tue'-this_v_u) >= 0;
    sTemp = psiPrimeInv(repmat(pMat(:,time)',size(ve_bigger,1),1).*(this_v_e*Tue' - this_v_u)).*ve_bigger;
    this_s = min(repmat(1./(pMat(:,1)'),size(ve_bigger,1),1),sTemp);
    s(zThisIterInd,:,time) = this_s(zThisIterIndMapping,:);   
    
    % Iterate backward
    zLastIter = zThisIter(zbarPrevInd:end);
    v_eLastIter = this_v_e(zbarPrevInd:end,:);
    v_uLastIter = this_v_u(zbarPrevInd:end,:);
    c_eLastIter = this_c_e(zbarPrevInd:end,:);
    c_uLastIter = this_c_u(zbarPrevInd:end,:);
    time = time-1;   
end

% Check that initial z covers all agents at the initial stationary
% distribution
if zThisIter(length(zThisIter)) < z(zLastPosDens)
    display([tabs,'Policies solved backwards do not cover all agents in the initial period!']);
    zLastPol = find(find(zThisIter,1,'last')<=origGridInds,1)-1;
else
    zLastPol = find(zLastPosDens<=origGridInds,1);
end

% Select subset over original gridpoints
v_e = v_e(origGridInds,:,:);
v_u = v_u(origGridInds,:,:);
s = s(origGridInds,:,:);
c_e = c_e(origGridInds,:,:);
znext_e = znext_e(origGridInds,:,:);
c_u = c_u(origGridInds,:,:);
znext_u = znext_u(origGridInds,:,:);

timeTaken = toc;
% Display performance if toDisplay = 'on'
if strcmp(toDisplay,'on')
    disp([tabs,'Backward policy iteration complete: ',num2str(timeTaken),...
        's, ',num2str(Tbar-1), ' iterations, ',...
        num2str(timeTaken/(Tbar-1)), ' s/iteration.  d(v1,vTbar) = ',...
        num2str(d([v_e(:,:,1),v_u(:,:,1)],...
        [v_e(:,:,Tbar),v_u(:,:,Tbar)]))]);
end