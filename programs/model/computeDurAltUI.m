function [durAltUI,sStatAltUI,c_eStatAltUI,znext_eStatAltUI,...
    c_uStatAltUI,znext_uStatAltUI] = ...
    computeDurAltUI(zStat,phi_u_tilde0,states0,m,pVec,deltaVec,Tee,Teu,...
    Tue,Tuu,y_e,y_uAltUI,betas_e,betas_u,eps_v,toDisplay,tol)
%% computeDurAltUI.m
% Computes avg duration given alternative UI path, and thus alternative
%  y_uAltUI.

global zbar numBetas numA

zAltUI = zStat;
gridPoints = length(zAltUI);

% Characterize optimal policies and value functions until
% znext_e doesn't hit the artifical upper bound
upperBoundHit = 1;
while upperBoundHit == 1
    [~,~,~,sAltUI,c_eAltUI,znext_eAltUI,c_uAltUI,znext_uAltUI] = ...
        statPoliciesEGM(m,pVec,deltaVec,Tee,Teu,Tue,Tuu,...
        y_e,y_uAltUI,betas_e,betas_u,eps_v,toDisplay,zAltUI);

    eSaver = numBetas*numA;
    % theoretically, this holds the index of the employed
    %  agent who saves the most at least at high levels of wealth (highest
    %  beta and highest w)
    if zAltUI(gridPoints) - znext_eAltUI(gridPoints,eSaver) > tol
        upperBoundHit = 0;
    else
        slope = (znext_eAltUI(gridPoints,eSaver) - znext_eAltUI(gridPoints-1,eSaver))/...
            (zAltUI(gridPoints)-zAltUI(gridPoints-1));
        interpolatedUpper = (znext_eAltUI(gridPoints,eSaver) - slope*zAltUI(gridPoints))/...
            (1-slope); % so that by linear interpolation, 
        % znext_e(interpolatedUpper,eSaver) should = interpolatedUpper

        % Then add 10% more gridPoints to z up to 1.1*interpolatedUpper
        newZAltUI = linspace(zAltUI(gridPoints),1.1*interpolatedUpper,round(0.1*gridPoints)+1);
        zAltUI = [zAltUI;newZAltUI(2:length(newZAltUI))'];
        gridPoints = length(zAltUI);

        disp(['Iterative algorithm upper bound hit: expanding it to ',...
            num2str(1.1*interpolatedUpper)]);            
    end
end

% Alert if p(theta)s = 1 for some unemployed agents
if sum(sum(repmat(pVec,gridPoints,1).*sAltUI >= 1)) > 0
    disp('Some agents at max effort!');
    pause;
end

% Interpolate key policies over zStat
sStatAltUI = interp1(zAltUI,sAltUI,zStat);
c_eStatAltUI = interp1(zAltUI,c_eAltUI,zStat);
znext_eStatAltUI = interp1(zAltUI,znext_eAltUI,zStat);
c_uStatAltUI = interp1(zAltUI,c_uAltUI,zStat);
znext_uStatAltUI = interp1(zAltUI,znext_uAltUI,zStat);

% Again, to avoid error
znext_eStatAltUI(znext_eStatAltUI <= zbar) = zbar;
znext_uStatAltUI(znext_uStatAltUI <= zbar) = zbar;

durAltUI = simDurationAlt(phi_u_tilde0,states0,zStat,znext_uStatAltUI,...
    sStatAltUI,pVec,Tuu,toDisplay);