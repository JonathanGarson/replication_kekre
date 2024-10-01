function [b,bUI] = genTransfers(numBetas,numP,numT,dbarbar,dbar,rr,...
    w,wP,omegaVar1,omegaVar2,zeta,uimax)
%% genTransfers.m
% Defines transfers to unemployed

%% Define UI

% UI state ordered according to:
%  eligibility+take-up (if applicable)
%  permanent level of income

wP_for_b = reshape(repmat(reshape(reshape(wP,numP,numBetas),...
    1,numP*numBetas),numT,1),1,numT*numP*numBetas);
    % is a numT*numP*numBetas vector of wages at the average value
    %  of transitory productivity

if zeta == 1
    bUI = min(reshape([rr*repmat(wP_for_b,dbar,1);...
        sa*ones(dbarbar-dbar,numT*numP*numBetas)],1,(dbarbar+1)*numT*numP*numBetas),uimax);
else
    bUI = min(reshape([zeros((dbarbar+1)*numT*numP,numBetas);...
        reshape([rr*repmat(wP_for_b,dbar,1);...
        zeros(dbarbar+1-dbar,numT*numP*numBetas)],(dbarbar+1)*numT*numP,numBetas)],...
        1,(dbarbar+1)*numT*numP*2*numBetas),uimax);
end
    
%% Define endowment when unemployed

endow1 = omegaVar1*w;
endow2 = omegaVar2*w;
if zeta == 1
    endow = reshape([repmat(endow1,dbar,1);repmat(endow2,dbarbar+1-dbar,1)],1,(dbarbar+1)*numT*numP*numBetas);
else
    endow1 = reshape(endow1,numT*numP,numBetas);
    endow2 = reshape(endow2,numT*numP,numBetas);
    endow = reshape([repmat(reshape([endow2;endow1],1,numT*numP*2*numBetas),dbar,1);...
        repmat(reshape([endow2;endow2],1,numT*numP*2*numBetas),dbarbar+1-dbar,1)],1,(dbarbar+1)*numT*numP*2*numBetas);
end

%% Define total transfers to unemployed

b = bUI + endow;