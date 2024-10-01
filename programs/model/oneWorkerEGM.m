function [snext,z_e_endog,c_e_endog,v_e_endog,gothicv_e,...
    z_u_endog,c_u_endog,v_u_endog,gothicv_u] = ...
    oneWorkerEGM(znext,v_enext,v_unext,c_enext,c_unext,...
    m,pVecnext,deltaVec,Tee,Teu,Tuenext,Tuu,y_e,y_u,betas_e,betas_u)
%% oneWorkerEGM.m
% Performs one iteration of endogenous-grid-point-augmented 
%  value function iteration for workers
% Note that if c_enext, c_unext are set to -99, this means that
%  numerical derivatives of gothic v's are computed, rather than Envelope

global xi sigma;

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
uPrimeInv = @(val) val.^(-1/sigma);
uPrime = @(c) c.^(-sigma); % used for alt maximization step (using Envelope)
psi = @(eval) (eval).^(xi+1);
psiPrimeInv = @(val) (val./(xi+1)).^(1/xi);

% Define dimensions of interest
gridPoints = length(znext);
eStates = size(v_enext,2);
uStates = size(v_unext,2);

% Next period's labor supply maximization
ve_bigger = (v_enext*Tuenext' - v_unext) >= 0;
snextTemp = psiPrimeInv(repmat(pVecnext,gridPoints,1).*(v_enext*Tuenext' - v_unext)).*ve_bigger;
snext = min(repmat(1./pVecnext,gridPoints,1),snextTemp);   

% Implied job-finding probabilities for the initially unemployed next
% period (as a function of next period's state variables)
Tempnext = repmat(pVecnext,gridPoints,1).*snext;

% Implied separation rates at the end of this period
Tsep = repmat(deltaVec,gridPoints,1);

% Gothic v computation (a function of next period's assets but this 
%  period's other state variables)
gothicv_u = ...
    ((Tempnext.*(v_enext*Tuenext')) + ...
    (1-Tempnext).*v_unext - ...
    psi(snext))*Tuu';
gothicv_e = ...
    (1-Tsep).*(v_enext*Tee') + ...
    Tsep.*((Tempnext.*(v_enext*Tuenext') + ...
           (1-Tempnext).*v_unext - ...
           psi(snext))*Teu');

% Gothic v slope computation 
if c_enext == -99
    gothicv_u_slopes = getSlopes(gothicv_u,znext);    
    gothicv_e_slopes = getSlopes(gothicv_e,znext);
else
    gothicv_u_slopes = ...
        (Tempnext.*(uPrime(c_enext)*Tuenext') + ...
        (1-Tempnext).*uPrime(c_unext))*Tuu';
    gothicv_e_slopes = ...
        (1-Tsep).*(uPrime(c_enext)*Tee') + ...
        Tsep.*((Tempnext.*(uPrime(c_enext)*Tuenext') + ...
               (1-Tempnext).*uPrime(c_unext))*Teu');
end

% Invert optimal consumption using Euler
c_e_endog = uPrimeInv(repmat(betas_e,gridPoints,1).*(1/m).*gothicv_e_slopes);
c_u_endog = uPrimeInv(repmat(betas_u,gridPoints,1).*(1/m).*gothicv_u_slopes);

% Compute inherited debt (endogenous gridpoints)
z_e_endog = c_e_endog + repmat(m*znext,1,eStates) - repmat(y_e,gridPoints,1);
z_u_endog = c_u_endog + repmat(m*znext,1,uStates) - repmat(y_u,gridPoints,1);

% Define new value functions at endogenous gridpoints
v_e_endog = u(c_e_endog) + repmat(betas_e,gridPoints,1).*gothicv_e;
v_u_endog = u(c_u_endog) + repmat(betas_u,gridPoints,1).*gothicv_u; 