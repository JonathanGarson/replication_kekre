function [a,aP,numA,Ta,numUI,Taui,Tua,pEpsilon] = genIncomeProcess(abar,...
    rhoP,sigmaP,numP,sigmaT,numT,zeta)
%% genIncomeProcess.m
%  Defines idiosyncratic productivity process

%% Start with productivity 

numA = numP*numT;

% Discretize permanent process (Rouwenhorst) -- following Kopecky-Suen (2010)

if numP >= 2
    bounds = sqrt(numP-1)*sqrt(sigmaP^2/(1-rhoP^2));
    nuP = linspace(-bounds,bounds,numP);
    val = (1+rhoP)/2;
    Tp = [val,1-val;1-val,val];
    for n=3:numP
        Tp = val*[Tp,zeros(n-1,1);zeros(1,n)] + ...
            (1-val)*[zeros(n-1,1),Tp;zeros(1,n)] + ...
            (1-val)*[zeros(1,n);Tp,zeros(n-1,1)] + ...
            val*[zeros(1,n);zeros(n-1,1),Tp];
        Tp(2:(n-1),:) = Tp(2:(n-1),:)/2;
    end
else
    nuP = 0;
    Tp = 1;
end
    
% Discretize transitory component (Gauss-Hermite) -- following Wouter Den Haan's notes

%  - first use code available at
%    https://www.mathworks.com/matlabcentral/fileexchange/26737-legendre-laguerre-and-hermite-gauss-quadrature/content/GaussHermite.m
%    to compute nodes and weights
i   = 1:numT-1;
a   = sqrt(i/2);
CM  = diag(a,1) + diag(a,-1);
[V L]   = eig(CM);
[x ind] = sort(diag(L));
V       = V(:,ind)';
weights = sqrt(pi) * V(:,1).^2;

% Now scale appropriately
epsilon = sqrt(2)*sigmaT*x';
pEpsilon = (1/sqrt(pi))*weights';

% Stack to form a and Tw
%   a ordered according to:
%    permanent shock
%    transitory shock

a = nan(1,numA);
aP = nan(1,numA);
Ta = nan(numA,numA);
for iP=1:numP
    for iT=1:numT
        i = (iP-1)*numT + iT;
        a(i) = exp(log(abar) + nuP(iP) + epsilon(iT)-(1/2)*sigmaT^2);
        aP(i) = exp(log(abar) + nuP(iP));
        Ta(i,:) = reshape(repmat(Tp(iP,:),numT,1),1,numP*numT).*...
            repmat(reshape(pEpsilon,1,numT),1,numP);
        Ta(i,:) = Ta(i,:)/sum(Ta(i,:)); % to eliminate miniscule computational error
    end
end

%% Now define relevant transitions

if zeta == 1
    numUI = numP;
else
    numUI = numP*2;
end

Tuia = zeros(numUI,numA);
Tua = zeros(numUI*numT,numA);
if zeta == 1
    Taui = reshape(repmat(reshape(eye(numP),1,numP^2),numT*1),numA,numP);
    Tua = eye(numA,numA);
else
    Taui = [(1-zeta)*reshape(repmat(reshape(eye(numP),1,numP^2),numT,1),numA,numP),...
        zeta*reshape(repmat(reshape(eye(numP),1,numP^2),numT,1),numA,numP)];
    Tua = [eye(numA,numA);eye(numA,numA)];
end