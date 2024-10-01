function dur = simDurationAlt(phi_u0,states0,z,znext_u,s,pVec,Tuu,toDisplay)
%% simDurationAlt.m
% Simulates the average duration of unemployment in an economy with 
%  steady-state distribution of unemployed agents at *beginning* of
%  unemployment spell phi_u0, states0 denoting the indices among 1:uStates
%  at which those unemployed agents live, job-finding rate pVec, 
%  {u}->{u} transition matrix Tuu, and policies znext_u and s defined over z. 
% Lower bound on duration is 1 (regain employment in the very next period) 
%  and upper bound is infinity. 

timeTaken = tic;

% Define dimensions of interest
gridPoints = length(z);
uStates = size(znext_u,2);

% Job-finding probabilities for the initially unemployed 
Temp = repmat(pVec,gridPoints,1).*s;

% Transition matrices for assets
Tznext_u = cell(1,uStates);

for uState=1:uStates
    Tznext_u{1,uState} = genSparseTransition(znext_u(:,uState),z);
end

% Define simulation parameters
T = 100; % a number > dbarbar
tol = 1e-10; % measure of remaining agents below which algorithm stops 

% Initialize objects simulated forward
fracUnemp = nan(1,T);
fracUnemp(1) = 1;

initUnempDist = zeros(gridPoints,uStates);
initUnempDist(:,states0) = phi_u0;
initUnempDist_znext = nan(gridPoints,uStates);
initUnempDist_next = nan(gridPoints,uStates);

% Simulate!
t = 1;
while fracUnemp(t) > tol && t < T
    for uState=1:uStates
        initUnempDist_znext(:,uState) = (initUnempDist(:,uState)'*Tznext_u{1,uState})';
    end
    initUnempDist_next = (1-Temp).*(initUnempDist_znext*Tuu)*fracUnemp(t);
    fracUnemp(t+1) = sum(sum(initUnempDist_next));
    initUnempDist = (1/fracUnemp(t+1))*initUnempDist_next;
    t = t+1;
end

% Compute fraction by duration
fracByDur = fracUnemp(1:T-1)-fracUnemp(2:T);

% Plot and return
dur = nansum(fracByDur.*(1:(T-1)));
if toDisplay == 1
    figure;
    bar(fracByDur(1:12));
    title(['Fraction of agents by duration, avg = ',num2str(dur)]);
    drawnow;
end