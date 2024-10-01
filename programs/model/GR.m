function GR(Tbar,T,saveAll,initShock,shockPersist,...
    deltaBarShocks,deltaBarPersist,dbarPathMonths,dbarShockMonths,Hss,dH,urdata,numTabs,...
    transParams,...
    RCEName,toName,invariantParams,varyingParams,equityShares_e,equityShares_u,...
    eps,epsMacro,znextDiff,goodsDiff,tDiff,...
    mDiff,equityValDiff,dvacpostDiff,dwageDiff,cTransSS,...
    startT,startName,predetShocks)
%% GR.m
% Simulates Great Recession using discount factor shocks
%  Tbar = number of total periods until economy assumed back in SS 
%  T = number of periods in Great Recession simulation
%  saveAll = 1 if additional output should be saved
%  initShock = initial value of beta shock
%  shockPersist = persistence of beta process
%  deltaBarShocks = aggregate separation rate shocks
%  deltaBarPersist = persistence of agg separation rate process
%  dbarPathMonths = matrix of path of UI duration along simulation
%  dbarShockMonths = vector of months featuring UI shocks in simulation
%  Hss = steady-state Jacobian matrix
%  dH = change in Jacobian from steady-state to time 1
%  urdata = observed path of unemployment rate vs. April 2008 used to
%   calibrate simulation
%  numTabs = number of tabs in output
%  transParams = parameters specific to transitional dynamics
%  RCEName = filename from which stationary RCE should be called
%  toName = filename to which simulation output should be saved
%  invariantParams = parameters assumed to be invariant over transitions
%  varyingParams = parameters which may vary over transitions
%  equityShares_{i}SS = assumed fraction of portfolio held in firm equity 
%   in the stationary RCE for i \in {e,u} agents.  A gridpointsStat x iStates matrix.
%  eps = convergence criteria to solve for equilibrium in response to each shock
%  epsMacro = convergence criteria to calibrate shocks to match urdata
%  znextDiff = znextTrans - zsupplyTrans used as benchmark for converg
%  goodsDiff = cTrans - prodTrans used as benchmark for converg
%  tDiff = tTrans - t used as benchmark for converg
%  mDiff = mTrans - m used as benchmark for converg
%  equityValDiff = equityValTrans - equityVal used as benchmark for converg
%  dvacpostDiffSS = vacpostTrans - sfNewHiresTrans used as a benchmark for converg
%  dwageDiffSS = wTrans - w used as a benchmark for converg
%  cTransSS = cTrans used as a benchmark for cTrans in Taylor rule 
%   (only relevant if phiC > 0)
%  startT = starting period of simulation (between 0 and T-1)
%  startName = (optional) starting dynamics of simulation if startT > 0
%  predetShocks = (optional) pre-determined beta shocks to use, rather than calibrate

tabs = '.';
for i=1:numTabs
    tabs = sprintf('   %s',tabs);
end

load(RCEName,'m','p_e','eStates','uStates',...
    'numBetas','xi','sigma','numA','numP','rhoP','sigmaP','numT','sigmaT',...
    'numUI','Tua','dbarbar',...
    'k','eta','pEpsilon','phi',...
    'betas','Ta','Taui','zeta','dbar','rr','deltaBar',...
    'eps_delta_a','eps_delta_beta','mbar','lambda',...
    'abar','a','aP','zbar','znext_g','omega0','omega1','omega2',...
    'qStat','z','p_e_tilde','phi_e_tilde','phi_u_tilde','gridpointsStat',...
    'dbarbar','theta','t','w','mu','RCEName');

shocks = zeros(1,T);
MPath = zeros(1,Tbar);
mPath = zeros(1,Tbar);
p_ePath = zeros(1,Tbar);
p_uPath = zeros(dbarbar+1,Tbar);
thetaPath = zeros(1,Tbar);
p_e_tildePath = zeros(1,Tbar);
sbarPath = zeros(1,Tbar);
avgWPath = zeros(1,Tbar);
wPath = zeros(eStates,Tbar);
cPath = zeros(1,Tbar);
PiPath = zeros(1,Tbar);
MForward = zeros(Tbar,T);
PiForward = zeros(Tbar,T);
znextPath = zeros(1,Tbar);

if startT == 0
    if isempty(predetShocks)
        shocks(1) = initShock;
    else
        shocks(1) = predetShocks(1);
    end
    lastBetasTrans = repmat(betas',1,Tbar);
    lastDeltaBarTrans = repmat(deltaBar',1,Tbar);
    starting = [];
    starting0 = struct('m',m*ones(1,Tbar),'theta',theta*ones(1,Tbar),...
        't',t*ones(1,Tbar),'equityVal',(1/m)*qStat,...
        'w',repmat(w',1,Tbar),'mu',mu*ones(1,Tbar));
    Hsticky = Hss;
    dHsticky = dH;
    Hzlb = Hsticky;
    lastAtZLB = 0;
else
    load([startName,'.mat']);
    time = startT;
    lastAtZLB = MTrans(1) == 1;
    if isempty(predetShocks)
        if time <= 8 && ~lastAtZLB
            shocks(time+1) = 0.25*(urdata(time+1) - (p_e-p_eTrans(2)))/dunemp;
        else
            shocks(time+1) = (urdata(time+1) - (p_e-p_eTrans(2)))/dunemp;
        end
    else
        shocks(time+1) = predetShocks(time+1);
    end
    lastBetasTrans = betasTrans(:,2:end);        
    lastDeltaBarTrans = deltaBarTrans(2:end);
    starting0 = struct('m',this_m(2:end),'theta',this_theta(2:end),...
        't',this_t(2:end),'equityVal',dist2.initEquityVal,...
        'w',this_w(:,2:end),'mu',this_mu(2:end));
    starting = struct('m',starting0.m + shocks(time+1)*dm(1:end-1),...
        'theta',starting0.theta + shocks(time+1)*dtheta(1:end-1),...
        't',starting0.t + shocks(time+1)*dt(1:end-1),...
        'equityVal',starting0.equityVal + shocks(time+1)*dequityVal,...
        'w',starting0.w + shocks(time+1)*dw(:,1:end-1),...
        'mu',starting0.mu + shocks(time+1)*dmu(1:end-1));                
    Hsticky = Hss;
    dHsticky = dH;            
    Hzlb = Hsticky + dHsticky*(p_e-p_ePath(time));         
end

warning('off','all');
isError = 0;
for time=startT+1:T
    disp([tabs,'*** time = ',num2str(time)]);
    converged = 0;
    if time == 1
        dist1 = [];
        p_e1 = p_e;      
        p_e_tildePath(1) = p_e_tilde;
    else
        dist1 = dist2;
        p_e1 = p_eTrans(2);   
        % Set znextDiff, goodsDiff, tDiff, mDiff, equityValDiff, dvacpostDiff, 
        %  and dwageDiff from the prior iteration in absence of any more shocks
        znextDiff = znextTrans(2:end) - zsupplyTrans(2:end);
        goodsDiff = cTrans(2:end) - prodTrans(2:end);
        tDiff = tTrans(2:end) - this_t(2:end);
        mDiff = mTrans(2:end) - this_m(2:end);
        equityValDiff = 0;
        dvacpostDiff = vacpostTrans(2:end) - sfNewHiresTrans(2:end);  
        dwageDiff = wTrans(:,2:end)-this_w(:,2:end);
        p_e_tildePath(time) = dist2.p_e_tilde1;

        % Update equity holding rule
        [equityShares] = ...
            initializeEquity(0.35,...
            dist1.initEquityVal,dist1.z,[dist1.p_e_tilde1*dist1.phi_e_tilde1,(1-dist1.p_e_tilde1)*dist1.phi_u_tilde1]);
        equityShares_e = repmat(equityShares,1,eStates);
        equityShares_u = repmat(equityShares,1,uStates);
    end
    % Simulate UI shock and separation rate shock (if there are any), and 
    %  we will then be solving for a fixed point in a beta shock).  
    deltaBarTrans = lastDeltaBarTrans;
    if (dbarShockMonths(time) && isempty(predetShocks)) || deltaBarShocks(time) ~= 0
        betasTrans = lastBetasTrans; 
        deltaBarAR1 = zeros(1,Tbar-time+1);
        deltaBarAR1(1) = deltaBarShocks(time);
        for time2=2:180
            deltaBarAR1(time2) = deltaBarPersist*deltaBarAR1(time2-1);
        end
        deltaBarTrans = lastDeltaBarTrans + deltaBarAR1;
        varyingParams = setfield(varyingParams,'betas',betasTrans);
        varyingParams = setfield(varyingParams,'deltaBar',deltaBarTrans);
        varyingParams = setfield(varyingParams,'Tbar',Tbar-time+1);
        varyingParams = setfield(varyingParams,'dbar',dbarPathMonths(1:Tbar-time+1,time)');
        try
            [this_m,this_t,this_theta,this_equityVal,this_w,this_mu,~,...
                ~,~,~,~,~,~,~,~,~,...
                ~,~,~,~,~,~,p_eTrans,~,~,~,~] = ...
                Transition(transParams,...
                RCEName,[],invariantParams,varyingParams,...
                equityShares_e,equityShares_u,eps,numTabs,znextDiff,goodsDiff,...
                tDiff,mDiff,equityValDiff,dvacpostDiff,dwageDiff,cTransSS(time:end),starting0,Hzlb,dist1,0,0);        
            starting0 = struct('m',this_m,'theta',this_theta,...
                't',this_t,'equityVal',this_equityVal,...
                'w',this_w,'mu',this_mu);
            p_e1 = p_eTrans(1);
            if time > 1
                if isempty(predetShocks)
                    if time <= 8 && ~lastAtZLB
                        shocks(time) = 0.25*(urdata(time) - (p_e-p_eTrans(1)))/dunemp;
                    else
                        shocks(time) = (urdata(time) - (p_e-p_eTrans(1)))/dunemp;
                    end 
                end
                starting = struct('m',starting0.m + shocks(time)*dm(1:end-1),...
                    'theta',starting0.theta + shocks(time)*dtheta(1:end-1),...
                    't',starting0.t + shocks(time)*dt(1:end-1),...
                    'equityVal',starting0.equityVal + shocks(time)*dequityVal,...
                    'w',starting0.w + shocks(time)*dw(:,1:end-1),...
                    'mu',starting0.mu + shocks(time)*dmu(1:end-1));            
            end
        catch err
            isError = 1;
            break;
        end
    end    
    if isError
        break;
    end
    while ~converged                    
        % Simulate beta shock
        betaAR1 = zeros(1,Tbar-time+1);
        betaAR1(1) = shocks(time);
        for time2=2:180
            betaAR1(time2) = shockPersist*betaAR1(time2-1);
        end
        betasTrans = lastBetasTrans + repmat(betaAR1,numBetas,1);
        varyingParams = setfield(varyingParams,'betas',betasTrans);
        varyingParams = setfield(varyingParams,'deltaBar',deltaBarTrans);
        varyingParams = setfield(varyingParams,'Tbar',Tbar-time+1);
        varyingParams = setfield(varyingParams,'dbar',dbarPathMonths(1:Tbar-time+1,time)');

        try
            if time <= 8 && ~lastAtZLB
                [this_m,this_t,this_theta,this_equityVal,this_w,this_mu,~,...
                    znextTrans,zsupplyTrans,cTrans,prodTrans,tTrans,mTrans,~,vacpostTrans,sfNewHiresTrans,...
                    wTrans,~,~,~,PiTrans,MTrans,p_eTrans,p_uTrans,sbarTrans,avgW,dist2] = ...
                    Transition(transParams,...
                    RCEName,[],invariantParams,varyingParams,...
                    equityShares_e,equityShares_u,eps,numTabs,znextDiff,goodsDiff,...
                    tDiff,mDiff,equityValDiff,dvacpostDiff,dwageDiff,cTransSS(time:end),starting,Hzlb,dist1,0.5,0,0);
            else
                [this_m,this_t,this_theta,this_equityVal,this_w,this_mu,~,...
                    znextTrans,zsupplyTrans,cTrans,prodTrans,tTrans,mTrans,~,vacpostTrans,sfNewHiresTrans,...
                    wTrans,~,~,~,PiTrans,MTrans,p_eTrans,p_uTrans,sbarTrans,avgW,dist2] = ...
                    Transition(transParams,...
                    RCEName,[],invariantParams,varyingParams,...
                    equityShares_e,equityShares_u,eps,numTabs,znextDiff,goodsDiff,...
                    tDiff,mDiff,equityValDiff,dvacpostDiff,dwageDiff,cTransSS(time:end),starting,Hzlb,dist1,0,0,0);                
            end
        catch err
            isError = 1;
            break;
        end
        if isError
            break;
        end            
        dunemp = (p_e1 - p_eTrans(1))/shocks(time);
        dm = (this_m - starting0.m)/shocks(time);
        dt = (this_t - starting0.t)/shocks(time);
        dtheta = (this_theta - starting0.theta)/shocks(time);
        dequityVal = (this_equityVal - starting0.equityVal)/shocks(time);
        dw = (this_w - starting0.w)/shocks(time);
        dmu = (this_mu - starting0.mu)/shocks(time);

        % Assess consistency with observed unemployment rate data and update
        if isempty(predetShocks) && abs(p_e - p_eTrans(1) - urdata(time)) >= epsMacro
            if time <= 8 && ~lastAtZLB
                shocks(time) = 0.5*shocks(time) + 0.5*(urdata(time) - (p_e-p_e1))/dunemp;
            else
                shocks(time) = (urdata(time) - (p_e-p_e1))/dunemp;
            end
            starting = struct('m',starting0.m + shocks(time)*dm,...
                'theta',starting0.theta + shocks(time)*dtheta,...
                't',starting0.t + shocks(time)*dt,...
                'equityVal',starting0.equityVal + shocks(time)*dequityVal,...
                'w',starting0.w + shocks(time)*dw,...
                'mu',starting0.mu + shocks(time)*dmu);                
        else
            converged = 1;
            MPath(time) = MTrans(1);
            mPath(time) = this_m(1);
            p_ePath(time) = p_eTrans(1);
            p_uPath(:,time) = p_uTrans(:,1);
            thetaPath(time) = this_theta(1);
            sbarPath(time) = sbarTrans(1);
            avgWPath(time) = avgW(1);
            wPath(:,time) = this_w(:,1);
            cPath(time) = cTrans(1);
            PiPath(time) = PiTrans(1);
            MForward(1:(Tbar-time+1),time) = MTrans';
			PiForward(1:(Tbar-time+1),time) = PiTrans';
			znextPath(time) = znextTrans(1);

            if time < T
                save([toName,'.mat'],'shocks',...
                    'MPath','mPath','p_ePath','p_uPath','thetaPath',...
                    'p_e_tildePath','sbarPath','avgWPath','wPath','cPath','PiPath','MForward',...
					'PiForward','znextPath',...
                    'p_eTrans','p_uTrans','sbarTrans','betasTrans','deltaBarTrans','MTrans',...
                    'dunemp','dm','dtheta','dt','dequityVal','dw','dmu',...
                    'this_m','this_theta','this_t','dist2','this_w','this_mu',...
                    'starting0','Hsticky','dHsticky','dist2',...
                    'znextTrans','zsupplyTrans','cTrans','prodTrans','tTrans','mTrans',...
                    'vacpostTrans','sfNewHiresTrans','wTrans','lastBetasTrans','lastDeltaBarTrans');
                lastAtZLB = MTrans(1) == 1;
                if isempty(predetShocks)
                    if time <= 8 && ~lastAtZLB
                        shocks(time+1) = 0.25*(urdata(time+1) - (p_e-p_eTrans(2)))/dunemp;
                    else
                        shocks(time+1) = (urdata(time+1) - (p_e-p_eTrans(2)))/dunemp;
                    end
                else
                    shocks(time+1) = predetShocks(time+1);
                end
                lastBetasTrans = betasTrans(:,2:end);   
                lastDeltaBarTrans = deltaBarTrans(2:end);
                starting0 = struct('m',this_m(2:end),'theta',this_theta(2:end),...
                    't',this_t(2:end),'equityVal',dist2.initEquityVal,...
                    'w',this_w(:,2:end),'mu',this_mu(2:end));
                starting = struct('m',starting0.m + shocks(time+1)*dm(1:end-1),...
                    'theta',starting0.theta + shocks(time+1)*dtheta(1:end-1),...
                    't',starting0.t + shocks(time+1)*dt(1:end-1),...
                    'equityVal',starting0.equityVal + shocks(time+1)*dequityVal,...
                    'w',starting0.w + shocks(time+1)*dw(:,1:end-1),...
                    'mu',starting0.mu + shocks(time+1)*dmu(1:end-1));
                % update Jacobian for next iteration
                selectCols = ones(1,(4+eStates)*(Tbar-time)+1);
                selectRows = ones(1,(7+eStates)*(Tbar-time)+1);
                for i=1:(4+eStates)
                    selectCols(1+i*(Tbar-time)) = 0;
                end
                for i=1:(7+eStates)
                    selectRows(1+i*(Tbar-time)) = 0;
                end 
                Hsticky = Hsticky(:,selectCols==1);
                Hsticky = Hsticky(selectRows==1,:);
                if time > 1
                    dHsticky = dHsticky(:,selectCols==1);
                    dHsticky = dHsticky(selectRows==1,:);
                end            
                Hzlb = Hsticky + dHsticky*(p_e-p_ePath(time));            
            end
        end            
    end
    if isError
        break;
    end    
end

if ~isError
    MPath(T+1:Tbar) = MTrans(2:end);
    mPath(T+1:Tbar) = this_m(2:end);
    p_ePath(T+1:Tbar) = p_eTrans(2:end);
    p_uPath(:,T+1:Tbar) = p_uTrans(:,2:end);
    thetaPath(T+1:Tbar) = this_theta(2:end);
    p_e_tildePath(T+1:Tbar) = dist2.p_e_tilde(2:end);
    sbarPath(T+1:Tbar) = sbarTrans(2:end);
    avgWPath(T+1:Tbar) = avgW(2:end);
    wPath(:,T+1:Tbar) = this_w(:,2:end);
    cPath(T+1:Tbar) = cTrans(2:end);
    PiPath(T+1:Tbar) = PiTrans(2:end);
	znextPath(T+1:Tbar) = znextTrans(2:end);	
	
    
    if ~saveAll
        save([toName,'.mat'],'shocks',...
            'MPath','mPath','p_ePath','p_uPath','thetaPath',...
            'p_e_tildePath','sbarPath','avgWPath','wPath','cPath','PiPath','MForward',...
			'PiForward','znextPath');
    else
        save([toName,'.mat'],'shocks',...
            'MPath','mPath','p_ePath','p_uPath','thetaPath',...
            'p_e_tildePath','sbarPath','avgWPath','wPath','cPath','PiPath','MForward',...
			'PiForward','znextPath',...
            'p_eTrans','p_uTrans','sbarTrans','betasTrans','deltaBarTrans','MTrans',...
            'dunemp','dm','dtheta','dt','dequityVal','dw','dmu',...
            'this_m','this_theta','this_t','dist2','this_w','this_mu',...
            'starting0','Hsticky','dHsticky','dist2',...
            'znextTrans','zsupplyTrans','cTrans','prodTrans','tTrans','mTrans',...
            'vacpostTrans','sfNewHiresTrans','wTrans','lastBetasTrans','lastDeltaBarTrans');
    end
else
    disp([tabs,'Error: ',err.identifier]);
end