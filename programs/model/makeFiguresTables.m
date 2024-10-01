%% makeFiguresTables.m
% Produces all figures and tables involving model simulations

tic;

%% Define names of stationary RCEs and transitional dynamics used to produce figures/tables

% Directory containing all stationary RCEs and transitional dynamics
outputName = 'output';
if ~exist(outputName,'dir')
    mkdir(outputName);
end

% Directory in which all figures and tables are produced
if ~exist([outputName,'/figures_and_tables'],'dir')
    mkdir([outputName,'/figures_and_tables']);
end

% Directory containing selected data series
dataName = 'output';
macroGR = 'macrodata.csv';
fvUI = 'ui_weeks_fv.csv';
sepratedata = 'seprate.csv';
seprateshocks = 'seprate05081214.csv';

% Names of stationary RCEs and targeted moments
calNames = {'RCE_baseline.mat','RCE_noeps.mat','RCE_durElast.mat','RCE_liquid.mat'};
calTables = {'Table1.txt','TableA5.txt','TableA6.txt','TableA7.txt'};
r_target = [0.02,0.02,0.02,0.02];
wealth_to_moincome_target = [66.0,66.0,66.0,3.7];
dwealth_ue_to_moincome_target = [-47.5,nan,-47.5,-2.7];
share_wealth_lt0_target = [0.08,0.08,0.08,0.26];
mpc_lt75k_target = [0.21,0.21,0.21,0.21];
frac_unemp_UI_target = [0.39,0.39,0.39,0.39];
max_rr_target = [0.6,0.6,0.6,0.6];
hhinc_ue_UI_target = [0.76,0.76,0.76,0.76];
hhinc_ue_NoUI_target = [0.55,0.55,0.55,0.55];
ur_target = [0.05,0.05,0.05,0.05];
share_ur_lt_target = [0.17,0.17,0.17,0.17];
weuCoeff_target = [-0.012,nan,-0.012,-0.012];
elast_dur_target = [0.1,0.1,0.4,0.1];
conv_tightness_target = [0.634,0.634,0.634,0.634];
mowage_per_vac_target = [0.108,0.108,0.108,0.108];
RCE_baselineName = 'RCE_baseline.mat';
RCE_baseline_dbarbar22Name = 'RCE_baseline_dbarbar22.mat';

% Identification of selected parameters in baseline calibration
identifyFileNames = {'ident_Delta',...
    'ident_eps_delta_beta','ident_eps_delta_a',...
    'ident_lambda'};
identifyParamNames = {'Delta','eps_delta_beta','eps_delta_a','lambda'};
identifyParams = {[0.0033 0.0037 0.0041 0.0045 0.0049],...
    [-3.95 -4.15 -4.35 -4.55 -4.75],...
    [-0.008 -0.009 -0.01 -0.011 -0.012],...
    [-0.08 -0.10 -0.12 -0.14 -0.16]};
identifyParamsChosen = {0.0045,-4.55,-0.011,-0.14};
identifyParamLabels = {'$\Delta^{\beta}$','$\epsilon^{\delta}_{\beta}$',...
    '$\epsilon^{\delta}_a$','$\lambda$'};
identifyTargetNames = {'mpc_lt75k','dwealth_ue_to_moincome',...
    'weuCoeff','share_ur_lt'};
identifyTargets = {0.21,-47.5,-0.012,0.17};
identifyTargetLabels = {'\textbf{Mean quarterly MPC to \$500 rebate$^\ast$}',...
    '\textbf{Mean (U-E) wealth / monthly HH income}',...
    '\textbf{EU probability on log wage}',...
    '\textbf{Fraction unemployed w/ duration $>$ 6 mos}'};

% Baseline effects of UI starting from stationary RCE
initName = 'trans_baseline_noshocks_sticky';
peName = 'trans_baseline_pe_sticky';
finName_flex = 'trans_baseline_flex';
finName_sticky = 'trans_baseline_sticky';
finName_m = 'trans_baseline_fixedr_sticky';
finName_zlb = 'trans_baseline_fixedi_sticky';

% Sensitivity analysis starting from stationary RCE
RCE_Name_noEps = 'RCE_noEps';
initName_noEps = 'trans_noEps_noshocks_sticky';
finName_zlb_noEps = 'trans_noEps_fixedi_sticky';

RCE_Name_liquid = 'RCE_liquid';
initName_liquid = 'trans_liquid_noshocks_sticky';
finName_zlb_liquid = 'trans_liquid_fixedi_sticky';

finName_flex_iota0p9 = 'trans_baseline_iota0p9_flex';
finName_sticky_iota0p9 = 'trans_baseline_iota0p9_sticky';
finName_zlb_iota0p9 = 'trans_baseline_iota0p9_fixedi_sticky';

RCE_Name_durElast = 'RCE_durElast';
initName_durElast = 'trans_durElast_noshocks_sticky';
finName_flex_durElast = 'trans_durElast_flex';
finName_sticky_durElast = 'trans_durElast_sticky';
finName_zlb_durElast = 'trans_durElast_fixedi_sticky';

durLabels = [3,6,9,12,15,18,21,24,27,30,33];
durExitZLB = 18;
durFinNames = {'trans_baseline_dur_3_fixedi_sticky',...
    'trans_baseline_dur_6_fixedi_sticky',...
    'trans_baseline_dur_9_fixedi_sticky',...
    'trans_baseline_dur_12_fixedi_sticky',...
    'trans_baseline_dur_15_fixedi_sticky',...
    'trans_baseline_dur_18_fixedi_sticky',...
    'trans_baseline_dur_21_fixedi_sticky',...
    'trans_baseline_dur_24_fixedi_sticky',...
    'trans_baseline_dur_27_fixedi_sticky',...
    'trans_baseline_dur_30_fixedi_sticky',...
    'trans_baseline_dur_33_fixedi_sticky'};

dbarLabels = [3,6,9,12];
dbarInitName = 'trans_baseline_dbarbar22_noshocks_sticky';
dbarFinNames = {'trans_baseline_dbar_3_fixedi_sticky',...
    'trans_baseline_dbar_6_fixedi_sticky',...
    'trans_baseline_dbar_9_fixedi_sticky',...
    'trans_baseline_dbar_12_fixedi_sticky'};
dbarFinNames1 = {'trans_baseline_dbar_3_1_fixedi_sticky',...
    'trans_baseline_dbar_6_1_fixedi_sticky',...
    'trans_baseline_dbar_9_1_fixedi_sticky',...
    'trans_baseline_dbar_12_1_fixedi_sticky'};
dbarFinNames12 = {'trans_baseline_dbar_3_12_fixedi_sticky',...
    'trans_baseline_dbar_6_12_fixedi_sticky',...
    'trans_baseline_dbar_9_12_fixedi_sticky',...
    'trans_baseline_dbar_12_12_fixedi_sticky'};

initName_zeta = 'trans_baseline_zeta_noui_fixedi_sticky';
finName_zlb_zeta = 'trans_baseline_zeta_fixedi_sticky';
finName_zlb_debt = 'trans_baseline_debt_fixedi_sticky';
finName_zlb_rr = 'trans_baseline_rr_fixedi_sticky';

finName_sticky_fiscal = 'trans_baseline_fiscal_sticky';
finName_zlb_fiscal = 'trans_baseline_fiscal_fixedi_sticky';
finName_sticky_fiscal_debt = 'trans_baseline_fiscal_debt_sticky';
finName_zlb_fiscal_debt = 'trans_baseline_fiscal_debt_fixedi_sticky';

% Macro shocks starting from stationary RCE
macroShockFileNames = {'trans_baseline_iota0p94_beta_sticky',...
    'trans_baseline_iota0p94_zbar_sticky',...
    'trans_baseline_iota0p94_abar_sticky',...
    'trans_baseline_iota0p94_delta_sticky',...
    'trans_baseline_iota0p94_mbar_sticky'};
macroShockNames = {'betas','zbar','abar','deltaBar','mbar'};
macroShockLabels = {'Average discount factor',...
    'Borrowing constraint',...
    'Aggregate productivity',...
    'Aggregate separation rate',...
    'Aggregate match efficiency'};

% Great Recession simulation
finName_ui = 'GR_iota0p94_ui';
finName_noui = 'GR_iota0p94_noui';
finName_nozlb = 'GR_iota0p94_nozlb';
finName_nozlbui = 'GR_iota0p94_nozlbui';

periods_mult = [3,8,11,13,19,21,23,24,28,33,45,47,57];
finName_ui_mult = {'trans_GR_iota0p94_time3_sticky',...
    'trans_GR_iota0p94_time8_sticky',...
    'trans_GR_iota0p94_time11_sticky',...
    'trans_GR_iota0p94_time13_sticky',...
    'trans_GR_iota0p94_time19_sticky',...
    'trans_GR_iota0p94_time21_sticky',...
    'trans_GR_iota0p94_time23_sticky',...
    'trans_GR_iota0p94_time24_sticky',...
    'trans_GR_iota0p94_time28_sticky',...
    'trans_GR_iota0p94_time33_sticky',...
    'trans_GR_iota0p94_time45_sticky',...
    'trans_GR_iota0p94_time47_sticky',...
    'trans_GR_iota0p94_time57_sticky'};
finName_noui_mult = {'trans_GR_iota0p94_time3_noui_sticky',...
    'trans_GR_iota0p94_time8_noui_sticky',...
    'trans_GR_iota0p94_time11_noui_sticky',...
    'trans_GR_iota0p94_time13_noui_sticky',...
    'trans_GR_iota0p94_time19_noui_sticky',...
    'trans_GR_iota0p94_time21_noui_sticky',...
    'trans_GR_iota0p94_time23_noui_sticky',...
    'trans_GR_iota0p94_time24_noui_sticky',...
    'trans_GR_iota0p94_time28_noui_sticky',...
    'trans_GR_iota0p94_time33_noui_sticky',...
    'trans_GR_iota0p94_time45_noui_sticky',...
    'trans_GR_iota0p94_time47_noui_sticky',...
    'trans_GR_iota0p94_time57_noui_sticky'};

finName_ui_iotaHi = 'GR_iota1p0_ui';
finName_ui_iotaLo = 'GR_iota0p88_ui';

finName_ui_iotaHi_mult = {'trans_GR_iota1p0_time3_sticky',...
    'trans_GR_iota1p0_time8_sticky',...
    'trans_GR_iota1p0_time11_sticky',...
    'trans_GR_iota1p0_time13_sticky',...
    'trans_GR_iota1p0_time19_sticky',...
    'trans_GR_iota1p0_time21_sticky',...
    'trans_GR_iota1p0_time23_sticky',...
    'trans_GR_iota1p0_time24_sticky',...
    'trans_GR_iota1p0_time28_sticky',...
    'trans_GR_iota1p0_time33_sticky',...
    'trans_GR_iota1p0_time45_sticky',...
    'trans_GR_iota1p0_time47_sticky',...
    'trans_GR_iota1p0_time57_sticky'};
finName_noui_iotaHi_mult = {'trans_GR_iota1p0_time3_noui_sticky',...
    'trans_GR_iota1p0_time8_noui_sticky',...
    'trans_GR_iota1p0_time11_noui_sticky',...
    'trans_GR_iota1p0_time13_noui_sticky',...
    'trans_GR_iota1p0_time19_noui_sticky',...
    'trans_GR_iota1p0_time21_noui_sticky',...
    'trans_GR_iota1p0_time23_noui_sticky',...
    'trans_GR_iota1p0_time24_noui_sticky',...
    'trans_GR_iota1p0_time28_noui_sticky',...
    'trans_GR_iota1p0_time33_noui_sticky',...
    'trans_GR_iota1p0_time45_noui_sticky',...
    'trans_GR_iota1p0_time47_noui_sticky',...
    'trans_GR_iota1p0_time57_noui_sticky'};

finName_ui_iotaLo_mult = {'trans_GR_iota0p88_time3_sticky',...
    'trans_GR_iota0p88_time8_sticky',...
    'trans_GR_iota0p88_time11_sticky',...
    'trans_GR_iota0p88_time13_sticky',...
    'trans_GR_iota0p88_time19_sticky',...
    'trans_GR_iota0p88_time21_sticky',...
    'trans_GR_iota0p88_time23_sticky',...
    'trans_GR_iota0p88_time24_sticky',...
    'trans_GR_iota0p88_time28_sticky',...
    'trans_GR_iota0p88_time33_sticky',...
    'trans_GR_iota0p88_time45_sticky',...
    'trans_GR_iota0p88_time47_sticky',...
    'trans_GR_iota0p88_time57_sticky'};
finName_noui_iotaLo_mult = {'trans_GR_iota0p88_time3_noui_sticky',...
    'trans_GR_iota0p88_time8_noui_sticky',...
    'trans_GR_iota0p88_time11_noui_sticky',...
    'trans_GR_iota0p88_time13_noui_sticky',...
    'trans_GR_iota0p88_time19_noui_sticky',...
    'trans_GR_iota0p88_time21_noui_sticky',...
    'trans_GR_iota0p88_time23_noui_sticky',...
    'trans_GR_iota0p88_time24_noui_sticky',...
    'trans_GR_iota0p88_time28_noui_sticky',...
    'trans_GR_iota0p88_time33_noui_sticky',...
    'trans_GR_iota0p88_time45_noui_sticky',...
    'trans_GR_iota0p88_time47_noui_sticky',...
    'trans_GR_iota0p88_time57_noui_sticky'};

finName_ui_seprate = 'GR_seprate_iota0p94_ui';
finName_noui_seprate = 'GR_seprate_iota0p94_noui';

finName_expiration = 'GR_iota0p94_expiration';
finName_noexpiration = 'GR_iota0p94_noexpiration';

%% PRODUCE ALL TABLES

%% Stationary RCE (sections 3, D.2, and D.3) 

for i=1:length(calNames)
    thisCal = load([outputName,'/',calNames{i}]);
    
    delete([outputName,'/figures_and_tables/',calTables{i}]); % clear any old output
    diary([outputName,'/figures_and_tables/',calTables{i}]);
   
    fprintf(strcat('Real interest rate','\t\t\t\t\t\t\t\t\t\t',sprintf('%.3f',r_target(i)),'\t',sprintf('%.3f',thisCal.r),'\t','znext_g=',sprintf('%.2f',thisCal.znext_g),'\n'));  
    fprintf(strcat('Wealth to mo income','\t\t\t\t\t\t\t\t\t\t',sprintf('%.1f',wealth_to_moincome_target(i)),'\t',sprintf('%.1f',thisCal.wealth_to_moincome),'\t','mean(betas)=',sprintf('%.5f',mean(thisCal.betas)),'\n'));
    if ~isnan(dwealth_ue_to_moincome_target(i))
        fprintf(strcat('(U - E) wealth to mo income','\t\t\t\t\t\t\t\t',sprintf('%.1f',dwealth_ue_to_moincome_target(i)),'\t',sprintf('%.1f',thisCal.dwealth_ue_to_moincome),'\t','eps_delta_beta=',sprintf('%.2f',thisCal.eps_delta_beta),'\n'));
    end
    fprintf(strcat('Frac agents with negative wealth','\t\t\t\t\t\t',sprintf('%.2f',share_wealth_lt0_target(i)),'\t',sprintf('%.2f',thisCal.share_wealth_lt0),'\t','zbar=',sprintf('%.2f',thisCal.zbar),'\n'));
    fprintf(strcat('Avg quarterly MPC out of $500 rebate, income <= 75k','\t\t',sprintf('%.2f',mpc_lt75k_target(i)),'\t',sprintf('%.2f',thisCal.mpc_lt75k),'\t','Delta=',sprintf('%.5f',thisCal.betas(length(thisCal.betas))-mean(thisCal.betas)),'\n'));

    fprintf(strcat('Frac unemp receiving UI','\t\t\t\t\t\t\t\t\t',sprintf('%.2f',frac_unemp_UI_target(i)),'\t',sprintf('%.2f',thisCal.frac_unemp_UI),'\t','zeta=',sprintf('%.1f',thisCal.zeta),'\n'));
    fprintf(strcat('Max UI relative to avg wage of UI recipients','\t\t\t',sprintf('%.2f',max_rr_target(i)),'\t',sprintf('%.2f',thisCal.max_rr),'\t','uimax=',sprintf('%.2f',thisCal.uimax),'\n'));
    fprintf(strcat('Avg income unemp HH with UI / pre-job-loss','\t\t\t\t',sprintf('%.2f',hhinc_ue_UI_target(i)),'\t',sprintf('%.2f',thisCal.hhinc_ue_UI),'\t','omega1=',sprintf('%.2f',thisCal.omega1),'\n'));
    fprintf(strcat('Avg income unemp HH without UI / pre-job-loss','\t\t\t',sprintf('%.2f',hhinc_ue_NoUI_target(i)),'\t',sprintf('%.2f',thisCal.hhinc_ue_NoUI),'\t','omega2=',sprintf('%.2f',thisCal.omega2),'\n'));

    fprintf(strcat('Unemp rate','\t\t\t\t\t\t\t\t\t\t\t\t',sprintf('%.3f',ur_target(i)),'\t',sprintf('%.3f',thisCal.ur),'\t','phi=',sprintf('%.3f',thisCal.phi),'\n'));
    fprintf(strcat('Share long-term unemp','\t\t\t\t\t\t\t\t\t',sprintf('%.2f',share_ur_lt_target(i)),'\t',sprintf('%.2f',thisCal.share_ur_lt),'\t','lambda=',sprintf('%.2f',thisCal.lambda),'\n'));   
    if ~isnan(weuCoeff_target(i))
        fprintf(strcat('EU on log wage','\t\t\t\t\t\t\t\t\t\t\t',sprintf('%.3f',weuCoeff_target(i)),'\t',sprintf('%.3f',thisCal.weuCoeff(2)),'\t','eps_delta_a=',sprintf('%.3f',thisCal.eps_delta_a),'\n'));
    end
    
    fprintf(strcat('Duration elast to benefit duration','\t\t\t\t\t\t',sprintf('%.2f',elast_dur_target(i)),'\t',sprintf('%.2f',thisCal.elast_dur),'\t','xi=',sprintf('%.1f',thisCal.xi),'\n'));
    fprintf(strcat('Conventional tightness','\t\t\t\t\t\t\t\t\t',sprintf('%.3f',conv_tightness_target(i)),'\t',sprintf('%.3f',thisCal.conv_tightness),'\t','mbar=',sprintf('%.2f',thisCal.mbar(1)),'\n'));
    fprintf(strcat('Frac of mo wage to recruit worker','\t\t\t\t\t\t',sprintf('%.3f',mowage_per_vac_target(i)),'\t',sprintf('%.3f',thisCal.mowage_per_vac),'\t','k=',sprintf('%.3f',thisCal.k),'\n'));

    diary off;
end

RCE = load([outputName,'/',RCE_baselineName]);
    
delete([outputName,'/figures_and_tables/Table2.txt']); % clear any old output
diary([outputName,'/figures_and_tables/Table2.txt']);
    
fprintf(strcat('Avg annual MPC out of mo income rebate, unemp-emp','\t\t',sprintf('%.2f',0.25),'\t',sprintf('%.2f',RCE.mpc_uAAvg-RCE.mpc_eAAvg),'\n'));
fprintf(strcat('Ganong-Noel 19 dC/dY at exhaustion','\t\t\t\t\t\t',sprintf('%.2f',0.20),'\t',sprintf('%.2f',RCE.gn17dCdY2),'\n'));

fprintf(strcat('EU on wealth/mo income','\t\t\t\t\t\t\t\t\t',sprintf('%.4f',-0.0002),'\t',sprintf('%.4f',RCE.wealtheuCoeff(2)),'\n'));
fprintf(strcat('Ganong-Noel 19 cons during receipt','\t\t\t\t\t\t',sprintf('%.2f',0.91),'\t',sprintf('%.2f',RCE.gn17receipt),'\n'));
fprintf(strcat('Ganong-Noel 19 cons after exhaustion','\t\t\t\t\t',sprintf('%.2f',0.80),'\t',sprintf('%.2f',RCE.gn17exhaust),'\n'));

fprintf(strcat('Aggregate consumption share, wealth Q5','\t\t\t\t\t',sprintf('%.3f',0.372),'\t',sprintf('%.3f',RCE.aggc_ZQ5/RCE.aggc),'\n'));  
fprintf(strcat('Aggregate consumption share, wealth Q4','\t\t\t\t\t',sprintf('%.3f',0.224),'\t',sprintf('%.3f',RCE.aggc_ZQ4/RCE.aggc),'\n'));  
fprintf(strcat('Aggregate consumption share, wealth Q3','\t\t\t\t\t',sprintf('%.3f',0.168),'\t',sprintf('%.3f',RCE.aggc_ZQ3/RCE.aggc),'\n'));
fprintf(strcat('Aggregate consumption share, wealth Q2','\t\t\t\t\t',sprintf('%.3f',0.124),'\t',sprintf('%.3f',RCE.aggc_ZQ2/RCE.aggc),'\n'));
fprintf(strcat('Aggregate consumption share, wealth Q1','\t\t\t\t\t',sprintf('%.3f',0.113),'\t',sprintf('%.3f',RCE.aggc_ZQ1/RCE.aggc),'\n'));

diary off;
  
delete([outputName,'/figures_and_tables/Table3.txt']); % clear any old output
diary([outputName,'/figures_and_tables/Table3.txt']);

fprintf(strcat('Avg quarterly MPC out of $500 rebate','\t\t\t\t\t',sprintf('%.2f',RCE.mpc_QAvg),'\n'));    
fprintf(strcat('Avg quarterly MPC out of $500 rebate, wealth Q5','\t\t\t',sprintf('%.2f',RCE.mpc_ZQ5),'\n'));            
fprintf(strcat('Avg quarterly MPC out of $500 rebate, wealth Q4','\t\t\t',sprintf('%.2f',RCE.mpc_ZQ4),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, wealth Q3','\t\t\t',sprintf('%.2f',RCE.mpc_ZQ3),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, wealth Q2','\t\t\t',sprintf('%.2f',RCE.mpc_ZQ2),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, wealth Q1','\t\t\t',sprintf('%.2f',RCE.mpc_ZQ1),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, emp','\t\t\t\t',sprintf('%.2f',RCE.mpc_eQAvg),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, unemp','\t\t\t\t',sprintf('%.2f',RCE.mpc_uQAvg),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, ST unemp','\t\t\t',sprintf('%.2f',RCE.mpc_uSTQAvg),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, MT unemp','\t\t\t',sprintf('%.2f',RCE.mpc_uMTQAvg),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, LT unemp','\t\t\t',sprintf('%.2f',RCE.mpc_uLTQAvg),'\n'));

diary off;

%% Baseline effects of UI starting from stationary RCE (section 4)

[avgDUnemp_flex,avgDR_flex,UImult_flex,~,UImult_pe] = ...
    computeUIMultiplier(RCE.b,[outputName,'/',initName],[outputName,'/',peName],[outputName,'/',finName_flex]);
[avgDUnemp_sticky,avgDR_sticky,UImult_sticky,~,~] = ...
    computeUIMultiplier(RCE.b,[outputName,'/',initName],[],[outputName,'/',finName_sticky]);
[avgDUnemp_m,avgDR_m,UImult_m,~,~] = ...
    computeUIMultiplier(RCE.b,[outputName,'/',initName],[],[outputName,'/',finName_m]);
[avgDUnemp_zlb,avgDR_zlb,UImult_zlb,~,~] = ...
    computeUIMultiplier(RCE.b,[outputName,'/',initName],[],[outputName,'/',finName_zlb]);

delete([outputName,'/figures_and_tables/Table4.txt']); % clear any old output
diary([outputName,'/figures_and_tables/Table4.txt']);
    
fprintf(strcat('Output multiplier','\t\t\t\t',...
    sprintf('%.1f',UImult_flex),'\t',sprintf('%.1f',UImult_sticky),'\t\t',...
    sprintf('%.1f',UImult_m),'\t\t',sprintf('%.1f',UImult_zlb),'\n'));
fprintf(strcat('Avg change in unemp rate','\t\t',...
    sprintf('%.2f',100*avgDUnemp_flex),'\t',sprintf('%.2f',100*avgDUnemp_sticky),'\t',...
    sprintf('%.2f',100*avgDUnemp_m),'\t',sprintf('%.2f',100*avgDUnemp_zlb),'\n'));

diary off;

%% Sensitivity analysis starting from stationary RCE (sections 4, D.1, D.4, and D.5)

RCE_noEps = load([outputName,'/',RCE_Name_noEps]);
[avgDUnemp_zlb_noEps,~,UImult_zlb_noEps,~,~] = ...
    computeUIMultiplier(RCE_noEps.b,[outputName,'/',initName_noEps],[],[outputName,'/',finName_zlb_noEps]);

delete([outputName,'/figures_and_tables/Table5.txt']); % clear any old output
diary([outputName,'/figures_and_tables/Table5.txt']);
    
fprintf(strcat('Avg quarterly MPC out of $500 rebate, emp','\t\t\t\t',sprintf('%.2f',RCE.mpc_eQAvg),'\t',sprintf('%.2f',RCE_noEps.mpc_eQAvg),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, ST unemp','\t\t\t',sprintf('%.2f',RCE.mpc_uSTQAvg),'\t',sprintf('%.2f',RCE_noEps.mpc_uSTQAvg),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, MT unemp','\t\t\t',sprintf('%.2f',RCE.mpc_uMTQAvg),'\t',sprintf('%.2f',RCE_noEps.mpc_uMTQAvg),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, LT unemp','\t\t\t',sprintf('%.2f',RCE.mpc_uLTQAvg),'\t',sprintf('%.2f',RCE_noEps.mpc_uLTQAvg),'\n'));
fprintf(strcat('Output multiplier','\t\t\t\t\t\t\t\t\t\t',sprintf('%.1f',UImult_zlb),'\t\t',sprintf('%.1f',UImult_zlb_noEps),'\n'));
fprintf(strcat('Avg change in unemp rate','\t\t\t\t\t\t\t\t',sprintf('%.2f',100*avgDUnemp_zlb),'\t',sprintf('%.2f',100*avgDUnemp_zlb_noEps),'\n'));

diary off;

RCE_liquid = load([outputName,'/',RCE_Name_liquid]);
[avgDUnemp_zlb_liquid,~,UImult_zlb_liquid,~,~] = ...
    computeUIMultiplier(RCE_liquid.b,[outputName,'/',initName_liquid],[],[outputName,'/',finName_zlb_liquid]);

delete([outputName,'/figures_and_tables/TableA8.txt']); % clear any old output
diary([outputName,'/figures_and_tables/TableA8.txt']);
    
fprintf(strcat('Avg quarterly MPC out of $500 rebate, emp','\t\t\t\t',sprintf('%.2f',RCE.mpc_eQAvg),'\t',sprintf('%.2f',RCE_liquid.mpc_eQAvg),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, ST unemp','\t\t\t',sprintf('%.2f',RCE.mpc_uSTQAvg),'\t',sprintf('%.2f',RCE_liquid.mpc_uSTQAvg),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, MT unemp','\t\t\t',sprintf('%.2f',RCE.mpc_uMTQAvg),'\t',sprintf('%.2f',RCE_liquid.mpc_uMTQAvg),'\n'));
fprintf(strcat('Avg quarterly MPC out of $500 rebate, LT unemp','\t\t\t',sprintf('%.2f',RCE.mpc_uLTQAvg),'\t',sprintf('%.2f',RCE_liquid.mpc_uLTQAvg),'\n'));
fprintf(strcat('Output multiplier','\t\t\t\t\t\t\t\t\t\t',sprintf('%.1f',UImult_zlb),'\t\t',sprintf('%.1f',UImult_zlb_liquid),'\n'));
fprintf(strcat('Avg change in unemp rate','\t\t\t\t\t\t\t\t',sprintf('%.2f',100*avgDUnemp_zlb),'\t',sprintf('%.2f',100*avgDUnemp_zlb_liquid),'\n'));

diary off;

[avgDUnemp_flex_iota0p9,avgDR_flex_iota0p9,UImult_flex_iota0p9,~,~] = ...
    computeUIMultiplier(RCE.b,[outputName,'/',initName],[],[outputName,'/',finName_flex_iota0p9]);
[avgDUnemp_sticky_iota0p9,avgDR_sticky_iota0p9,UImult_sticky_iota0p9,~,~] = ...
    computeUIMultiplier(RCE.b,[outputName,'/',initName],[],[outputName,'/',finName_sticky_iota0p9]);
[avgDUnemp_zlb_iota0p9,avgDR_zlb_iota0p9,UImult_zlb_iota0p9,~,~] = ...
    computeUIMultiplier(RCE.b,[outputName,'/',initName],[],[outputName,'/',finName_zlb_iota0p9]);

RCE_durElast = load([outputName,'/',RCE_Name_durElast]);
[avgDUnemp_flex_durElast,avgDR_flex_durElast,UImult_flex_durElast,~,~] = ...
    computeUIMultiplier(RCE_durElast.b,[outputName,'/',initName_durElast],[],[outputName,'/',finName_flex_durElast]);
[avgDUnemp_sticky_durElast,avgDR_sticky_durElast,UImult_sticky_durElast,~,~] = ...
    computeUIMultiplier(RCE_durElast.b,[outputName,'/',initName_durElast],[],[outputName,'/',finName_sticky_durElast]);
[avgDUnemp_zlb_durElast,avgDR_zlb_durElast,UImult_zlb_durElast,~,~] = ...
    computeUIMultiplier(RCE_durElast.b,[outputName,'/',initName_durElast],[],[outputName,'/',finName_zlb_durElast]);

delete([outputName,'/figures_and_tables/Table6.txt']); % clear any old output
diary([outputName,'/figures_and_tables/Table6.txt']);
 
fprintf(strcat('Avg change in real rate (ann.)','\t\t\t',sprintf('%.2f',100*avgDR_flex),'\t',sprintf('%.2f',100*avgDR_flex_iota0p9),'\t',sprintf('%.2f',100*avgDR_flex_durElast),'\n'));
fprintf(strcat('Output multiplier','\t\t\t\t\t\t',sprintf('%.1f',UImult_flex),'\t',sprintf('%.1f',UImult_flex_iota0p9),'\t',sprintf('%.1f',UImult_flex_durElast),'\n'));
fprintf(strcat('Avg change in unemp rate','\t\t\t\t',sprintf('%.2f',100*avgDUnemp_flex),'\t',sprintf('%.2f',100*avgDUnemp_flex_iota0p9),'\t',sprintf('%.2f',100*avgDUnemp_flex_durElast),'\n'));

fprintf(strcat('Avg change in real rate (ann.)','\t\t\t',sprintf('%.2f',100*avgDR_sticky),'\t',sprintf('%.2f',100*avgDR_sticky_iota0p9),'\t',sprintf('%.2f',100*avgDR_sticky_durElast),'\n'));
fprintf(strcat('Output multiplier','\t\t\t\t\t\t',sprintf('%.1f',UImult_sticky),'\t\t',sprintf('%.1f',UImult_sticky_iota0p9),'\t\t',sprintf('%.1f',UImult_sticky_durElast),'\n'));
fprintf(strcat('Avg change in unemp rate','\t\t\t\t',sprintf('%.2f',100*avgDUnemp_sticky),'\t',sprintf('%.2f',100*avgDUnemp_sticky_iota0p9),'\t',sprintf('%.2f',100*avgDUnemp_sticky_durElast),'\n'));

fprintf(strcat('Avg change in real rate (ann.)','\t\t\t',sprintf('%.2f',100*avgDR_zlb),'\t',sprintf('%.2f',100*avgDR_zlb_iota0p9),'\t',sprintf('%.2f',100*avgDR_zlb_durElast),'\n'));
fprintf(strcat('Output multiplier','\t\t\t\t\t\t',sprintf('%.1f',UImult_zlb),'\t\t',sprintf('%.1f',UImult_zlb_iota0p9),'\t\t',sprintf('%.1f',UImult_zlb_durElast),'\n'));
fprintf(strcat('Avg change in unemp rate','\t\t\t\t',sprintf('%.2f',100*avgDUnemp_zlb),'\t',sprintf('%.2f',100*avgDUnemp_zlb_iota0p9),'\t',sprintf('%.2f',100*avgDUnemp_zlb_durElast),'\n'));

diary off;

[avgDUnemp_zlb_zeta,~,UImult_zlb_zeta,~,~] = ...
    computeUIMultiplier(RCE.b,[outputName,'/',initName_zeta],[],[outputName,'/',finName_zlb_zeta]);
[avgDUnemp_zlb_debt,~,UImult_zlb_debt,~,~] = ...
    computeUIMultiplier(RCE.b,[outputName,'/',initName],[],[outputName,'/',finName_zlb_debt]);
[avgDUnemp_zlb_rr,~,UImult_zlb_rr,~,~] = ...
    computeUIMultiplier(RCE.b,[outputName,'/',initName],[],[outputName,'/',finName_zlb_rr]);

delete([outputName,'/figures_and_tables/TableA9.txt']); % clear any old output
diary([outputName,'/figures_and_tables/TableA9.txt']);
    
fprintf(strcat('Output multiplier','\t\t\t\t\t\t',sprintf('%.1f',UImult_zlb),'\t\t',sprintf('%.1f',UImult_zlb_zeta),'\t\t',sprintf('%.1f',UImult_zlb_debt),'\t\t',sprintf('%.1f',UImult_zlb_rr),'\n'));
fprintf(strcat('Avg change in unemp rate','\t\t\t\t',sprintf('%.2f',100*avgDUnemp_zlb),'\t',sprintf('%.2f',100*avgDUnemp_zlb_noEps),'\t',sprintf('%.2f',100*avgDUnemp_zlb_debt),'\t',sprintf('%.2f',100*avgDUnemp_zlb_rr),'\n'));

diary off;

[~,~,Gmult_sticky,~,~] = ...
    computeGMultiplier([outputName,'/',initName],[],[outputName,'/',finName_sticky_fiscal]);
[~,~,Gmult_zlb,~,~] = ...
    computeGMultiplier([outputName,'/',initName],[],[outputName,'/',finName_zlb_fiscal]);
[~,~,Gmult_debt_sticky,~,~] = ...
    computeGMultiplier([outputName,'/',initName],[],[outputName,'/',finName_sticky_fiscal_debt]);
[~,~,Gmult_debt_zlb,~,~] = ...
    computeGMultiplier([outputName,'/',initName],[],[outputName,'/',finName_zlb_fiscal_debt]);

delete([outputName,'/figures_and_tables/TableA10.txt']); % clear any old output
diary([outputName,'/figures_and_tables/TableA10.txt']);
    
fprintf(strcat('Budget-balanced','\t\t\t',sprintf('%.1f',Gmult_sticky),'\t\t',sprintf('%.1f',Gmult_zlb),'\n'));
fprintf(strcat('Deficit-financed','\t\t',sprintf('%.1f',Gmult_debt_sticky),'\t\t',sprintf('%.1f',Gmult_debt_zlb),'\n'));

diary off;

%% Great Recession simulation (section 6)

avgDUunemp_ui = nan(1,length(periods_mult));
UImult_ui = nan(1,length(periods_mult));
RCE_dbarbar22 = load([outputName,'/',RCE_baseline_dbarbar22Name]);
for i=1:length(periods_mult) 
    [thisAvgDUunemp_ui,~,thisUImult_ui,~,~] = ...
        computeUIMultiplier(RCE_dbarbar22.b,[outputName,'/',finName_noui_mult{i}],[],[outputName,'/',finName_ui_mult{i}]);
    avgDUunemp_ui(i) = thisAvgDUunemp_ui;
    UImult_ui(i) = thisUImult_ui;
end

delete([outputName,'/figures_and_tables/Table8.txt']); % clear any old output
diary([outputName,'/figures_and_tables/Table8.txt']);

for i=1:length(periods_mult)
    fprintf(strcat(num2str(periods_mult(i)),'\t\t\t',sprintf('%.1f',UImult_ui(i)),'\t\t',sprintf('%.2f',100*avgDUunemp_ui(i)),'\n'));
end

diary off;

UImult_iotaHi_ui = nan(1,length(periods_mult));
UImult_iotaLo_ui = nan(1,length(periods_mult));
for i=1:length(periods_mult) 
    [~,~,thisUImult_iotaHi_ui,~,~] = ...
        computeUIMultiplier(RCE_dbarbar22.b,[outputName,'/',finName_noui_iotaHi_mult{i}],[],[outputName,'/',finName_ui_iotaHi_mult{i}]);
    UImult_iotaHi_ui(i) = thisUImult_iotaHi_ui;
    
    [~,~,thisUImult_iotaLo_ui,~,~] = ...
        computeUIMultiplier(RCE_dbarbar22.b,[outputName,'/',finName_noui_iotaLo_mult{i}],[],[outputName,'/',finName_ui_iotaLo_mult{i}]);
    UImult_iotaLo_ui(i) = thisUImult_iotaLo_ui;    
end

delete([outputName,'/figures_and_tables/Table9.txt']); % clear any old output
diary([outputName,'/figures_and_tables/Table9.txt']);

for i=1:length(periods_mult)
    fprintf(strcat(num2str(periods_mult(i)),'\t\t\t',sprintf('%.1f',UImult_iotaLo_ui(i)),'\t\t',sprintf('%.1f',UImult_ui(i)),'\t\t',sprintf('%.1f',UImult_iotaHi_ui(i)),'\n'));
end

diary off;

%% PRODUCE ALL FIGURES

%% Set plotting format
set(0,'defaultaxescolororder',[0 0 0.5]);
set(0,'defaultaxeslinestyleorder',{'-','--','-.',':'});
set(0,'defaultLineLineWidth',2);
set(0,'defaultAxesLineWidth',0.5);
Interpreter = 'latex';
%set(0,'DefaultAxesFontName','CMU Serif Roman');

%% Stationary RCE (sections 3, D.2, and D.3)

for identifyParamInd=1:length(identifyParamNames)
    figure;    
    hold on;
    thisParamRange = identifyParams{identifyParamInd};
    for i=1:5
        try
            thisStruct = load([outputName,'\',identifyFileNames{identifyParamInd},'_',num2str(i),'.mat'],identifyTargetNames{identifyParamInd});
            thisMoment = extractfield(thisStruct,identifyTargetNames{identifyParamInd});
            if find(strcmp(identifyTargetNames{identifyParamInd},{'weuCoeff','wealtheuCoeff'}),1)
                thisMoment = thisMoment(2);
            end            
            if thisParamRange(i) ~= identifyParamsChosen{identifyParamInd}
                plot(thisParamRange(i),thisMoment,'o','MarkerSize',6);
            else
                plot(thisParamRange(i),thisMoment,'o','MarkerFaceColor',[0 0 0.5],'MarkerSize',6);
            end
            clearvars thisStruct thisMoment;
        catch
        end
    end
    hline = refline(0,identifyTargets{identifyParamInd});
    hline.LineWidth = 0.5;  
    hold off;
    title(identifyTargetLabels{identifyParamInd},'Interpreter',Interpreter);
    fig=gcf;  
    yl = ylim;    
    axis([min(thisParamRange),max(thisParamRange),-inf,inf]);
    sortedRange = sort(thisParamRange);
    set(gca,'xtick',sortedRange);
    xl = xlim;
    text(max(thisParamRange)-(xl(2)-xl(1))/80,identifyTargets{identifyParamInd}+(yl(2)-yl(1))/30,...
        strcat('target=',num2str(identifyTargets{identifyParamInd})),...
        'HorizontalAlignment','right','Interpreter',Interpreter);
    box on;
    curtick = get(gca, 'XTick');
    set(gca, 'XTickLabel', cellstr(num2str(curtick(:),6)));
    set(gca,'TickLabelInterpreter',Interpreter);
    set(findall(fig,'-property','FontSize'),'FontSize',11);     
    set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
    xlabel(identifyParamLabels{identifyParamInd},'Interpreter','latex',...
        'FontSize',11);
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
    print(gcf, strcat([outputName,'/figures_and_tables/Figure1_',num2str(identifyParamInd)]), '-depsc', '-r600','-loose');
end

avgByInd = @(values,dist,inds) sum(values(:,inds).*repmat(dist(inds),size(values,1),1),2)/...
    sum(dist(inds));

lb = 1;
rbVal = 10;
rb = find(RCE.zStat >= rbVal,1);

% Compute distribution of unemployed at each value of
%  (uBeta,uType,uUI,uEndow), but *summed* over uDur.  Used to aggregate
%  policy functions and distributions below to depict how these change
%  with unemployment duration alone.
phi_u_tildeoverdur = zeros(1,RCE.uStates);
phi_u_overdur = zeros(1,RCE.uStates);
for i_u = 1:(RCE.uStates/(RCE.dbarbar+1))
    phi_u_tildeoverdur(1+(i_u-1)*(RCE.dbarbar+1):RCE.dbarbar+1+(i_u-1)*(RCE.dbarbar+1)) = sum(sum(RCE.phi_u_tilde(:,1+(i_u-1)*(RCE.dbarbar+1):RCE.dbarbar+1+(i_u-1)*(RCE.dbarbar+1))))/(RCE.dbarbar+1);
    phi_u_overdur(1+(i_u-1)*(RCE.dbarbar+1):RCE.dbarbar+1+(i_u-1)*(RCE.dbarbar+1)) = sum(sum(RCE.phi_u(:,1+(i_u-1)*(RCE.dbarbar+1):RCE.dbarbar+1+(i_u-1)*(RCE.dbarbar+1))))/(RCE.dbarbar+1);
end

figure;
hax=axes; 
plot(RCE.zStat(lb:rb),avgByInd(RCE.c_eStat(lb:rb,:),sum(RCE.phi_e,1),logical(ones(1,RCE.eStates))));
hold on;
plot(RCE.zStat(lb:rb),avgByInd(RCE.c_uStat(lb:rb,:),phi_u_overdur,RCE.uDur==1));
plot(RCE.zStat(lb:rb),avgByInd(RCE.c_uStat(lb:rb,:),phi_u_overdur,RCE.uDur==4));
plot(RCE.zStat(lb:rb),avgByInd(RCE.c_uStat(lb:rb,:),phi_u_overdur,RCE.uDur==7));
line([RCE.zbar RCE.zbar],get(hax,'YLim'),'LineStyle',':','LineWidth',0.5,'Color','k');
hold off;
legend('e','u(d=0)','u(d=3)','u(d=6)','Location','southeast');        
xlabel('Assets','Interpreter',Interpreter);
title('\textbf{Consumption}','Interpreter',Interpreter);
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11);
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
axis([floor(RCE.zbar-0.5),floor(rbVal),-inf,inf]);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure2_1']), '-depsc', '-r600','-loose');

figure;
hax=axes; 
phi_eAvg = sum(RCE.phi_e,2);
phi_uAvgDur0 = sum(RCE.phi_u(:,RCE.uDur==1),2)./sum(sum(RCE.phi_u(:,RCE.uDur==1)));
phi_uAvgDur3 = sum(RCE.phi_u(:,RCE.uDur==4),2)./sum(sum(RCE.phi_u(:,RCE.uDur==4)));
phi_uAvgDur6 = sum(RCE.phi_u(:,RCE.uDur==7),2)./sum(sum(RCE.phi_u(:,RCE.uDur==7)));
plot(RCE.zStat(lb:rb),phi_eAvg(lb:rb));
hold on;
plot(RCE.zStat(lb:rb),phi_uAvgDur0(lb:rb));
plot(RCE.zStat(lb:rb),phi_uAvgDur3(lb:rb));
plot(RCE.zStat(lb:rb),phi_uAvgDur6(lb:rb));
line([RCE.zbar RCE.zbar],get(hax,'YLim'),'LineStyle',':','LineWidth',0.5,'Color','k');
hold off;       
xlabel('Assets','Interpreter',Interpreter);
title('\textbf{Invariant marginal distributions}','Interpreter',Interpreter);
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11);
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
axis([floor(RCE.zbar-0.5),floor(rbVal),-inf,0.3]); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure2_2']), '-depsc', '-r600','-loose');

%% Baseline effects of UI starting from stationary RCE (sections 4 and D.1)

init = load([outputName,'\',initName]);
pe = load([outputName,'\',peName]);
fin_flex = load([outputName,'\',finName_flex]);
fin_sticky = load([outputName,'\',finName_sticky]);
fin_m = load([outputName,'\',finName_m]);
fin_zlb = load([outputName,'\',finName_zlb]);
T = 24;

cStatSS = RCE.cStat;
mSS = RCE.m;
p_eSS = RCE.p_e;
phi_eSS = RCE.phi_e;
phi_uSS = RCE.phi_u;
uDurSS = RCE.uDur;
dbarSS = RCE.dbar;
wSS = RCE.w;
aSS = RCE.a;
zbarSS = RCE.zbar;
betasSS = RCE.betas;
thetaSS = RCE.theta;
sbarSS = RCE.sbar;
p_e_tildeSS = RCE.p_e_tilde;

cTransInit = init.cTrans;
mInit = init.m;
MTransInit = init.MTrans;
p_eInit = init.p_e;
p_e_tildeInit = init.p_e_tilde;
p_uInit = init.p_u;
thetaInit = init.theta;
sbarTransInit = init.sbarTrans;
muInit = init.mu;
PiTransInit = init.PiTrans;
cpiInit = cumprod(1+PiTransInit)-1;

cTransPE = pe.cTrans;

cTransFinFlex = fin_flex.cTrans;
mFinFlex = fin_flex.m;
p_eFlex = fin_flex.p_e;
p_e_tildeFlex = fin_flex.p_e_tilde;
p_uFlex = fin_flex.p_u;
thetaFlex = fin_flex.theta;
sbarTransFlex = fin_flex.sbarTrans;
muFlex = fin_flex.mu;
PiTransFlex = fin_flex.PiTrans;

cTransFinSticky = fin_sticky.cTrans;
mFinSticky = fin_sticky.m;
MTransFinSticky = fin_sticky.MTrans;
p_eSticky = fin_sticky.p_e;
p_e_tildeSticky = fin_sticky.p_e_tilde;
p_uSticky = fin_sticky.p_u;
thetaSticky = fin_sticky.theta;
sbarTransSticky = fin_sticky.sbarTrans;
muSticky = fin_sticky.mu;
PiTransSticky = fin_sticky.PiTrans;
cpiSticky = cumprod(1+PiTransSticky)-1;

cTransFinM = fin_m.cTrans;
mFinM = fin_m.m;
MTransFinM = fin_m.MTrans;
p_eM = fin_m.p_e;
p_e_tildeM = fin_m.p_e_tilde;
p_uM = fin_m.p_u;
thetaM = fin_m.theta;
sbarTransM = fin_m.sbarTrans;
muM = fin_m.mu;
PiTransM = fin_m.PiTrans;
cpiM = cumprod(1+PiTransM)-1;

cTransFinZLB = fin_zlb.cTrans;
mFinZLB = fin_zlb.m;
MTransFinZLB = fin_zlb.MTrans;
p_eZLB = fin_zlb.p_e;
p_e_tildeZLB = fin_zlb.p_e_tilde;
p_uZLB = fin_zlb.p_u;
thetaZLB = fin_zlb.theta;
sbarTransZLB = fin_zlb.sbarTrans;
muZLB = fin_zlb.mu;
PiTransZLB = fin_zlb.PiTrans;
cpiZLB = cumprod(1+PiTransZLB)-1;

figure;
plot(-4:1:(T-1),[zeros(1,4),100*(cTransPE(1:T) - cTransInit(1:T))/cStatSS],'LineWidth',1);
hold on
plot(-4:1:(T-1),[zeros(1,4),100*(cTransFinFlex(1:T) - cTransInit(1:T))/cStatSS],'-');
plot(-4:1:(T-1),[zeros(1,4),100*(cTransFinSticky(1:T) - cTransInit(1:T))/cStatSS],'--');
plot(-4:1:(T-1),[zeros(1,4),100*(cTransFinM(1:T) - cTransInit(1:T))/cStatSS],'-.');
plot(-4:1:(T-1),[zeros(1,4),100*(cTransFinZLB(1:T) - cTransInit(1:T))/cStatSS],':');
hold off;
axis([-4 T-1 -0.1 0.07]);
title(['\textbf{Aggregate consumption}'],'Interpreter',Interpreter);
ylabel('\% of SS consumption','Interpreter',Interpreter);
xlabel('Month','Interpreter',Interpreter);
legend('PE','GE, flex prices','GE, sticky prices','GE, sticky prices + fixed $r$','GE, sticky prices + fixed $i$','Location','Southeast');
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure3_1']), '-depsc', '-r600','-loose');

figure;
plot(-4:1:(T-1),[zeros(1,4),(1./mFinFlex(1:T)-1)*1200 - (1./mInit(1:T)-1)*1200]);
hold on;
plot(-4:1:(T-1),[zeros(1,4),(1./mFinSticky(1:T)-1)*1200 - (1./mInit(1:T)-1)*1200]);
plot(-4:1:(T-1),[zeros(1,4),(1./mFinM(1:T)-1)*1200 - (1./mInit(1:T)-1)*1200]);
plot(-4:1:(T-1),[zeros(1,4),(1./mFinZLB(1:T)-1)*1200 - (1./mInit(1:T)-1)*1200]);
hold off;
axis([-4 T-1 -0.015 0.05]);
title(['\textbf{Real interest rate (ann.)}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
xlabel('Month','Interpreter',Interpreter);
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure3_2']), '-depsc', '-r600','-loose');

figure;
plot(-4:1:(T-1),[zeros(1,4),100*(1-p_eFlex(1:T)) - 100*(1-p_eInit(1:T))]);
hold on;
plot(-4:1:(T-1),[zeros(1,4),100*(1-p_eSticky(1:T)) - 100*(1-p_eInit(1:T))]);
plot(-4:1:(T-1),[zeros(1,4),100*(1-p_eM(1:T)) - 100*(1-p_eInit(1:T))]);
plot(-4:1:(T-1),[zeros(1,4),100*(1-p_eZLB(1:T)) - 100*(1-p_eInit(1:T))]);
hold off;
axis([-4 T-1 -inf inf]);
title(['\textbf{Unemployment rate}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
xlabel('Month','Interpreter',Interpreter);
legend('Flex prices','Sticky prices','Sticky prices + fixed $r$','Sticky prices + fixed $i$','Location','Southeast');
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA2_1']), '-depsc', '-r600','-loose');

vacSS = thetaSS * sbarSS;
vacInit = thetaInit .* sbarTransInit;
vacFlex = thetaFlex .* sbarTransFlex;
vacSticky = thetaSticky .* sbarTransSticky;
vacM = thetaM .* sbarTransM;
vacZLB = thetaZLB .* sbarTransZLB;
figure;
plot(-4:1:(T-1),[zeros(1,4),100*(vacFlex(1:T) - vacInit(1:T))/vacSS]);
hold on;
plot(-4:1:(T-1),[zeros(1,4),100*(vacSticky(1:T) - vacInit(1:T))/vacSS]);
plot(-4:1:(T-1),[zeros(1,4),100*(vacM(1:T) - vacInit(1:T))/vacSS]);
plot(-4:1:(T-1),[zeros(1,4),100*(vacZLB(1:T) - vacInit(1:T))/vacSS]);
hold off;
axis([-4 T-1 -inf inf]);
title(['\textbf{Vacancies}'],'Interpreter',Interpreter);
ylabel('\% of SS','Interpreter',Interpreter);
xlabel('Month','Interpreter',Interpreter);
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA2_3']), '-depsc', '-r600','-loose');     

avgSbarSS = sbarSS / (1-p_e_tildeSS);
avgSbarInit = sbarTransInit ./ (1-p_e_tildeInit);
avgSbarFlex = sbarTransFlex ./ (1-p_e_tildeFlex);
avgSbarSticky = sbarTransSticky ./ (1-p_e_tildeSticky);
avgSbarM = sbarTransM ./ (1-p_e_tildeM);
avgSbarZLB = sbarTransZLB ./ (1-p_e_tildeZLB);
figure;
plot(-4:1:(T-1),[zeros(1,4),100*(avgSbarFlex(1:T) - avgSbarInit(1:T))/avgSbarSS]);
hold on;
plot(-4:1:(T-1),[zeros(1,4),100*(avgSbarSticky(1:T) - avgSbarInit(1:T))/avgSbarSS]);
plot(-4:1:(T-1),[zeros(1,4),100*(avgSbarM(1:T) - avgSbarInit(1:T))/avgSbarSS]);
plot(-4:1:(T-1),[zeros(1,4),100*(avgSbarZLB(1:T) - avgSbarInit(1:T))/avgSbarSS]);
hold off;
axis([-4 T-1 -inf inf]);
title(['\textbf{Average search effort}'],'Interpreter',Interpreter);
ylabel('\% of SS','Interpreter',Interpreter);
xlabel('Month','Interpreter',Interpreter);
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA2_4']), '-depsc', '-r600','-loose');     

figure;
plot(-4:1:(T-1),[zeros(1,4),100*(thetaFlex(1:T) - thetaInit(1:T))/thetaSS]);
hold on;
plot(-4:1:(T-1),[zeros(1,4),100*(thetaSticky(1:T) - thetaInit(1:T))/thetaSS]);
plot(-4:1:(T-1),[zeros(1,4),100*(thetaM(1:T) - thetaInit(1:T))/thetaSS]);
plot(-4:1:(T-1),[zeros(1,4),100*(thetaZLB(1:T) - thetaInit(1:T))/thetaSS]);
hold off;
axis([-4 T-1 -inf inf]);
title(['\textbf{Labor market tightness}'],'Interpreter',Interpreter);
ylabel('\% of SS','Interpreter',Interpreter);
xlabel('Month','Interpreter',Interpreter);
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA2_2']), '-depsc', '-r600','-loose');     

figure;
plot(-4:1:(T-1),[zeros(1,4),(1./MTransFinSticky(1:T)-1)*1200 - (1./MTransInit(1:T)-1)*1200],'--');
hold on;
plot(-4:1:(T-1),[zeros(1,4),(1./MTransFinM(1:T)-1)*1200 - (1./MTransInit(1:T)-1)*1200],'-.');
plot(-4:1:(T-1),[zeros(1,4),(1./MTransFinZLB(1:T)-1)*1200 - (1./MTransInit(1:T)-1)*1200],':');
hold off;
axis([-4 T-1 -inf inf]);
title(['\textbf{Nominal interest rate (ann.)}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
xlabel('Month','Interpreter',Interpreter);
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA2_5']), '-depsc', '-r600','-loose');   

figure;
plot(-4:1:(T-1),[zeros(1,4),100*(cpiSticky(1:T) - cpiInit(1:T))],'--');
hold on;
plot(-4:1:(T-1),[zeros(1,4),100*(cpiM(1:T) - cpiInit(1:T))],'-.');
plot(-4:1:(T-1),[zeros(1,4),100*(cpiZLB(1:T) - cpiInit(1:T))],':');
axis([-4 T-1 -inf inf]);
title(['\textbf{Final good prices}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
xlabel('Month','Interpreter',Interpreter);
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA2_6']), '-depsc', '-r600','-loose');   

%% Sensitivity analysis starting from stationary RCE (sections 4, D.4, and D.5)

durMultipliersHoriz = nan(1,length(durLabels));
for i=1:length(durLabels)
    [~,~,thisUImult,~,~] = ...
        computeUIMultiplier(RCE.b,[outputName,'/',initName],[],[outputName,'/',durFinNames{i}]);    
    durMultipliersHoriz(i) = thisUImult;
end

figure;
hold on;
for i=1:length(durLabels)
    if durLabels(i) ~= 12
        plot(durLabels(i),durMultipliersHoriz(i),'o','MarkerSize',6);
    else
        plot(durLabels(i),durMultipliersHoriz(i),'o','MarkerFaceColor',[0 0 0.5],'MarkerSize',6);
    end
end
hold off;
title('\textbf{Output multiplier}','Interpreter',Interpreter);
xlabel('Horizon of UI extension (months)','Interpreter',Interpreter);
fig=gcf;  
axis([durLabels(1),durLabels(length(durLabels)),0,2]);
yl = ylim;
set(gca,'xtick',durLabels);
xl = xlim;
box on;
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:),6)));
set(gca,'TickLabelInterpreter',Interpreter);
line([durExitZLB durExitZLB],yl,'LineWidth',0.5,'Color','k');
text(durExitZLB+(xl(2)-xl(1))/80,yl(1)+(yl(2)-yl(1))/20,...
    'i can adjust','HorizontalAlignment','left','Interpreter',Interpreter);
set(findall(fig,'-property','FontSize'),'FontSize',11);     
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure4_1']), '-depsc', '-r600','-loose');

dbarMultipliers = nan(1,length(dbarLabels));
dbarMultipliers1 = nan(1,length(dbarLabels));
dbarMultipliers12 = nan(1,length(dbarLabels));
for i=1:length(dbarLabels)
    [~,~,~,thisUImult,~] = ...
        computeUIMultiplier(RCE_dbarbar22.b,[outputName,'/',dbarInitName],[],[outputName,'/',dbarFinNames{i}]);
    [~,~,~,thisUImult1,~] = ...
        computeUIMultiplier(RCE_dbarbar22.b,[outputName,'/',dbarInitName],[],[outputName,'/',dbarFinNames1{i}]);  
    [~,~,~,thisUImult12,~] = ...
        computeUIMultiplier(RCE_dbarbar22.b,[outputName,'/',dbarInitName],[],[outputName,'/',dbarFinNames12{i}]);      
    dbarMultipliers(i) = thisUImult;
    dbarMultipliers1(i) = thisUImult1;
    dbarMultipliers12(i) = thisUImult12;
end

figure;
hold on;
for i=1:length(dbarLabels)
    if dbarLabels(i) ~= 3
        plot(dbarLabels(i),dbarMultipliers(i),'o','MarkerSize',6);
    else
        plot(dbarLabels(i),dbarMultipliers(i),'o','MarkerFaceColor',[0 0 0.5],'MarkerSize',6);
    end
end
for i=1:length(dbarLabels)
    dbar1 = plot(dbarLabels(i),dbarMultipliers1(i),'s','MarkerSize',3);
    dbar12 = plot(dbarLabels(i),dbarMultipliers12(i),'d','MarkerSize',3);
end
hold off;
title('\textbf{Output multiplier}','Interpreter',Interpreter,'FontSize',11);
xlabel('Magnitude of UI extension (months)','Interpreter',Interpreter);
fig=gcf;  
yl = ylim; 
axis([dbarLabels(1),dbarLabels(length(dbarLabels)),0,2]); 
set(gca,'xtick',dbarLabels);
xl = xlim;
box on;
curtick = get(gca, 'XTick');
set(gca, 'XTickLabel', cellstr(num2str(curtick(:),6)));
set(gca,'TickLabelInterpreter',Interpreter);
legend([dbar1 dbar12],'extension in 1st month only','extension in 12th month only','Location','southeast');
legend('boxoff');
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure4_2']), '-depsc', '-r600','-loose');
	
%% Macro shocks starting from stationary RCE (section E.1)

T = 24;

% Load init paths
init = load([outputName,'\',initName]);
iInit = (1./init.MTrans-1)*12;
rInit = (1./init.m-1)*12;
p_eInit = init.p_e;
p_uInit = init.p_u;
thetaInit = init.theta;
p_e_tildeInit = init.p_e_tilde;
sbarTransInit = init.sbarTrans;
tightnessInit = thetaInit.*sbarTransInit./(1-p_e_tildeInit);
avgWInit = init.avgW;
cTransInit = init.cTrans;
PiTransInit = init.PiTrans;
cpiInit = cumprod(1+PiTransInit)-1;
nomWInit = avgWInit.*(cpiInit+1);

betasInit = mean(init.varyingParams.betas,1);
zbarInit = init.varyingParams.zbar;
deltaBarInit = init.varyingParams.deltaBar;
abarInit = init.varyingParams.abar;
mbarInit = init.varyingParams.mbar;
lambdaInit = init.varyingParams.lambda;

for macroInd=1:length(macroShockNames)
    % Load paths following macro shock
    fin_macro = load([outputName,'\',macroShockFileNames{macroInd}]);    
    iFin_macro = (1./fin_macro.MTrans-1)*12;
    p_eFin_macro = fin_macro.p_e;
    p_uFin_macro = fin_macro.p_u;
    avgWFin_macro = fin_macro.avgW;
    cTransFin_macro = fin_macro.cTrans;
    PiTransFin_macro = fin_macro.PiTrans;
    cpiFin_macro = cumprod(1+PiTransFin_macro)-1;
    nomWFin_macro = avgWFin_macro.*(cpiFin_macro+1);
    if strcmp(macroShockNames{macroInd},'betas')
        shock_macro = mean(fin_macro.varyingParams.betas,1);
    else
        shock_macro = extractfield(fin_macro.varyingParams,macroShockNames{macroInd});
    end
    
    % Compute differentials and plot normalized shock achieving abs(urDiff_macro(1)) = 0.0005
    shockDiff_macro = shock_macro - shock_macro(end); % assumes that value at end is steady-state!
    iDiff_macro = iFin_macro-iInit;
    urDiff_macro = p_eInit - p_eFin_macro;
    ltuDiff_macro = sum(p_uFin_macro(7:end,:),1) - sum(p_uInit(7:end,:),1);
    cDiff_macro = cTransFin_macro - cTransInit;
    cpiDiff_macro = cpiFin_macro - cpiInit;
    nomWDiff_macro = nomWFin_macro - nomWInit;
    scale_macro = 0.0005/urDiff_macro(1);

    figure;
    plot(-4:1:(T-1),[zeros(1,4),shockDiff_macro(1:T)*scale_macro]);
    axis([-4 T-1 -inf inf]);
    title(['\textbf{',macroShockLabels{macroInd},'}'],'Interpreter',Interpreter);
    xlabel('Month','Interpreter',Interpreter);
    set(gca,'TickLabelInterpreter',Interpreter);
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',11); 
    set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
    print(gcf, strcat([outputName,'/figures_and_tables/FigureA',num2str(macroInd+2),'_1']), '-depsc', '-r600','-loose');       
    
    figure;
    plot(-4:1:(T-1),[zeros(1,4),iDiff_macro(1:T)*100*scale_macro]);
    if ~strcmp(macroShockNames{macroInd},'abar')
        axis([-4 T-1 -0.15 .15]);
    else
        axis([-4 T-1 -inf inf]); 
    end
    title(['\textbf{Nominal interest rate (ann.)}'],'Interpreter',Interpreter);
    ylabel('pp','Interpreter',Interpreter);
    xlabel('Month','Interpreter',Interpreter);
    set(gca,'TickLabelInterpreter',Interpreter);
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',11); 
    set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
    print(gcf, strcat([outputName,'/figures_and_tables/FigureA',num2str(macroInd+2),'_2']), '-depsc', '-r600','-loose');       
    
    figure;
    plot(-4:1:(T-1),[zeros(1,4),100*urDiff_macro(1:T)*scale_macro]); 
    axis([-4 T-1 -0.01 0.06]);
    title(['\textbf{Unemployment rate}'],'Interpreter',Interpreter);
    ylabel('pp','Interpreter',Interpreter);
    xlabel('Month','Interpreter',Interpreter);
    set(gca,'TickLabelInterpreter',Interpreter);
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',11); 
    set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);   
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
    print(gcf, strcat([outputName,'/figures_and_tables/FigureA',num2str(macroInd+2),'_3']), '-depsc', '-r600','-loose');           
    
    figure;
    plot(-4:1:(T-1),[zeros(1,4),100*(cDiff_macro(1:T)/cStatSS)*scale_macro]); 
    if ~strcmp(macroShockNames{macroInd},'abar')
        axis([-4 T-1 -0.06 0.01]);
    else
        axis([-4 T-1 -inf inf]);
    end
    title(['\textbf{Consumption per capita}'],'Interpreter',Interpreter);
    ylabel('pp','Interpreter',Interpreter);
    xlabel('Month','Interpreter',Interpreter);
    set(gca,'TickLabelInterpreter',Interpreter);
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',11); 
    set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);    
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
    print(gcf, strcat([outputName,'/figures_and_tables/FigureA',num2str(macroInd+2),'_4']), '-depsc', '-r600','-loose');               

    figure;
    plot(-4:1:(T-1),[zeros(1,4),100*(nomWDiff_macro(1:T)/(wSS*sum(phi_eSS,1)'))*scale_macro]); 
    if ~strcmp(macroShockNames{macroInd},'abar')
        axis([-4 T-1 -0.05 0.2]);
    else
        axis([-4 T-1 -inf inf]);
    end
    title(['\textbf{Average nominal wage}'],'Interpreter',Interpreter);
    ylabel('pp','Interpreter',Interpreter);
    xlabel('Month','Interpreter',Interpreter);
    set(gca,'TickLabelInterpreter',Interpreter);
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',11); 
    set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
    print(gcf, strcat([outputName,'/figures_and_tables/FigureA',num2str(macroInd+2),'_5']), '-depsc', '-r600','-loose');               
    
    figure;
    plot(-4:1:(T-1),[zeros(1,4),100*cpiDiff_macro(1:T)*scale_macro]); 
    if ~strcmp(macroShockNames{macroInd},'abar')
        axis([-4 T-1 -0.05 0.2]);
    else
        axis([-4 T-1 -inf inf]);
    end
    title(['\textbf{Final good prices}'],'Interpreter',Interpreter);
    ylabel('pp','Interpreter',Interpreter);
    xlabel('Month','Interpreter',Interpreter);
    set(gca,'TickLabelInterpreter',Interpreter);
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',11); 
    set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
    print(gcf, strcat([outputName,'/figures_and_tables/FigureA',num2str(macroInd+2),'_6']), '-depsc', '-r600','-loose');
end

%% UI extensions simulated in model (and comparison to Farber-Valletta)

ui_weeks_fv = importdata([dataName,'/',fvUI]);
Tbar = 481;
[dbarPathWeeks,dbarPathMonths,~] = genUIGreatRecession(Tbar);

T = 80; % number of calendar periods in Great Recession simulation to plot (greater than or equal to 80!)

figure;
hax=axes; 
plot(1:(T+4),[6*ones(1,4),dbarPathMonths(1,1:80),6*ones(1,T-80)]);
hold on;
plot(1:(T+4),[6*ones(1,4),round(ui_weeks_fv'/4.5),6*ones(1,T-length(ui_weeks_fv))]);
axis([1 T 0 24]);
line([5 5],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Duration of UI}'],'Interpreter',Interpreter);
ylabel('Months','Interpreter',Interpreter);
legend('Model','Data (median state)','Location','Southeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure5']), '-depsc', '-r600','-loose');

%% Great Recession simulation (sections 6 and E): model vs. data

% Load data to be replicated (relative to April 2008)
macrodata = importdata([dataName,'/',macroGR],',',1);
urdata = macrodata.data(1:end,1);
ltushdata = macrodata.data(1:end,2);
idata = macrodata.data(1:end,3);
corepcedata = macrodata.data(1:end,4);
nomWdtdata = macrodata.data(1:end,5);
cdtdata = macrodata.data(1:end,6);
tightnessdata = macrodata.data(1:end,7);

% Load Great Recession simulation using beta shocks and compare to data
load([outputName,'\',finName_ui],'shocks');
figure;
hax=axes; 
plot(1:(T+4),[zeros(1,4),shocks*1000,zeros(1,T-length(shocks))]);
hold on;
axis([1 T+4 -1 1]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Shocks to average discount factor ($\times 10^{-3}$)}'],'Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure6_1']), '-depsc', '-r600','-loose');

% Load fin paths
fin_ui = load([outputName,'\',finName_ui]);
iFin = (1./fin_ui.MPath-1)*12;
rFin = (1./fin_ui.mPath-1)*12;
p_eFin = fin_ui.p_ePath;
p_uFin = fin_ui.p_uPath;
thetaFin = fin_ui.thetaPath;
p_e_tildeFin = fin_ui.p_e_tildePath;
sbarTransFin = fin_ui.sbarPath;
tightnessFin = thetaFin.*sbarTransFin./(1-p_e_tildeFin);
avgWFin = fin_ui.avgWPath;
cTransFin = fin_ui.cPath;
PiTransFin = fin_ui.PiPath;
cpiFin = cumprod(1+PiTransFin)-1;
nomWFin = avgWFin.*(cpiFin+1);
atZLBFin = nan(1,size(fin_ui.MForward,2));
for ind=1:size(fin_ui.MForward,2)
    atZLBFin(ind) = find(fin_ui.MForward(:,ind)~=1,1,'first')-1;
end

% Compute differentials and plot
iDiff = iFin-iInit;
urDiff = p_eInit - p_eFin;
ltuDiff = sum(p_uFin(7:end,:),1) - sum(p_uInit(7:end,:),1);
tightnessDiff = thetaFin.*sbarTransFin./(1-p_e_tildeFin) - ...
    thetaInit.*sbarTransInit./(1-p_e_tildeInit);
cDiff = cTransFin - cTransInit;
cpiDiff = cpiFin - cpiInit;
nomWDiff = nomWFin - nomWInit;

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*urDiff(1:T)]); 
hold on;
plot(1:(T+4),100*urdata(1:T+4));
axis([1 T+4 -0.5 7]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Unemployment rate}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model','Data','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure6_2']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),iDiff(1:T)*100]); 
hold on;
plot(1:(T+4),ceil(idata(1:T+4)*10000/25)*25/100);
axis([1 T+4 -2.2 0]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Nominal interest rate (ann.)}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model','Data','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11);
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure7_1']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(tightnessDiff(1:T)/(thetaSS*sbarSS/(1-p_e_tildeSS)))]); 
hold on;
plot(1:(T+4),100*tightnessdata(1:T+4));
axis([1 T+4 -100 20]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Vacancies / unemployed}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure7_2']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*((ltuDiff(1:T)+sum(sum(phi_uSS(:,uDurSS>dbarSS),1))*(1-p_eSS))./(urDiff(1:T)+(1-p_eSS))-sum(sum(phi_uSS(:,uDurSS>dbarSS),1)))]); 
hold on;
plot(1:(T+4),100*ltushdata(1:T+4));
axis([1 T+4 -5 30]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Fraction long-term unemployed}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure7_3']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(cDiff(1:T)/cStatSS)]); 
hold on;
plot(1:(T+4),100*cdtdata(1:T+4));
axis([1 T+4 -15 5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Consumption per capita}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FIgure7_4']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(nomWDiff(1:T)/(wSS*sum(phi_eSS,1)'))]); 
hold on;
plot(1:(T+4),100*nomWdtdata(1:T+4));
axis([1 T+4 -5 2]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Average nominal wage}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure7_5']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*cpiDiff(1:T)]); 
hold on;
plot(1:(T+4),100*corepcedata(1:T+4));
axis([1 T+4 -5 2]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Final good prices}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure7_6']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),atZLBFin(1:T)]); 
hold on;
axis([1 T+4 0 15]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Expected horizon at ZLB}'],'Interpreter',Interpreter);
ylabel('Months','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure8']), '-depsc', '-r600','-loose');

%% Great Recession simulation (sections 6 and E): effects of UI and ZLB

% Load paths
fin_noui = load([outputName,'\',finName_noui]);
iFin_noui = (1./fin_noui.MPath-1)*12;
rFin_noui = (1./fin_noui.mPath-1)*12;
p_eFin_noui = fin_noui.p_ePath;
p_uFin_noui = fin_noui.p_uPath;
thetaFin_noui = fin_noui.thetaPath;
p_e_tildeFin_noui = fin_noui.p_e_tildePath;
sbarTransFin_noui = fin_noui.sbarPath;
tightnessFin_noui = thetaFin_noui.*sbarTransFin_noui./(1-p_e_tildeFin_noui);
avgWFin_noui = fin_noui.avgWPath;
cTransFin_noui = fin_noui.cPath;
PiTransFin_noui = fin_noui.PiPath;
cpiFin_noui = cumprod(1+PiTransFin_noui)-1;
nomWFin_noui = avgWFin_noui.*(cpiFin_noui+1);

fin_nozlb = load([outputName,'\',finName_nozlb]);
iFin_nozlb = (1./fin_nozlb.MPath-1)*12;
rFin_nozlb = (1./fin_nozlb.mPath-1)*12;
p_eFin_nozlb = fin_nozlb.p_ePath;
p_uFin_nozlb = fin_nozlb.p_uPath;
thetaFin_nozlb = fin_nozlb.thetaPath;
p_e_tildeFin_nozlb = fin_nozlb.p_e_tildePath;
sbarTransFin_nozlb = fin_nozlb.sbarPath;
tightnessFin_nozlb = thetaFin_nozlb.*sbarTransFin_nozlb./(1-p_e_tildeFin_nozlb);
avgWFin_nozlb = fin_nozlb.avgWPath;
cTransFin_nozlb = fin_nozlb.cPath;
PiTransFin_nozlb = fin_nozlb.PiPath;
cpiFin_nozlb = cumprod(1+PiTransFin_nozlb)-1;
nomWFin_nozlb = avgWFin_nozlb.*(cpiFin_nozlb+1);

fin_nozlbui = load([outputName,'\',finName_nozlbui]);
iFin_nozlbui = (1./fin_nozlbui.MPath-1)*12;
rFin_nozlbui = (1./fin_nozlbui.mPath-1)*12;
p_eFin_nozlbui = fin_nozlbui.p_ePath;
p_uFin_nozlbui = fin_nozlbui.p_uPath;
thetaFin_nozlbui = fin_nozlbui.thetaPath;
p_e_tildeFin_nozlbui = fin_nozlbui.p_e_tildePath;
sbarTransFin_nozlbui = fin_nozlbui.sbarPath;
tightnessFin_nozlbui = thetaFin_nozlbui.*sbarTransFin_nozlbui./(1-p_e_tildeFin_nozlbui);
avgWFin_nozlbui = fin_nozlbui.avgWPath;
cTransFin_nozlbui = fin_nozlbui.cPath;
PiTransFin_nozlbui = fin_nozlbui.PiPath;
cpiFin_nozlbui = cumprod(1+PiTransFin_nozlbui)-1;
nomWFin_nozlbui = avgWFin_nozlbui.*(cpiFin_nozlbui+1);

% Compute differentials and plot
iDiff_noui = iFin_noui-iInit;
urDiff_noui = p_eInit - p_eFin_noui;
ltuDiff_noui = sum(p_uFin_noui(7:end,:),1) - sum(p_uInit(7:end,:),1);
tightnessDiff_noui = thetaFin_noui.*sbarTransFin_noui./(1-p_e_tildeFin_noui) - ...
    thetaInit.*sbarTransInit./(1-p_e_tildeInit);
cDiff_noui = cTransFin_noui - cTransInit;
cpiDiff_noui = cpiFin_noui - cpiInit;
nomWDiff_noui = nomWFin_noui - nomWInit;

iDiff_nozlb = iFin_nozlb-iInit;
urDiff_nozlb = p_eInit - p_eFin_nozlb;
ltuDiff_nozlb = sum(p_uFin_nozlb(7:end,:),1) - sum(p_uInit(7:end,:),1);
tightnessDiff_nozlb = thetaFin_nozlb.*sbarTransFin_nozlb./(1-p_e_tildeFin_nozlb) - ...
    thetaInit.*sbarTransInit./(1-p_e_tildeInit);
cDiff_nozlb = cTransFin_nozlb - cTransInit;
cpiDiff_nozlb = cpiFin_nozlb - cpiInit;
nomWDiff_nozlb = nomWFin_nozlb - nomWInit;

iDiff_nozlbui = iFin_nozlbui-iInit;
urDiff_nozlbui = p_eInit - p_eFin_nozlbui;
ltuDiff_nozlbui = sum(p_uFin_nozlbui(7:end,:),1) - sum(p_uInit(7:end,:),1);
tightnessDiff_nozlbui = thetaFin_nozlbui.*sbarTransFin_nozlbui./(1-p_e_tildeFin_nozlbui) - ...
    thetaInit.*sbarTransInit./(1-p_e_tildeInit);
cDiff_nozlbui = cTransFin_nozlbui - cTransInit;
cpiDiff_nozlbui = cpiFin_nozlbui - cpiInit;
nomWDiff_nozlbui = nomWFin_nozlbui - nomWInit;

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*urDiff(1:T)]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*urDiff_noui(1:T)]); 
axis([1 T+4 -0.5 5.5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Unemployment rate}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model','No UI extensions','Location','Southeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11);
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure9_1']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(tightnessDiff(1:T)/(thetaSS*sbarSS/(1-p_e_tildeSS)))]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*(tightnessDiff_noui(1:T)/(thetaSS*sbarSS/(1-p_e_tildeSS)))]); 
axis([1 T+4 -60 10]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Vacancies / unemployed}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure9_2']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*((ltuDiff(1:T)+sum(sum(phi_uSS(:,uDurSS>dbarSS),1))*(1-p_eSS))./(urDiff(1:T)+(1-p_eSS))-sum(sum(phi_uSS(:,uDurSS>dbarSS),1)))]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*((ltuDiff_noui(1:T)+sum(sum(phi_uSS(:,uDurSS>dbarSS),1))*(1-p_eSS))./(urDiff_noui(1:T)+(1-p_eSS))-sum(sum(phi_uSS(:,uDurSS>dbarSS),1)))]); 
axis([1 T+4 -3 20]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Fraction long-term unemployed}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure9_3']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(cDiff(1:T)/cStatSS)]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*(cDiff_noui(1:T)/cStatSS)]); 
axis([1 T+4 -5 1]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Consumption per capita}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure9_4']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(nomWDiff(1:T)/(wSS*sum(phi_eSS,1)'))]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*(nomWDiff_noui(1:T)/(wSS*sum(phi_eSS,1)'))]); 
axis([1 T+4 -3.5 0.5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Average nominal wage}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure9_5']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*cpiDiff(1:T)]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*cpiDiff_noui(1:T)]); 
axis([1 T+4 -3.5 0.5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Final good prices}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure9_6']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*urDiff(1:T)]); 
hold on; 
plot(1:(T+4),[zeros(1,4),100*urDiff_nozlb(1:T)]); 
plot(1:(T+4),[zeros(1,4),100*urDiff_nozlbui(1:T)]); 
axis([1 T+4 -0.5 7.5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Unemployment rate}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model','No ZLB','No ZLB, no UI extensions','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure10_1']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),iDiff(1:T)*100]); 
hold on;
plot(1:(T+4),[zeros(1,4),iDiff_nozlb(1:T)*100]); 
plot(1:(T+4),[zeros(1,4),iDiff_nozlbui(1:T)*100]); 
axis([1 T+4 -4.5 0]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Nominal interest rate (ann.)}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model','No ZLB','No ZLB, no UI extensions','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11);
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure10_2']), '-depsc', '-r600','-loose');

%% Great Recession simulation (sections 6 and E): role of real wage rigidity

% Load fin paths under high alternative iota
fin_ui_iotaHi = load([outputName,'\',finName_ui_iotaHi]);
avgWFin_iotaHi = fin_ui_iotaHi.avgWPath;
PiTransFin_iotaHi = fin_ui_iotaHi.PiPath;
cpiFin_iotaHi = cumprod(1+PiTransFin_iotaHi)-1;
nomWFin_iotaHi = avgWFin_iotaHi.*(cpiFin_iotaHi+1);

% Load fin paths under low alternative iota
fin_ui_iotaLo = load([outputName,'\',finName_ui_iotaLo]);
avgWFin_iotaLo = fin_ui_iotaLo.avgWPath;
PiTransFin_iotaLo = fin_ui_iotaLo.PiPath;
cpiFin_iotaLo = cumprod(1+PiTransFin_iotaLo)-1;
nomWFin_iotaLo = avgWFin_iotaLo.*(cpiFin_iotaLo+1);

% Compute differentials and plot
cpiDiff_iotaHi = cpiFin_iotaHi - cpiInit;
nomWDiff_iotaHi = nomWFin_iotaHi - nomWInit;

cpiDiff_iotaLo = cpiFin_iotaLo - cpiInit;
nomWDiff_iotaLo = nomWFin_iotaLo - nomWInit;

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(nomWDiff(1:T)/(wSS*sum(phi_eSS,1)'))]); 
hold on;
plot(1:(T+4),100*nomWdtdata(1:T+4));
plot(1:(T+4),[zeros(1,4),100*(nomWDiff_iotaHi(1:T)/(wSS*sum(phi_eSS,1)'))]); 
plot(1:(T+4),[zeros(1,4),100*(nomWDiff_iotaLo(1:T)/(wSS*sum(phi_eSS,1)'))]); 
axis([1 T+4 -5 2]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Average nominal wage}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model','Data','$\iota = 1$','$\iota = 0.88$','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure11_1']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*cpiDiff(1:T)]); 
hold on;
plot(1:(T+4),100*corepcedata(1:T+4));
plot(1:(T+4),[zeros(1,4),100*cpiDiff_iotaHi(1:T)]);
plot(1:(T+4),[zeros(1,4),100*cpiDiff_iotaLo(1:T)]);
axis([1 T+4 -5 2]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Final good prices}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/Figure11_2']), '-depsc', '-r600','-loose');

%% Great Recession simulation (sections 6 and E): with separation rate shocks

fin_ui_seprate = load([outputName,'\',finName_ui_seprate]);
figure;
hax=axes; 
plot(1:(T+4),[zeros(1,4),1000*fin_ui_seprate.shocks,zeros(1,T-length(fin_ui_seprate.shocks))]);
hold on;
axis([1 T+4 -1 1]);
plot(1:(T+4),[zeros(1,4),1000*shocks,zeros(1,T-length(shocks))],'-','LineWidth',0.5);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Shocks to average discount factor ($\times 10^{-3}$)}'],'Interpreter',Interpreter);
legend('Model with separation rate shocks','Model','Location','Northeast');
%xlabel('Month');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA9_1']), '-depsc', '-r600','-loose');

deltaBarShocks = importdata([dataName,'/',seprateshocks]);
figure;
hax=axes; 
plot(1:(T+4),[zeros(1,4),100*deltaBarShocks',zeros(1,T-length(deltaBarShocks))]);
hold on;
axis([1 T+4 -0.8 1.2]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Shocks to aggregate separation rate ($\times 10^{-2}$)}'],'Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA8_2']), '-depsc', '-r600','-loose');

deltaBar9019 = importdata([dataName,'/',sepratedata]);
figure;
plot(deltaBar9019(:,1)+deltaBar9019(:,2)/12,deltaBar9019(:,3));
hold on;
plot(deltaBar9019(:,1)+deltaBar9019(:,2)/12,deltaBar9019(:,3)-deltaBar9019(:,4),'k-','LineWidth',2);
hold off;      
title('\textbf{Aggregate separation rate}','Interpreter',Interpreter);
set(gca,'TickLabelInterpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11);
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
axis([-inf,inf,0,0.06]);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA8_1']), '-depsc', '-r600','-loose');     

% Load fin paths
iFin_seprate = (1./fin_ui_seprate.MPath-1)*12;
rFin_seprate = (1./fin_ui_seprate.mPath-1)*12;
p_eFin_seprate = fin_ui_seprate.p_ePath;
p_uFin_seprate = fin_ui_seprate.p_uPath;
thetaFin_seprate = fin_ui_seprate.thetaPath;
p_e_tildeFin_seprate = fin_ui_seprate.p_e_tildePath;
sbarTransFin_seprate = fin_ui_seprate.sbarPath;
tightnessFin_seprate = thetaFin_seprate.*sbarTransFin_seprate./(1-p_e_tildeFin_seprate);
avgWFin_seprate = fin_ui_seprate.avgWPath;
cTransFin_seprate = fin_ui_seprate.cPath;
PiTransFin_seprate = fin_ui_seprate.PiPath;
cpiFin_seprate = cumprod(1+PiTransFin_seprate)-1;
nomWFin_seprate = avgWFin_seprate.*(cpiFin_seprate+1);

% Compute differentials and plot
iDiff_seprate = iFin_seprate-iInit;
urDiff_seprate = p_eInit - p_eFin_seprate;
ltuDiff_seprate = sum(p_uFin_seprate(7:end,:),1) - sum(p_uInit(7:end,:),1);
tightnessDiff_seprate = thetaFin_seprate.*sbarTransFin_seprate./(1-p_e_tildeFin_seprate) - ...
    thetaInit.*sbarTransInit./(1-p_e_tildeInit);
cDiff_seprate = cTransFin_seprate - cTransInit;
cpiDiff_seprate = cpiFin_seprate - cpiInit;
nomWDiff_seprate = nomWFin_seprate - nomWInit;

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*urDiff_seprate(1:T)]); 
hold on;
plot(1:(T+4),100*urdata(1:T+4));
axis([1 T+4 -0.5 7]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Unemployment rate}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model with separation rate shocks','Data','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA9_2']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),iDiff_seprate(1:T)*100]); 
hold on;
plot(1:(T+4),ceil(idata(1:T+4)*10000/25)*25/100);
axis([1 T+4 -2.2 0]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Nominal interest rate (ann.)}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model with separation rate shocks','Data','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11);
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA10_1']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(tightnessDiff_seprate(1:T)/(thetaSS*sbarSS/(1-p_e_tildeSS)))]); 
hold on;
plot(1:(T+4),100*tightnessdata(1:T+4));
axis([1 T+4 -100 20]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Vacancies / unemployed}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA10_2']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*((ltuDiff_seprate(1:T)+sum(sum(phi_uSS(:,uDurSS>dbarSS),1))*(1-p_eSS))./(urDiff_seprate(1:T)+(1-p_eSS))-sum(sum(phi_uSS(:,uDurSS>dbarSS),1)))]); 
hold on;
plot(1:(T+4),100*ltushdata(1:T+4));
axis([1 T+4 -5 30]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Fraction long-term unemployed}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA10_3']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(cDiff_seprate(1:T)/cStatSS)]); 
hold on;
plot(1:(T+4),100*cdtdata(1:T+4));
axis([1 T+4 -15 5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Consumption per capita}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA10_4']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(nomWDiff_seprate(1:T)/(wSS*sum(phi_eSS,1)'))]); 
hold on;
plot(1:(T+4),100*nomWdtdata(1:T+4));
axis([1 T+4 -5 2]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Average nominal wage}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA10_5']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*cpiDiff_seprate(1:T)]); 
hold on;
plot(1:(T+4),100*corepcedata(1:T+4));
axis([1 T+4 -5 2]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Final good prices}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA10_6']), '-depsc', '-r600','-loose');

% Compare wages and prices vs. baseline calibration

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(nomWDiff_seprate(1:T)/(wSS*sum(phi_eSS,1)'))]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*(nomWDiff(1:T)/(wSS*sum(phi_eSS,1)'))]); 
axis([1 T+4 -3.5 0.5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Average nominal wage}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model with separation rate shocks','Model','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA11_1']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*cpiDiff_seprate(1:T)]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*cpiDiff(1:T)]); 
axis([1 T+4 -3.5 0.5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Final good prices}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA11_2']), '-depsc', '-r600','-loose');

% Effect of UI extensions on unemployment 

fin_noui_seprate = load([outputName,'\',finName_noui_seprate]);
iFin_noui_seprate = (1./fin_noui_seprate.MPath-1)*12;
rFin_noui_seprate = (1./fin_noui_seprate.mPath-1)*12;
p_eFin_noui_seprate = fin_noui_seprate.p_ePath;
p_uFin_noui_seprate = fin_noui_seprate.p_uPath;
thetaFin_noui_seprate = fin_noui_seprate.thetaPath;
p_e_tildeFin_noui_seprate = fin_noui_seprate.p_e_tildePath;
sbarTransFin_noui_seprate = fin_noui_seprate.sbarPath;
tightnessFin_noui_seprate = thetaFin_noui_seprate.*sbarTransFin_noui_seprate./(1-p_e_tildeFin_noui_seprate);
avgWFin_noui_seprate = fin_noui_seprate.avgWPath;
cTransFin_noui_seprate = fin_noui_seprate.cPath;
PiTransFin_noui_seprate = fin_noui_seprate.PiPath;
cpiFin_noui_seprate = cumprod(1+PiTransFin_noui_seprate)-1;
nomWFin_noui_seprate = avgWFin_noui_seprate.*(cpiFin_noui_seprate+1);

iDiff_noui_seprate = iFin_noui_seprate-iInit;
urDiff_noui_seprate = p_eInit - p_eFin_noui_seprate;
ltuDiff_noui_seprate = sum(p_uFin_noui_seprate(7:end,:),1) - sum(p_uInit(7:end,:),1);
tightnessDiff_noui_seprate = thetaFin_noui_seprate.*sbarTransFin_noui_seprate./(1-p_e_tildeFin_noui_seprate) - ...
    thetaInit.*sbarTransInit./(1-p_e_tildeInit);
cDiff_noui_seprate = cTransFin_noui_seprate - cTransInit;
cpiDiff_noui_seprate = cpiFin_noui_seprate - cpiInit;
nomWDiff_noui_seprate = nomWFin_noui_seprate - nomWInit;

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*urDiff_seprate(1:T)]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*urDiff_noui_seprate(1:T)]); 
axis([1 T+4 -0.5 5.5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Unemployment rate}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model with separation rate shocks','No UI extensions','Location','Southeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11);
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA12_1']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(tightnessDiff_seprate(1:T)/(thetaSS*sbarSS/(1-p_e_tildeSS)))]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*(tightnessDiff_noui_seprate(1:T)/(thetaSS*sbarSS/(1-p_e_tildeSS)))]); 
axis([1 T+4 -60 10]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Vacancies / unemployed}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA12_2']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*((ltuDiff_seprate(1:T)+sum(sum(phi_uSS(:,uDurSS>dbarSS),1))*(1-p_eSS))./(urDiff_seprate(1:T)+(1-p_eSS))-sum(sum(phi_uSS(:,uDurSS>dbarSS),1)))]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*((ltuDiff_noui_seprate(1:T)+sum(sum(phi_uSS(:,uDurSS>dbarSS),1))*(1-p_eSS))./(urDiff_noui_seprate(1:T)+(1-p_eSS))-sum(sum(phi_uSS(:,uDurSS>dbarSS),1)))]); 
axis([1 T+4 -3 20]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Fraction long-term unemployed}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA12_3']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(cDiff_seprate(1:T)/cStatSS)]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*(cDiff_noui_seprate(1:T)/cStatSS)]); 
axis([1 T+4 -5 1]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Consumption per capita}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA12_4']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*(nomWDiff_seprate(1:T)/(wSS*sum(phi_eSS,1)'))]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*(nomWDiff_noui_seprate(1:T)/(wSS*sum(phi_eSS,1)'))]); 
axis([1 T+4 -3.5 0.5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Average nominal wage}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA12_5']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(1:(T+4),[zeros(1,4),100*cpiDiff_seprate(1:T)]); 
hold on;
plot(1:(T+4),[zeros(1,4),100*cpiDiff_noui_seprate(1:T)]); 
axis([1 T+4 -3.5 0.5]);
line([4 4],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = 12*(1:7);
ax.XTickLabel = {'12/08','12/09','12/10','12/11','12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Final good prices}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA12_6']), '-depsc', '-r600','-loose');

%% Great Recession simulation (sections 6 and E): unexpected expiration in Dec 2013

startT = 53; % initial period to start plots
endT = 80; % end period of plot

% Load paths with expiration
fin_expiration = load([outputName,'\',finName_expiration]);
shocksFin_expiration = fin_expiration.shocks;
p_eFin_expiration = fin_expiration.p_ePath;
thetaFin_expiration = fin_expiration.thetaPath;
p_e_tildeFin_expiration = fin_expiration.p_e_tildePath;
sbarTransFin_expiration = fin_expiration.sbarPath;
tightnessFin_expiration = thetaFin_expiration.*sbarTransFin_expiration./(1-p_e_tildeFin_expiration);

% Compute differentials and plot
urDiff_expiration = p_eInit - p_eFin_expiration;
tightnessDiff_expiration = thetaFin_expiration.*sbarTransFin_expiration./(1-p_e_tildeFin_expiration) - ...
    thetaInit.*sbarTransInit./(1-p_e_tildeInit);

figure;
hax=axes; 
plot(startT:endT,1000*shocksFin_expiration(startT:endT));
hold on;
plot(startT:endT,1000*shocks(startT:endT),'-','LineWidth',0.5);
axis([startT endT -0.7 0.7]);
line([startT+3 startT+3],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = startT+[3,15,27];
ax.XTickLabel = {'12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Shocks to average discount factor ($\times 10^{-3}$)}'],'Interpreter',Interpreter);
legend('Model with expiration in 12/13','Model','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA13']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(startT:endT,100*urDiff_expiration(startT:endT)); 
hold on;
plot(startT:endT,100*urdata(4+startT:4+endT));
plot(startT:endT,100*urDiff(startT:endT),'-','LineWidth',0.5);
axis([startT endT -0.5 4.5]);
line([startT+3 startT+3],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = startT+[3,15,27];
ax.XTickLabel = {'12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Unemployment rate}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model with expiration in 12/13','Data','Model','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA14_1']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(startT:endT,100*(tightnessDiff_expiration(startT:endT)/(thetaSS*sbarSS/(1-p_e_tildeSS)))); 
hold on;
plot(startT:endT,100*tightnessdata(4+startT:4+endT));
plot(startT:endT,100*(tightnessDiff(startT:endT)/(thetaSS*sbarSS/(1-p_e_tildeSS))),'-','LineWidth',0.5); 
axis([startT endT -60 20]);
line([startT+3 startT+3],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = startT+[3,15,27];
ax.XTickLabel = {'12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Vacancies / unemployed}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA14_2']), '-depsc', '-r600','-loose');

% Load paths with no expiration
fin_noexpiration = load([outputName,'\',finName_noexpiration]);
p_eFin_noexpiration = fin_noexpiration.p_ePath;
thetaFin_noexpiration = fin_noexpiration.thetaPath;
p_e_tildeFin_noexpiration = fin_noexpiration.p_e_tildePath;
sbarTransFin_noexpiration = fin_noexpiration.sbarPath;
tightnessFin_noexpiration = thetaFin_noexpiration.*sbarTransFin_noexpiration./(1-p_e_tildeFin_noexpiration);

% Compute differentials and plot
urDiff_noexpiration = p_eInit - p_eFin_noexpiration;
tightnessDiff_noexpiration = thetaFin_noexpiration.*sbarTransFin_noexpiration./(1-p_e_tildeFin_noexpiration) - ...
    thetaInit.*sbarTransInit./(1-p_e_tildeInit);

figure;
hax=axes;
plot(startT:endT,100*urDiff_expiration(startT:endT)); 
hold on;
plot(startT:endT,100*urDiff_noexpiration(startT:endT)); 
axis([startT endT -0.5 4.5]);
line([startT+3 startT+3],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = startT+[3,15,27];
ax.XTickLabel = {'12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Unemployment rate}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
legend('Model with expiration in 12/13','No expiration in 12/13','Location','Northeast');
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11);
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA15_1']), '-depsc', '-r600','-loose');

figure;
hax=axes;
plot(startT:endT,100*(tightnessDiff_expiration(startT:endT)/(thetaSS*sbarSS/(1-p_e_tildeSS)))); 
hold on;
plot(startT:endT,100*(tightnessDiff_noexpiration(startT:endT)/(thetaSS*sbarSS/(1-p_e_tildeSS)))); 
axis([startT endT -60 20]);
line([startT+3 startT+3],get(hax,'YLim'),'LineWidth',0.5,'Color','k');
hold off;
ax = gca;
ax.XTick = startT+[3,15,27];
ax.XTickLabel = {'12/12','12/13','12/14'};
ax.TickLabelInterpreter = Interpreter;
title(['\textbf{Vacancies / unemployed}'],'Interpreter',Interpreter);
ylabel('pp','Interpreter',Interpreter);
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',11); 
set(findall(fig,'-property','Interpreter'),'Interpreter',Interpreter); 
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0 0 4 3.25]);
print(gcf, strcat([outputName,'/figures_and_tables/FigureA15_2']), '-depsc', '-r600','-loose');

disp(['Code completed in ',num2str(toc),' s.']);