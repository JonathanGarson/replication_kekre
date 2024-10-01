/* analyzeRV2017.do 
   Uses Rothstein-Valletta (2017) extract from SIPP to summarize income through unemployment */

#delimit ;

clear all;
set maxvar 10000;

timer clear 1;
timer on 1;

local sourceFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\data";
local dtaOutputFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\programs\empirical\output";
local figtabOutputFolder = "`dtaOutputFolder'\figures_and_tables";
local modelOutputFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\programs\model\output";

global restrictHHHeads = 1;
 * RV(2017) use all individuals, whereas I restrict attention to HH heads;

cd "`sourceFolder'\RV_2017";
use if uiexhaustsamp==1 & yl_preue<. using eventstudysample_3-20-17.dta, clear;

quietly keep if ~($restrictHHHeads* (errp ~= 1 & errp ~= 2));

*Define pre-separation, UI receipt, and post-exhaustion periods;
* (indicator period (0,1,2) defines these);
gen post1=(spelltime_start>0) if spelltime_start>=-4 & spelltime_start<=6 & (spellongoing==1 | spelltime_start<=0);
replace post1=. if (spelltime_start==-1 | spelltime_start==0);

gen post2=(timetilluiend>0) if timetilluiend>=-3 & timetilluiend<=1;
replace post2=. if mnum > spellenddate;
replace post2=. if timetilluiend==0;

gen period = 0 if post1 == 0;
replace period = 1 if post2 == 0;
replace period = 2 if post2 == 1;
drop if period==.;

*Summarize means and standard errors of income sources by period;
collapse (mean) rely_total rely_earn rely_ui rely_earno rely_snsoc rely_socsec ,  by(ssuid id spell period);
sort id spell period;
by id spell: keep if _N==3;
by id spell: assert period[1]==0 & period[2]==1 & period[3]==2;

matrix TableA4 = J(7,6,.);
matrix colnames TableA4 = "total" "own" "ui" "other hh" "snap + welf" "ss";
matrix rownames TableA4 = "mean, prior job loss" "se" "mean, during UI" "se" "mean, after UI" "se" "obs";

mean rely_total rely_earn rely_ui rely_earno rely_snsoc rely_socsec if period==0, cluster(ssuid);
forvalues col=1/6 {;
	matrix TableA4[1,`col'] = el(r(table),1,`col');
	matrix TableA4[2,`col'] = el(r(table),2,`col');
};
mean rely_total rely_earn rely_ui rely_earno rely_snsoc rely_socsec if period==1, cluster(ssuid);
forvalues col=1/6 {;
	matrix TableA4[3,`col'] = el(r(table),1,`col');
	matrix TableA4[4,`col'] = el(r(table),2,`col');
};
mean rely_total rely_earn rely_ui rely_earno rely_snsoc rely_socsec if period==2, cluster(ssuid);
forvalues col=1/6 {;
	matrix TableA4[5,`col'] = el(r(table),1,`col');
	matrix TableA4[6,`col'] = el(r(table),2,`col');
};
tabstat rely_total rely_earn rely_ui rely_earno rely_snsoc rely_socsec if period==0, statistics(N) save;
forvalues col=1/6 {;
	matrix TableA4[7,`col'] = el(r(StatTotal),1,`col');
};

matrix list TableA4;
putexcel set "`figtabOutputFolder'\TableA4.xlsx", replace;
putexcel A1 = matrix(TableA4), names;

timer off 1;
timer list 1;
