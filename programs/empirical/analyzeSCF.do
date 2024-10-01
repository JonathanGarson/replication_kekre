/* analyzeSCF.do 
   Uses 2004 SCF to estimate moments of wealth distribution */

#delimit ;

clear all;
set maxvar 10000;

timer clear 1;
timer on 1;

local sourceFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\data";
local dtaOutputFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\programs\empirical\output";
local figtabOutputFolder = "`dtaOutputFolder'\figures_and_tables";
local modelOutputFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\programs\model\output";

global makedata = 1;

if $makedata {;
	use "`sourceFolder'\2004_SCF\p04i6.dta", clear;
	keep YY1 Y1 X42000 X42001 X432 X7575 X6670 X6671 X6672 X6673 X6674 X6675 X6676 X6677 X6780 X6781;
	sort Y1;
	save "`dtaOutputFolder'\2004_SCF_1.dta", replace;

	use "`sourceFolder'\2004_SCF\rscfp2004.dta", clear;
	keep YY1 Y1 wgt liq cds nmmf stocks bond govtbnd retqliq ccbal lf income fin nfin debt networth othma othfin bus asset houses oresre mrthel resdbt equity;
	sort Y1;
	save "`dtaOutputFolder'\2004_SCF_2.dta", replace;

	merge Y1 using "`dtaOutputFolder'\2004_SCF_1.dta";
	save "`dtaOutputFolder'\2004_SCF.dta", replace;
	clear;
};

use "`dtaOutputFolder'\2004_SCF.dta", replace;

* Define employment status;

gen currentemp = 1 if X6670 == 1 | X6671 == 1 | X6672 == 1 | X6673 == 1 | X6674 == 1 | X6675 == 1 | X6676 == 1 | X6677 == 1;
replace currentemp = 0 if X6670 == 3 | X6671 == 3 | X6672 == 3 | X6673 == 3 | X6674 == 3 | X6675 == 3 | X6676 == 3 | X6677 == 3;
keep if currentemp == 1 | currentemp == 0;

* Characterize median and mean moments of wealth distribution;

* Scale by mean monthly earnings;
summarize income [aw=wgt] if ~missing(currentemp);
scalar monthlyInc = r(mean)/12;
scalar monthlyInc2004 = monthlyInc*188.9/232.957; 
 * in 2004 dollars (188.9 and 232.957 are 2004 and 2013 CPIs from BLS, respectively);
 * = $6,761;

gen liqS = liq/monthlyInc;
gen bondS = bond/monthlyInc;
gen otherfinS = (fin - liq - bond)/monthlyInc;
gen nfinS = nfin/monthlyInc;
gen ccbalS = -ccbal/monthlyInc;
gen otherdebtS = -(debt - ccbal)/monthlyInc;
gen laS = (liq + bond - ccbal)/monthlyInc;
gen networthS = networth/monthlyInc;

* For each of the main variables of interest, now compute summary stats of interest;

matrix TableA2 = J(9,5,.);
matrix colnames TableA2 = "median, emp" "median, unemp" "mean, emp" "mean, unemp" "mean";
matrix rownames TableA2 = "trans" "bonds" "oth fin" "non-fin" "cc" "oth debt" "liquid" "total" "hh";

* Transaction accounts;
tabstat liqS [aw=wgt], by(currentemp) statistics(median mean n) save;
matrix TableA2[1,1] = el(r(Stat2),1,1);
matrix TableA2[1,3] = el(r(Stat2),2,1);
matrix TableA2[1,2] = el(r(Stat1),1,1);
matrix TableA2[1,4] = el(r(Stat1),2,1);
matrix TableA2[1,5] = el(r(StatTotal),2,1);
matrix TableA2[9,1] = el(r(Stat2),3,1)/5;
matrix TableA2[9,3] = el(r(Stat2),3,1)/5;
matrix TableA2[9,2] = el(r(Stat1),3,1)/5;
matrix TableA2[9,4] = el(r(Stat1),3,1)/5;
matrix TableA2[9,5] = el(r(StatTotal),3,1)/5;

* Directly held bonds;
tabstat bondS [aw=wgt], by(currentemp) statistics(median mean) save;
matrix TableA2[2,1] = el(r(Stat2),1,1);
matrix TableA2[2,3] = el(r(Stat2),2,1);
matrix TableA2[2,2] = el(r(Stat1),1,1);
matrix TableA2[2,4] = el(r(Stat1),2,1);
matrix TableA2[2,5] = el(r(StatTotal),2,1);

* Other financial assets;
tabstat otherfinS [aw=wgt], by(currentemp) statistics(median mean) save;
matrix TableA2[3,1] = el(r(Stat2),1,1);
matrix TableA2[3,3] = el(r(Stat2),2,1);
matrix TableA2[3,2] = el(r(Stat1),1,1);
matrix TableA2[3,4] = el(r(Stat1),2,1);
matrix TableA2[3,5] = el(r(StatTotal),2,1);

* Non-financial assets;
tabstat nfinS [aw=wgt], by(currentemp) statistics(median mean) save;
matrix TableA2[4,1] = el(r(Stat2),1,1);
matrix TableA2[4,3] = el(r(Stat2),2,1);
matrix TableA2[4,2] = el(r(Stat1),1,1);
matrix TableA2[4,4] = el(r(Stat1),2,1);
matrix TableA2[4,5] = el(r(StatTotal),2,1);

* Credit card debt;
tabstat ccbalS [aw=wgt], by(currentemp) statistics(median mean) save;
matrix TableA2[5,1] = el(r(Stat2),1,1);
matrix TableA2[5,3] = el(r(Stat2),2,1);
matrix TableA2[5,2] = el(r(Stat1),1,1);
matrix TableA2[5,4] = el(r(Stat1),2,1);
matrix TableA2[5,5] = el(r(StatTotal),2,1);

* Other debt;
tabstat otherdebtS [aw=wgt], by(currentemp) statistics(median mean) save;
matrix TableA2[6,1] = el(r(Stat2),1,1);
matrix TableA2[6,3] = el(r(Stat2),2,1);
matrix TableA2[6,2] = el(r(Stat1),1,1);
matrix TableA2[6,4] = el(r(Stat1),2,1);
matrix TableA2[6,5] = el(r(StatTotal),2,1);

* Liquid wealth;
tabstat laS [aw=wgt], by(currentemp) statistics(median mean) save;
matrix TableA2[7,1] = el(r(Stat2),1,1);
matrix TableA2[7,3] = el(r(Stat2),2,1);
matrix TableA2[7,2] = el(r(Stat1),1,1);
matrix TableA2[7,4] = el(r(Stat1),2,1);
matrix TableA2[7,5] = el(r(StatTotal),2,1);

* Total net worth;
tabstat networthS [aw=wgt], by(currentemp) statistics(median mean) save;
matrix TableA2[8,1] = el(r(Stat2),1,1);
matrix TableA2[8,3] = el(r(Stat2),2,1);
matrix TableA2[8,2] = el(r(Stat1),1,1);
matrix TableA2[8,4] = el(r(Stat1),2,1);
matrix TableA2[8,5] = el(r(StatTotal),2,1);

matrix list TableA2;
putexcel set "`figtabOutputFolder'\TableA2.xlsx", replace;
putexcel A1 = matrix(TableA2), names;
 
* Characterize equity shares by net worth quantile;

gen corpexposure = (equity + bus) / networth;

quietly xtile networthquant = networth [aw=wgt], nq(20);
tabstat corpexposure [aw=wgt], by(networthquant) statistics(median) save;

matrix corpexposure_by_wealth = J(20,2,.); 
forvalues q=1/20 {;
	matrix corpexposure_by_wealth[`q',1] = `q';
	matrix corpexposure_by_wealth[`q',2] = el(r(Stat`q'),1,1);
};

matrix list corpexposure_by_wealth;
putexcel set "`dtaOutputFolder'\corpexposure_by_wealth.xlsx", replace;
putexcel A1 = matrix(corpexposure_by_wealth);
putexcel set "`modelOutputFolder'\corpexposure_by_wealth.xlsx", replace;
putexcel A1 = matrix(corpexposure_by_wealth); 

timer off 1;
timer list 1;
