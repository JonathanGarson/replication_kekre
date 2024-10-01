/* analyzeSHIW.do 
   Uses 2010 SHIW to estimate average MPCs of unemployed vs. employed */
 
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

if $makedata==1{;

	cd "`sourceFolder'\2010_SHIW\";

	use carcom10.dta;

	merge m:1 nquest using q10e.dta;

	keep if cfred == 1; // keep only head of households

	keep nquest nord pesofit apqual riscons1 riscons2;

	cd "`dtaOutputFolder'\";

	save 2010_SHIW.dta, replace;
	clear;

};

cd "`dtaOutputFolder'\";

use 2010_SHIW.dta;

* Define employment status;

gen empstatus = "employed" if apqual <= 10 | apqual == 20;
replace empstatus = "unemployed" if apqual == 12;
replace empstatus = "other" if apqual == 11 | apqual == 21 | (apqual > 12 & apqual < 20);

* Characterize mean MPC and number of observations by employment status;

matrix TableA1 = J(4,2,.);
matrix rownames TableA1 = "ann MPC emp" "se" "ann MPC unemp" "se";
matrix colnames TableA1 = "mean" "obs";

mean riscons2 if empstatus == "employed" [aweight=pesofit];
matrix TableA1[1,1] = el(r(table),1,1)/100;
matrix TableA1[2,1] = el(r(table),2,1)/100;
matrix TableA1[1,2] = el(r(table),7,1)+1;

mean riscons2 if empstatus == "unemployed" [aweight=pesofit];
matrix TableA1[3,1] = el(r(table),1,1)/100;
matrix TableA1[4,1] = el(r(table),2,1)/100;
matrix TableA1[3,2] = el(r(table),7,1)+1;

matrix list TableA1;
putexcel set "`figtabOutputFolder'\TableA1.xlsx", replace;
putexcel A1 = matrix(TableA1), names;

timer off 1;
timer list 1;
