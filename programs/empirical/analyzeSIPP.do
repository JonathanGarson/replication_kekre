/* analyzeSIPP.do 
   Uses 2004 SIPP panel to summarize EU-wealth and EU-earnings relationships */

#delimit ;

clear all;
set maxvar 10000;

timer clear 1;
timer on 1;

local sourceFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\data";
local dtaOutputFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\programs\empirical\output";
local figtabOutputFolder = "`dtaOutputFolder'\figures_and_tables";
local modelOutputFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\programs\model\output";

clear;
graph drop _all;

global makedata = 1;

if $makedata {; 

	/* Set filenames to load */;
	local lf1_3 = "sippl04puw3";
	local al1_3 = "sippp04putm3";
	local lf2a_3 = "sippl04puw5";
	local lf2b_3 = "sippl04puw6";
	local lf2c_3 = "sippl04puw7";

	local lf1_6 = "sippl04puw6";
	local al1_6 = "sippp04putm6";
	local lf2a_6 = "sippl04puw8";
	local lf2b_6 = "sippl04puw9";
	local lf2c_6 = "sippl04puw10";

	foreach module in 3 6 {;
		cd "`sourceFolder'\2004_SIPP";
		
		/* Load initial labor force status and own earnings */;
		use `lf1_`module'', replace;
		keep ssuseq ssuid srefmon rhcalyr rhcalmn shhadid eoutcome epppnum errp wpfinwgt rwkesr4 tpearn;
		rename rwkesr4 rwkesr_init;
		/* I only care about the fourth reference month, and for simplicity assume that
		 4th week summarizes employment status at the end of this month, to match 
		 with assets reported in topical module */;
		drop if srefmon < 4;
		/* Further focus only reference person in each household, with aim of focusing 
		 on household head */;
		keep if errp == 1 | errp == 2;
		sort ssuid shhadid eoutcome epppnum;
		save "`dtaOutputFolder'\2004_SIPP_`module'_1.dta", replace;

		/* Prep al1 to merge */;
		clear;
		use `al1_`module'', replace;
		keep ssuid shhadid eoutcome epppnum thhtnw;
		sort ssuid shhadid eoutcome epppnum;
		save "`dtaOutputFolder'\2004_SIPP_`module'_2.dta", replace;

		/* Merge */;
		use "`dtaOutputFolder'\2004_SIPP_`module'_1.dta", replace;
		merge 1:1 ssuid shhadid eoutcome epppnum using "`dtaOutputFolder'\2004_SIPP_`module'_2.dta", keep(1 3);
		drop _merge;
		save "`dtaOutputFolder'\2004_SIPP_`module'_3.dta", replace;

		/* Now merge in labor force status one year ahead from lf2a */;
		merge 1:1 ssuid shhadid eoutcome epppnum rhcalmn using `lf2a_`module'', keep(1 3) keepusing(rwkesr4);
		rename rwkesr4 rwkesr_fin_a;
		rename _merge _merge_a;

		/* Now merge in labor force status one year ahead from lf2b */;
		merge 1:1 ssuid shhadid eoutcome epppnum rhcalmn using `lf2b_`module'', keep(1 3) keepusing(rwkesr4);
		rename rwkesr4 rwkesr_fin_b;
		rename _merge _merge_b;

		/* Now merge in labor force status one year ahead from lf2c */;
		merge 1:1 ssuid shhadid eoutcome epppnum rhcalmn using `lf2c_`module'', keep(1 3) keepusing(rwkesr4);
		rename rwkesr4 rwkesr_fin_c;
		rename _merge _merge_c;

		/* Define final labor force status */;
		gen rwkesr_fin = rwkesr_fin_a if _merge_a == 3;
		replace rwkesr_fin = rwkesr_fin_b if _merge_b == 3;
		replace rwkesr_fin = rwkesr_fin_c if _merge_c == 3;

		cd "`dtaOutputFolder'";
		save 2004_SIPP_`module'.dta, replace;
		
		erase 2004_SIPP_`module'_1.dta;
		erase 2004_SIPP_`module'_2.dta;
		erase 2004_SIPP_`module'_3.dta;
	};

	clear;
	use 2004_SIPP_3.dta, replace;
	append using 2004_SIPP_6.dta;
	save 2004_SIPP.dta, replace;
};

use "`dtaOutputFolder'\2004_SIPP.dta", replace;

/* Define initial and final employment status */;
gen empl_init = 1 if rwkesr_init == 1 | rwkesr_init == 2;
replace empl_init = 0 if rwkesr_init == 3 | rwkesr_init == 4; 
gen unem_fin = 1 if rwkesr_fin == 3 | rwkesr_fin == 4;
replace unem_fin = 0 if rwkesr_fin == 1 | rwkesr_fin == 2;

gen logearn = log(tpearn);

/* Based on binned scatter, it appears that bottom and top 5% of observations are considerable outliers.  
 Hence, I define sample to exclude these.  Further exclude missing observations of logearn. */;
 
quietly summarize thhtnw if ~missing(unem_fin) & ~missing(thhtnw) & ~missing(logearn) & empl_init == 1 & eoutcome == 201, detail;
gen samplen = inrange(thhtnw, r(p5), r(p95)) & ~missing(unem_fin) & ~missing(thhtnw) & ~missing(logearn) & empl_init == 1 & eoutcome == 201;

/* Finally scale by average monthly household income of $6761 from 2004 SCF */;
scalar meanInc = 6761;
egen networthMean = mean(thhtnw) if samplen == 1;
gen networthS = (thhtnw - networthMean)/meanInc;

/* Define time FE */;

gen yyyymm = rhcalyr*100+rhcalmn;
tab yyyymm, gen(ym);
drop ym1 ;

/* Run EU regressions of interest */;

eststo: regress unem_fin ym* networthS if samplen == 1 [aw=wpfinwgt], cluster(ssuid);
eststo: regress unem_fin ym* networthS logearn if samplen == 1 [aw=wpfinwgt], cluster(ssuid);

esttab using "`figtabOutputFolder'\TableA3_2.csv", b(a2) se(a1) nostar keep(networthS logearn) replace;
eststo clear;

timer off 1;
timer list 1;
