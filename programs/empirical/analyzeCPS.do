/* analyzeCPS.do 
   Uses 2004-2007 CPS to estimate EU flows by weekly pay */
 
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
	/* Merges successive ORG samples from CEPR */;

	local first=2004;
	while `first' <=2006  {;
		local second = `first' + 1;
		
		cd "`sourceFolder'\2004_2007_CPS";
		use cepr_org_`second'.dta;
		keep hhid hhid2 lineno fnlwgt age female wbho empl unem nilf;

		tostring hhid2, replace;
		replace hhid2 = "" if hhid2 == ".";
		recast str5 hhid2;

		rename fnlwgt fnlwgt2;
		rename empl empl2;
		rename unem unem2;
		rename nilf nilf2;
		gen age1 = age - 1;
		rename age age2;
		sort hhid hhid2 lineno wbho female age1;
		quietly by hhid hhid2 lineno wbho female age1: gen dup = cond(_N==1,0,_n);
			*** Flag any duplicate records and delete.;
		drop if dup > 0;
		sort hhid hhid2 lineno wbho female age1;
		cd "`dtaOutputFolder'";
		save second.dta, replace;

		clear;
		cd "`sourceFolder'\2004_2007_CPS";
		use cepr_org_`first'.dta;
		keep year month minsamp hhid hhid2 lineno fnlwgt age female wbho empl unem nilf weekpay;

		tostring hhid2, replace;
		replace hhid2 = "" if hhid2 == ".";
		recast str5 hhid2;

		rename fnlwgt fnlwgt1;
		rename empl empl1;
		rename unem unem1;
		rename nilf nilf1;
		rename age age1;
		sort hhid hhid2 lineno wbho female age1;
		quietly by hhid hhid2 lineno wbho female age1: gen dup = cond(_N==1,0,_n);
			*** Flag any duplicate records and delete.;
		drop if dup > 0;

		cd "`dtaOutputFolder'";
		sort hhid hhid2 lineno wbho female age1;
		merge 1:1 hhid hhid2 lineno wbho female age1 using second.dta;
		   *** This creates a record of all individuals, matched if they agree on hh, line, race, sex, and exactly on age. And also neither record can be duplicated.;
		save `first'_CPS.dta, replace;
		erase second.dta;

		local first = `second';
	};
	
	local startYear = 2004;
	local endYear = 2006;
	
	cd "`dtaOutputFolder'";

	use `startYear'_CPS.dta, replace;
	local year=`startYear'+1;

	while `year'<=`endYear' {;
		quietly append using `year'_CPS.dta;
		local ++year;
	};
	
	save `startYear'_`endYear'_CPS.dta, replace;
};

cd "`dtaOutputFolder'";
use `startYear'_`endYear'_CPS.dta, replace;

quietly gen logweekpayRaw = log(weekpay);
quietly egen logweekpayMean = mean(logweekpayRaw) if ~missing(unem2) & ~missing(weekpay) & empl1==1 & nilf2~=1, by(year month);
quietly gen logweekpay = logweekpayRaw - logweekpayMean;

quietly gen yyyymm = year*100+month;
quietly tab yyyymm, gen(ym);
drop ym1 ; 

quietly gen uniqueHH = hhid+hhid2;
eststo: regress unem2 ym* logweekpay if empl1==1 & nilf2~=1 [aw=fnlwgt1], cluster(uniqueHH);
esttab using "`figtabOutputFolder'\TableA3_1.csv", b(a2) se(a1) nostar keep(logweekpay) replace;
eststo clear;

timer off 1;
timer list 1;
