/* estimateSepRate.do 
   Estimates separation rate following methodology of Shimer (2012) */

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
	* BLS, short-term unemployment;
	import excel "`sourceFolder'\sep_rate\BLS_unemployed_lt5wk.xlsx", sheet("LNS13008396") cellrange(A19:D901) case(lower) firstrow clear;
	destring year, replace;
	gen mon = substr(period,2,2);
	destring mon, gen(month);
	gen ym = ym(year, month);
	format ym %tm;
	rename observationvalue short_ue;
	drop year month period mon label;
	order ym;
	save "`dtaOutputFolder'\seprate_short_ue.dta", replace;

	* BLS, unemployment;
	import excel "`sourceFolder'\sep_rate\BLS_unemployed.xlsx", sheet("BLS Data Series") cellrange(A12:M86) case(lower) firstrow clear;
	* reshape dataset from wide to long;
	ds;
	local vars `r(varlist)';
	local year "year";
	local varlist: list vars - year;
	foreach var of local varlist{;
		rename `var' m_`var';
	};
	reshape long m_, i(year) j(mon, string);
	rename m_ unemploy;
	* generate year_month variable;
	gen month = month(date(mon, "M"));
	drop mon;
	gen ym = ym(year, month);
	format ym %tm;
	sort ym;
	order ym unemploy;
	drop if missing(unemploy);
	drop year month;
	save "`dtaOutputFolder'\seprate_unemploy.dta", replace;

	* BLS, employment;
	import excel "`sourceFolder'\sep_rate\BLS_employed.xlsx", sheet("BLS Data Series") cellrange(A12:M86) case(lower) firstrow clear;
	* reshape dataset from wide to long;
	ds;
	local vars `r(varlist)';
	local year "year";
	local varlist: list vars - year;
	foreach var of local varlist{;
		rename `var' m_`var';
	};
	reshape long m_, i(year) j(mon, string);
	rename m_ employ;
	* generate year_month variable;
	gen month = month(date(mon, "M"));
	drop mon;
	gen ym = ym(year, month);
	format ym %tm;
	sort ym;
	drop year month;
	order ym employ;
	drop if missing(employ);
	save "`dtaOutputFolder'\seprate_employ.dta", replace;

	/* Fraction of short-term unemployed workers among all unemployed workers 
	   from CPS Basic Monthly Files, seasonally adjusted using Census' X-13ARIMA-SEATS. */
	use "`sourceFolder'\sep_rate\CPS_basic_final.dta", replace;
	rename date_ym ym;
	drop year month;
	save "`dtaOutputFolder'\seprate_cpsb.dta", replace;

	* Merge;
	use "`dtaOutputFolder'\seprate_short_ue.dta", replace;
	merge 1:1 ym using "`dtaOutputFolder'\seprate_unemploy.dta", nogen;
	merge 1:1 ym using "`dtaOutputFolder'\seprate_employ.dta", nogen;
	merge 1:1 ym using "`dtaOutputFolder'\seprate_cpsb.dta", nogen;
	tsset ym;
	save "`dtaOutputFolder'\seprate.dta", replace;
};

use "`dtaOutputFolder'\seprate.dta", replace;

* Construct short-term unemployment series which corrects for 94 structural break;
gen short_ue_adj = short_ue if ym <= tm(1993m12);
replace short_ue_adj = unemploy * short_ue_frac_adj if ym >= tm(1994m1);

gen year = year(dofm(ym));
gen month = month(dofm(ym));
order year month ym;

* Output to csv to analyze in Matlab;
export delimited using "`dtaOutputFolder'\seprate_for_matlab.csv", replace;

timer off 1;
timer list 1;
