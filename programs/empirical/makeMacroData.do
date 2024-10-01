/* makeMacroData.do 
   Produces macro time series from variety of sources to compare vs. model */

#delimit ;

clear all;
set maxvar 10000;

timer clear 1;
timer on 1;

local sourceFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\data";
local dtaOutputFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\programs\empirical\output";
local figtabOutputFolder = "`dtaOutputFolder'\figures_and_tables";
local modelOutputFolder = "C:\Users\rkekre\Dropbox\UI in Macro Stabilization\Journals\Restud replication package\data_and_programs\programs\model\output";

* BLS, unemployment rate;
import excel "`sourceFolder'\macro_ts\LNS14000000.xlsx", sheet("LNS14000000") cellrange(A12:M35) case(lower) firstrow clear;
rename jan unemp1;
rename feb unemp2;
rename mar unemp3;
rename apr unemp4;
rename may unemp5;
rename jun unemp6;
rename jul unemp7;
rename aug unemp8;
rename sep unemp9;
rename oct unemp10;
rename nov unemp11;
rename dec unemp12;
reshape long unemp, i(year) j(month);
gen ym = ym(year,month);
format ym %tm;
keep ym unemp;
order ym unemp;
save "`dtaOutputFolder'\macro_unemp.dta", replace;

* BLS, fraction long-term unemployed;
import excel "`sourceFolder'\macro_ts\LNS13025703.xlsx", sheet("LNS13025703") cellrange(A14:M37) case(lower) firstrow clear;
rename jan ltu1;
rename feb ltu2;
rename mar ltu3;
rename apr ltu4;
rename may ltu5;
rename jun ltu6;
rename jul ltu7;
rename aug ltu8;
rename sep ltu9;
rename oct ltu10;
rename nov ltu11;
rename dec ltu12;
reshape long ltu, i(year) j(month);
gen ym = ym(year,month);
format ym %tm;
keep ym ltu;
order ym ltu;
save "`dtaOutputFolder'\macro_ltu.dta", replace;

* FRB, Fed Funds rate;
import excel "`sourceFolder'\macro_ts\FEDFUNDS.xls", sheet("FRED Graph") cellrange(A11:B775) case(lower) firstrow clear;
gen ym = mofd(observation_date);
format ym %tm;
order ym fedfunds;
drop observation_date;
save "`dtaOutputFolder'\macro_fedfunds.dta", replace;

* BEA, core PCE deflator;
import excel "`sourceFolder'\macro_ts\PCEPILFE.xls", sheet("FRED Graph") cellrange(A11:B762) case(lower) firstrow clear;
gen ym = mofd(observation_date);
format ym %tm;
rename pcepilfe corepce;
order ym corepce;
drop observation_date;
save "`dtaOutputFolder'\macro_corepce.dta", replace;

* BLS, nominal wage index;
import excel "`sourceFolder'\macro_ts\CES0500000003.xls", sheet("FRED Graph") cellrange(A11:B180) case(lower) firstrow clear;
gen ym = mofd(observation_date);
format ym %tm;
rename ces0500000003 nomwage;
order ym nomwage;
drop observation_date;
save "`dtaOutputFolder'\macro_nomwage.dta", replace;

* BLS, vacancies;
import excel "`sourceFolder'\macro_ts\JTSJOL.xls", sheet("FRED Graph") cellrange(A11:B242) case(lower) firstrow clear;
gen ym = mofd(observation_date);
format ym %tm;
rename jtsjol vacancies;
order ym vacancies;
drop observation_date;
save "`dtaOutputFolder'\macro_vacancies.dta", replace;

* BLS, unemployed;
import excel "`sourceFolder'\macro_ts\UNEMPLOY.xls", sheet("FRED Graph") cellrange(A11:B878) case(lower) firstrow clear;
gen ym = mofd(observation_date);
format ym %tm;
rename unemploy unemployed;
order ym unemployed;
drop observation_date;
save "`dtaOutputFolder'\macro_unemployed.dta", replace;

* BEA, nominal GDP, consumption nondurables, and consumption services;
import delim using "`sourceFolder'\macro_ts\US_GDP.csv", rowrange(7:12) colrange(3:202) clear;
xpose, clear;
gen year = floor((_n-1)/4)+1970;
gen month = 3*(mod(_n-1,4)+1);
rename v1 nomgdp;
rename v5 cndur;
rename v6 cserv;
gen ym = ym(year,month);
format ym %tm;
keep ym nomgdp cndur cserv;
order ym nomgdp cndur cserv;
save "`dtaOutputFolder'\macro_nom.dta", replace;

* BEA, real GDP;
import delim using "`sourceFolder'\macro_ts\US_Real_GDP.csv", rowrange(7:7) colrange(3:202) clear;
xpose, clear;
gen year = floor((_n-1)/4)+1970;
gen month = 3*(mod(_n-1,4)+1);
rename v1 realgdp;
gen ym = ym(year,month);
format ym %tm;
keep ym realgdp;
order ym realgdp;
save "`dtaOutputFolder'\macro_real.dta", replace;

* BLS, civilian population over 16;
import excel "`sourceFolder'\macro_ts\US_Civilian_Pop.xlsx", sheet("BLS Data Series") cellrange(B12:D212) case(lower) firstrow clear;
gen quarter = substr(period,3,1);
destring quarter, replace;
gen month = 3*quarter;
gen ym = ym(year,month);
format ym %tm;
rename value pop;
keep ym pop;
order ym pop;
save "`dtaOutputFolder'\macro_pop.dta", replace;

* Merge datasets;
use "`dtaOutputFolder'\macro_unemp.dta", replace;
merge 1:1 ym using "`dtaOutputFolder'\macro_ltu.dta", nogen;
merge 1:1 ym using "`dtaOutputFolder'\macro_fedfunds.dta", nogen;
merge 1:1 ym using "`dtaOutputFolder'\macro_corepce.dta", nogen;
merge 1:1 ym using "`dtaOutputFolder'\macro_nomwage.dta", nogen;
merge 1:1 ym using "`dtaOutputFolder'\macro_vacancies.dta", nogen;
merge 1:1 ym using "`dtaOutputFolder'\macro_unemployed.dta", nogen;
merge 1:1 ym using "`dtaOutputFolder'\macro_nom.dta", nogen;
merge 1:1 ym using "`dtaOutputFolder'\macro_real.dta", nogen;
merge 1:1 ym using "`dtaOutputFolder'\macro_pop.dta", nogen;
merge 1:1 ym using "`dtaOutputFolder'\macro_ui_weeks_fv.dta", nogen;
tsset ym;

* Compute new variables;
qui gen conspc = (cndur+cserv)/(pop*nomgdp/realgdp);
qui gen tightness = vacancies/unemployed;

* Compute average annual growth rates from 1990-2019;
qui gen corepce_growth = (corepce/L.corepce)^12-1;
qui gen nomwage_growth = (nomwage/L.nomwage)^12-1;
qui gen conspc_growth = (conspc/L3.conspc)^4-1;
tabstat corepce_growth nomwage_growth conspc_growth if ym >= ym(1990,1) & ym <= ym(2019,12), statistics(mean);
*  implies trend growth of 1.9%, 2.6%, and 1.7%, respectively;

* Interpolate real consumption per capita from quarterly to monthly;
qui ipolate conspc ym, gen(conspc_monthly);

* Compute series vs. April 2008 (detrended as needed);
gen long obsn = _n;
summ obsn if ym == ym(2008,4), meanonly;
local base = obsn[`r(min)'];
qui gen dunemp = (unemp - unemp[`base'])/100;
qui gen dltu = (ltu - ltu[`base'])/100;
qui gen dfedfunds = (fedfunds - fedfunds[`base'])/100;
qui gen dcorepce = (corepce/corepce[`base']-1) - (1.019^((ym-ym(2008,4))/12)-1);
qui gen dnomwage = (nomwage/nomwage[`base']-1) - (1.026^((ym-ym(2008,4))/12)-1);
qui gen dconspc = (conspc_monthly/conspc_monthly[`base']-1) - (1.017^((ym-ym(2008,4))/12)-1);
qui gen dtightness = tightness/tightness[`base']-1;
keep if ym >= ym(2008,1) & ym <= ym(2015,12);
keep ym dunemp dltu dfedfunds dcorepce dnomwage dconspc dtightness;

* Output to csv;
export delimited using "`dtaOutputFolder'\macrodata.csv", replace;
export delimited using "`modelOutputFolder'\macrodata.csv", replace;

* Transfer unemployment rate to csv in model folder;
keep if ym >= ym(2008,5) & ym <= ym(2014,12);
keep dunemp;
export delimited using "`modelOutputFolder'\ur05081214.csv", novarnames datafmt replace;

* Transfer Farber-Valletta median weeks of UI to csv in model folder;
use "`sourceFolder'\macro_ts\ui_weeks_fv.dta", replace;
keep ui_weeks_fv;
export delimited using "`modelOutputFolder'\ui_weeks_fv.csv", novarnames replace;

* Transfer estimated separation rate to csv in model folder;
import delim using "`dtaOutputFolder'\seprate.csv", asdouble clear;
export delimited using "`modelOutputFolder'\seprate.csv", novarnames replace;
rename v1 year;
rename v2 month;
rename v5 seprate_inno;
gen ym = ym(year,month);
keep if ym >= ym(2008,5) & ym <= ym(2014,12);
keep seprate_inno;
export delimited using "`modelOutputFolder'\seprate05081214.csv", novarnames replace;

timer off 1;
timer list 1;
