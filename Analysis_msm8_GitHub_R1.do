

*****************************************************
**#Cleaning data************************************* 
*****************************************************
*Original database - incident T2DM based on code only
cls
clear all
use ".../predmdb2.dta"
keep patid patid t2d_date t2d_censoring t2d_outcome
renames t2d_date t2d_censoring t2d_outcome, suffix(_R0)
tempfile dbr0
save `dbr0', replace 

*Updated database - incident T2DM based on earliest of code + HbA1c + FPG + meds
*The variable t2d_date is the min of t2d_date_a1c t2d_date_code t2d_date_drug t2d_date_fgt, and t2d_censoring
*The variable t2d_outcome is 0/1, if missing on all four variables (t2d_date_a1c t2d_date_code t2d_date_drug t2d_date_fgt) is 0, otherwise 1.
*The variable t2d_censoring is identical for all patients, which is censoring date.
cls
global off ".../Analysis"
cap use "$off/predmdb2.dta", clear
drop smk_date weight w_date monthday pracid yob hba1c fgt type ethnic hba1c_date fgt_date
order patid age_index gender bmi imd2010_5 white smk_status 			/*
	  */dob index 								/*
	  */t2d_date t2d_outcome t2d_censoring					/*
	  */cancersite cancer_date cancer cancer_censoring 			/*
	  */cause dodeath death death_censoring					/*
	  */t2d_date_a1c t2d_date_fgt t2d_date_drug t2d_date_code
renames age_index index cancersite cancer cause dodeath death smk_status \ ageind indexd cancer_site cancer_outcome death_cause death_date death_outcome smk
renames t2d_date t2d_censoring t2d_outcome, suffix(_R1)
describe, short 
distinct patid
mdesc
sum

merge 1:1 patid using `dbr0', nogen
tab t2d_outcome_R0 t2d_outcome_R1, m col
drop t2d_date_a1c-t2d_outcome_R0
renames t2d_date_R1 t2d_outcome_R1 t2d_censoring_R1 \ t2d_date t2d_outcome t2d_censoring

***************************************************************
*For Death until 14/01/2019 - For nonabsorbing until 30/11/2018
*Therefore you can only die but not develop a non-absorbing
*All times stopped at 30/11/2018
***************************************************************
replace death_censoring = date("30/11/2018", "DMY")
replace death_outcome = 0 if death_date>death_censoring & death_outcome == 1
replace death_date = death_censoring if death_outcome == 0
replace death_cause = "" if death_outcome == 0
mdesc death_cause if death_outcome == 1					
drop if death_cause == "" & death_outcome == 1					/*38 individuals DROPPED*/
distinct patid
replace cancer_site = "" if cancer_outcome == 0					/*Original databate already right-censored at 30/11/2018 for outcome but not removed cancer ICD*/

renames death_outcome death_cause \ survi surv_cause
gen surv = (death_date-indexd)/365.24
drop death_censoring

renames cancer_outcome \ canceri
gen cancer = cond(canceri, (cancer_date-indexd)/365.24, (death_date-indexd)/365.24) 
drop cancer_date cancer_censoring

renames t2d_outcome \ t2di
gen t2d = cond(t2di, (t2d_date-indexd)/365.24, (death_date-indexd)/365.24) 
drop t2d_date t2d_censoring

order t2d t2di cancer canceri cancer_site surv survi surv_cause, last

*********************************************************			
*19 individuals t=zero because diagnosis=last day in HES*
*For T2D identified until 30/11/2018
*For Cancer identified until 30/11/2018
*For Death until 14/01/2019 (amended to 30/11/2018)
*********************************************************
describe, short 
drop if cancer == 0 | t2d == 0 | surv == 0					/*19 individuals DROPPED*/
describe, short

sum cancer t2d surv										

*t2d vs cancer*
gen ddc = t2d - cancer
drop if ddc == 0 & t2di == 1 & canceri == 1					/*408 individuals DROPPED*/

*t2d vs survival
gen dds = t2d - surv
drop if dds>0									/*2 individuals DROPPED*/

*cancer vs survival
gen dcs = cancer - surv
drop if dcs>0									/*1 individuals DROPPED*/

distinct patid

*Same times*
count if dds == 0 & t2di == 1 & survi == 1
replace t2d = surv - (0.5/365.24) if dds == 0 & t2di == 1 & survi == 1		/*20 people t2d=death*/
drop dds

count if dcs == 0 & canceri == 1 & survi == 1					/*74 people cancer=death*/
replace cancer = surv - (0.5/365.24) if dcs == 0 & canceri == 1 & survi == 1

gen cdeath = substr(surv_cause,1,1)
tab cdeath canceri if cdeath == "C", row missing				/*1662 individuals with cancer death but not cancer event*/
gen flag = 1 if cdeath == "C" & canceri == 0
replace cancer = surv - (0.5/365.24) if flag == 1
replace canceri = 1 if flag == 1
replace cancer_site = surv_cause if flag == 1					/*For these individuals used the same ICD of death for incident cancer*/
drop dcs cdeath flag

*t2d after cancer*
count if ddc>0 & t2di == 1 & canceri == 1
replace t2d = surv if ddc>0 & t2di == 1 & canceri == 1				
replace t2di = 0   if ddc>0 & t2di == 1 & canceri == 1							
drop ddc

describe, short 

count if bmi <18.5
egen float bmig = cut(bmi), at(0 25 30 35 40 100) icodes label
label define bmig 0 `"<25"', modify
label define bmig 1 `"≥25 to <30"', modify
label define bmig 2 `"≥30 to <35"', modify
label define bmig 3 `"≥35 to <40"', modify
label define bmig 4 `"≥40"', modify
label variable bmi  `"Body mass index (kg/m2)"'
label variable bmig `"Body mass index category (kg/m2)"'

egen float ageg = cut(ageind), at(0 55 65 75 110) icodes label
label define ageg 0 `"<55"', modify
label define ageg 1 `"≥55 to <65"', modify
label define ageg 2 `"≥65 to <75"', modify
label define ageg 3 `"≥75"', modify
label variable ageind `"Age at prediabetes diagnosis (years)"'
label variable ageg   `"Age at prediabetes diagnosis group (years)"'

tostring white, replace
replace white = "White" if white == "1"
replace white = "Non-White" if white == "0"
replace white = "" if white == "."
rename white ethnic
codebook ethnic
label variable ethnic `"Ethnicity"'

groups imd2010_5, missing
label variable imd2010_5 `"IMD_2010"'

sum surv cancer t2d
mdesc
distinct patid
describe, short
cap save "$off/Datasets/db2.dta", replace 


**********************************************************
**#Set the data for Multistate for transition graph/rates*
**********************************************************
cls
clear all
global off ".../Analysis"
cap use "$off/Datasets/db2.dta", replace

*Non absorbing states*
gen T2D = t2d if t2d<surv & t2d<cancer & t2di == 1
gen T2Di = cond(T2D, 1, 0, 0)
replace T2D = surv if T2D == .
tab T2Di t2di, m
tab T2Di survi, m
drop t2d t2di

gen Kp = cancer if cancer<surv & cancer<T2D & T2Di == 0
gen Kpi = cond(Kp, 1, 0, 0)
replace Kp = surv if Kp == .
tab Kpi canceri, m
tab Kpi survi, m

gen Kt = cancer if cancer<surv & cancer>T2D & T2Di == 1
gen Kti = cond(Kt, 1, 0, 0)
replace Kt = surv if Kt == .
tab Kti canceri, m
tab Kti survi, m

*Absorbing states*
foreach nm in Dp DKp DKt Dt {
	gen `nm' = surv
	gen `nm'i = survi
}

mdesc
groups survi canceri T2Di, sepby(survi) missing

matrix tmat = (.,1,.,2,3,.,.,.\.,.,4,.,.,.,.,.\.,.,.,.,.,.,.,.\.,.,.,.,.,.,.,.\.,.,.,.,.,5,.,6\.,.,.,.,.,.,7,.\.,.,.,.,.,.,.,.\.,.,.,.,.,.,.,.)
mat colnames tmat = "pT2D" "Kp" "DKp" "Dp" "T2D" "Kt" "DKt" "Dt"
mat rownames tmat = "pT2D" "Kp" "DKp" "Dp" "T2D" "Kt" "DKt" "Dt"
mat list tmat

order Kp Kpi DKp DKpi Dp Dpi T2D T2Di Kt Kti DKt DKti Dt Dti, after(ageg)
msset, id(patid) states(Kpi DKpi Dpi T2Di Kti DKti Dti) times(Kp DKp Dp T2D Kt DKt Dt) transmatrix(tmat)
matrix list r(freqmatrix)
matrix ntmat = r(freqmatrix)
mat colnames ntmat = "pT2D" "Kp" "DKp" "Dp" "T2D" "Kt" "DKt" "Dt"
mat rownames ntmat = "pT2D" "Kp" "DKp" "Dp" "T2D" "Kt" "DKt" "Dt"
matrix list ntmat

cap save "$off/Datasets/multistate_structure.dta", replace 


****************************************************
**#Set the data for Pseudocounts********************
****************************************************
cls
clear all
global off ".../Analysis"
cap use "$off/Datasets/db2.dta", replace  

*-------------------------------*
*This is as above for multistate*
*-------------------------------*

*Non absorbing states*
gen T2D = t2d if t2d<surv & t2d<cancer & t2di == 1
gen T2Di = cond(T2D, 1, 0, 0)
replace T2D = surv if T2D == .
tab T2Di t2di, m
tab T2Di survi, m
drop t2d t2di

gen Kp = cancer if cancer<surv & cancer<T2D & T2Di == 0
gen Kpi = cond(Kp, 1, 0, 0)
replace Kp = surv if Kp == .
tab Kpi canceri, m
tab Kpi survi, m

gen Kt = cancer if cancer<surv & cancer>T2D & T2Di == 1
gen Kti = cond(Kt, 1, 0, 0)
replace Kt = surv if Kt == .
tab Kti canceri, m
tab Kti survi, m

*Absorbing states*
foreach nm in Dp DKp DKt Dt {
	gen `nm' = surv
	gen `nm'i = survi
}

*-------------------------------*
*Below is for pseudocounts------*
*-------------------------------*
stset DKp, failure(DKpi==1) id(patid)						/*Any among Dp DKp DKt Dt can be used*/

foreach var of varlist Kp DKp Dp T2D Kt DKt Dt {
	stsplit post`var' if `var'i == 1, at(0) after(`var')
	recode post`var' (0 = 1) (else = 0)
}

label define statelbl 1 "pT2D" 2 "Kp" 3 "DKp" 4 "Dp" 5 "T2D" 6 "Kt" 7 "DKt" 8 "Dt"

groups post*

*Starting and non-absorbing*
generate fromstate = 1
replace  fromstate = 2 if postKp
replace  fromstate = 5 if postT2D
replace  fromstate = 6 if postT2D & postKt

*Absorbing*
generate tostate = fromstate
replace tostate = 3 if DKpi == 1 & fromstate == 2								
replace tostate = 7 if DKti == 1 & fromstate == 6
replace tostate = 4 if Dpi  == 1 & fromstate == 1
replace tostate = 8 if Dti  == 1 & fromstate == 5
by patid (_t), sort: replace tostate = fromstate[_n + 1] if _n < _N

label values fromstate tostate statelbl
groups fromstate tostate if gender == 1, sepby(fromstate)
groups fromstate tostate if gender == 2, sepby(fromstate)

cap save "$off/Datasets/pseudocounts_structure.dta", replace  


********************************************
**#Generate pseudocounts SPLIT CPUs*********
********************************************
global off ".../Analysis"

clear
use "$off/Datasets/pseudocounts_structure.dta", clear
stpmstate ps_p1_=p(1) ps_p2_=p(2) ps_p3_=p(3) ps_p4_=p(4) ps_p5_=p(5) ps_p6_=p(6) ps_p7_=p(7) ps_p8_=p(8) 					///
		, by(gender ageg) at(0(0.25)10) from(fromstate) to(tostate) id(patid) atnumbers replace
save "$off/Datasets/ps_probs2.dta", replace

clear
use "$off/Datasets/pseudocounts_structure.dta", clear
stpmstate ps_los1_=los(1) ps_los2_=los(2) ps_los3_=los(3) ps_los4_=los(4) ps_los5_=los(5) ps_los6_=los(6) ps_los7_=los(7) ps_los8_=los(8)	///
		, by(gender ageg) at(0(0.25)10) from(fromstate) to(tostate) id(patid) atnumbers replace
save "$off/Datasets/ps_los2.dta", replace


******************************************************************
**#GLM estimates from pseudocounts -- ALL STATES -- SPLIT CPUs****
******************************************************************

*---------------------------------------------*
***Probs - OVERALL POPULATION [DESKTOP 1]-----*
*---------------------------------------------*
cls
clear all
global off ".../Analysis"
cap use "$off/Datasets/ps_probs2.dta", replace

rename ps_p?__25 ps_p?_0_25 
rename ps_p?__5  ps_p?_0_5 
rename ps_p?__75 ps_p?_0_75 

keep patid gender ageg ps_p*
egen float nmiss = rowmiss(ps_p*)
tab nmiss, m
drop if nmiss != 0
drop nmiss

cd "$off/Results"
forvalues g = 1/2 {																
	foreach t of varlist ps_p* {
		qui cap glm `t' ib(frequent).ageg if gender == `g', eform link(log) vce(robust)
		
		if _rc != 0 {							/*If there are no estimates*/
			preserve
			clear
			qui set obs 4
			qui gen estimate = . 
			qui gen stderr   = .
			qui gen min95    = .
			qui gen max95    = .
			qui gen ageg     = _n
			qui save "po_g`g'_`t'", replace
			restore
		}
	
		else {
			qui cap margins ageg, saving("po_g`g'_`t'", replace)
		}

		di "Error == " _rc " || Sex == `g' -- Var == `t' || $S_TIME"
	}
}

clear
tempfile respo
save `respo', emptyok replace
local files: dir "$off/Results/" files "po_g*.dta"
foreach file in `files' {
	qui use "$off/Results/`file'", clear
	qui gen filen = "`file'"
	qui append using `respo'
	qui save `respo', replace
}
use `respo', clear
replace _m1 = ageg-1 if ageg != .
keep _margin _se_margin _ci_lb-filen
split filen, p("_")
drop filen3
replace filen1 = "Probs"
replace filen2 = "Women" if filen2 == "g2"
replace filen2 = "Men" if filen2 == "g1"
replace filen5 = subinstr(filen5, ".dta", "", .)
rename filen5 t1
replace filen6 = subinstr(filen6, ".dta", "", .)
rename filen6 t2
replace filen4 = subinstr(filen4, "p", "", .)
rename filen4 state
gen time = t1 + "." + t2
destring time, replace
destring state, replace
drop t1 t2
renames filen1 filen2 _m1 \ metric sex ageg
order metric sex state time ageg, first
sencode sex, gsort(-sex) replace
sort sex state time ageg
foreach var of varlist _margin _ci_lb _ci_ub {
	replace `var' = 1 if time == 0 & state == 1
	replace `var' = 0 if time == 0 & state != 1
}
renames _margin _se_margin \ prob se_prob
mdesc
list if _ci_lb == . | _ci_ub == . | prob == ., sepby(filen)
replace _ci_lb = 0 if _ci_lb == . & prob != .
replace _ci_ub = 0 if _ci_ub == . & prob != .
foreach var of varlist prob _ci_lb _ci_ub {
	replace `var' = 0 if `var' == .
}
mdesc
list if se_prob == ., sepby(filen)
sort sex state time ageg 
tab time state if sex == 1
tab time state if sex == 2
save "$off/Results/States_Overall_probs2.dta", replace

local files: dir "$off/Results/" files "po_g*.dta"
foreach file in `files' {
	qui erase "$off/Results/`file'"
}


*---------------------------------------------*
***LOS - OVERALL POPULATION  [HPC 1]----------*
*---------------------------------------------*
cls
clear all
global off ".../Analysis"
cap use "$off/Datasets/ps_los2.dta", replace

rename ps_los?__25 ps_los?_0_25 
rename ps_los?__5  ps_los?_0_5 
rename ps_los?__75 ps_los?_0_75 

keep patid gender ageg ps_los*
egen float nmiss = rowmiss(ps_los*)
tab nmiss, m
drop if nmiss != 0
drop nmiss

cd "$off/Results"
forvalues g = 1/2 {																
	foreach t of varlist ps_los* {
		qui cap glm `t' ib(frequent).ageg if gender == `g', link(identity) vce(robust)
		
		if _rc != 0 {							/*If there are no estimates*/
			preserve
			clear
			qui set obs 4
			qui gen estimate = . 
			qui gen stderr   = .
			qui gen min95    = .
			qui gen max95    = .
			qui gen ageg     = _n
			qui save "lo_g`g'_`t'", replace
			restore
		}
	
		else {
			qui cap margins ageg, saving("lo_g`g'_`t'", replace)
		}

		di "Error == " _rc " || Sex == `g' -- Var == `t' || $S_TIME"
	}
}

clear
tempfile reslo
save `reslo', emptyok replace
local files: dir "$off/Results/" files "lo_g*.dta"
foreach file in `files' {
	qui use "$off/Results/`file'", clear
	qui gen filen = "`file'"
	qui append using `reslo'
	qui save `reslo', replace
}
use `reslo', clear
replace _m1 = ageg-1 if ageg != .
keep _margin _se_margin _ci_lb-filen
split filen, p("_")
drop filen3
replace filen1 = "LOS"
replace filen2 = "Women" if filen2 == "g2"
replace filen2 = "Men"   if filen2 == "g1"
replace filen5 = subinstr(filen5, ".dta", "", .)
rename filen5 t1
replace filen6 = subinstr(filen6, ".dta", "", .)
rename filen6 t2
replace filen4 = subinstr(filen4, "los", "", .)
rename filen4 state
gen time = t1 + "." + t2
destring time, replace
destring state, replace
drop t1 t2
renames filen1 filen2 _m1 \ metric sex ageg
order metric sex state time ageg, first
sencode sex, gsort(-sex) replace
sort sex state time ageg 
foreach var of varlist _margin _ci_lb _ci_ub {
	replace `var' = 0 if time == 0
}
renames _margin _se_margin \ los se_los
mdesc
list if _ci_lb == . | _ci_ub == . | los == ., sepby(filen)
replace _ci_lb = 0 if _ci_lb == . & los != .
replace _ci_ub = 0 if _ci_ub == . & los != .
foreach var of varlist los _ci_lb _ci_ub {
	replace `var' = 0 if `var' == .
}
mdesc
list if se_los == ., sepby(filen)
sort sex state time ageg 
tab time state if sex == 1
tab time state if sex == 2
save "$off/Results/States_Overall_los2.dta", replace

local files: dir "$off/Results/" files "lo_g*.dta"
foreach file in `files' {
	qui erase "$off/Results/`file'"
}


*--------------------------------------------*
***Probs - RELATIVE [DEATH COMBINED] [HPC 2]-*
*--------------------------------------------*
cls
clear all
global off ".../Analysis"
use "$off/Datasets/ps_probs2.dta", clear

rename ps_p?__25 ps_p?_0_25 
rename ps_p?__5  ps_p?_0_5 
rename ps_p?__75 ps_p?_0_75 

egen float nmiss = rowmiss(ps_p1_0-ps_p8_10)
tab nmiss, m
drop if nmiss != 0
drop nmiss
sencode ethnic, replace
gen ps_p99_10 = ps_p4_10 + ps_p8_10 + ps_p7_10 + ps_p3_10			/*State all death combined*/
keep patid ageg gender bmig imd2010_5 ethnic smk ps_p*_10
sdecode smk, replace
sencode smk, gsort(-smk) gen(smkn)
drop smk
order smkn, after(ethnic)

preserve
clear
tempfile mprobs
save `mprobs', emptyok replace
restore

forvalues g = 1/2 {
	foreach s in 1 2 5 6 99 {						/*Non-absorbing + death*/
		qui cap glm ps_p`s'_10 i.ageg i.bmig i.smkn i.imd2010_5 i.ethnic if gender == `g', eform link(log) vce(robust) allbaselevels family(gaussian)
		qui local np = e(N)
		preserve
		qui parmest, fast
		qui gen sex = `g'
		qui gen state = `s'
		qui gen npart = `np'
		qui append using `mprobs'
		qui save `mprobs', replace		
		restore
		di "RelProbs -- Sex = `g' (N = `np') | State = `s' -- $S_TIME"
	}
}
qui use `mprobs', clear
qui gen time = "10 years"
qui gen out = "States_Relative_Probs"
qui drop z p
qui drop if parm == "_cons"
qui sencode parm, gen(sparm)
qui sort sex state sparm
compress
save "$off/Results/States_Relative_prob2.dta", replace


*------------------------------------------*
***LOS - RELATIVE [DEATH COMBINED] [HPC 3]-*
*------------------------------------------*
cls
clear all
global off ".../Analysis"
use "$off/Datasets/ps_los2.dta", clear

rename ps_los?__25 ps_los?_0_25 
rename ps_los?__5  ps_los?_0_5 
rename ps_los?__75 ps_los?_0_75 

egen float nmiss = rowmiss(ps_los1_0-ps_los8_10)
tab nmiss, m
drop if nmiss != 0
drop nmiss
sencode ethnic, replace
gen ps_los99_10 = ps_los4_10 + ps_los8_10 + ps_los7_10 + ps_los3_10		/*State all death combined*/
keep patid ageg gender bmig imd2010_5 ethnic smk ps_los*_10
sdecode smk, replace
sencode smk, gsort(-smk) gen(smkn)
drop smk
order smkn, after(ethnic)

preserve
clear
tempfile mlos
save `mlos', emptyok replace
restore

forvalues g = 1/2 {
	foreach s in 1 2 5 6 99 {						/*Non-absorbing + death*/
		qui cap glm ps_los`s'_10 i.ageg i.bmig i.smkn i.imd2010_5 i.ethnic if gender == `g', link(identity) vce(robust) allbaselevels family(gaussian)
		qui local np = e(N)
		preserve
		qui parmest, fast
		qui gen sex = `g'
		qui gen state = `s'
		qui gen npart = `np'
		qui append using `mlos'
		qui save `mlos', replace
		restore
		di "RelLOS -- Sex = `g' (N = `np') | State = `s' -- $S_TIME"
	}
}
qui use `mlos', clear
qui gen time = "10 years"
qui gen out = "States_Relative_LOS"
qui drop z p
qui drop if parm == "_cons"
qui sencode parm, gen(sparm)
qui sort sex state sparm
compress
save "$off/Results/States_Relative_los2.dta", replace


******************************************************************
**#Rate comparisons after Prediabetes or T2DM*********************
******************************************************************
frame reset
global off ".../Analysis"
use "$off/Datasets/multistate_structure.dta", replace 
stset _stop, enter(_start) failure(_status==1)
keep patid ageind gender _stop _start _trans _st-_t0
keep if _trans == 1 | _trans == 5
sort patid _trans

gen ageup = ageind + _start
sum ageup, d
egen float ageupg = cut(ageup), at(0 55 65 75 110) icodes label
tabstat ageup, statistics(N min max) by(ageupg)
tab ageupg, nolab
label define ageupg 0 `"<55"', modify
label define ageupg 1 `"≥55 to <65"', modify
label define ageupg 2 `"≥65 to <75"', modify
label define ageupg 3 `"≥75"', modify
tabstat ageup, statistics(N min max) by(ageupg)
groups gender ageupg, missing sepby(gender)
label variable ageind `"Age at prediabetes diagnosis (years)"'
label variable ageup  `"Updated age (years)"'
label variable ageupg `"Updated age group (years)"'
rename _trans state
tostring state, replace
replace state = "PreD" if state == "1"
replace state = "T2D" if state == "5"
sencode state, gsort(state) replace

cls
strate gender if state == 1, per(1000)
strate gender if state == 2, per(1000)
gen event = _d
drop _st-_t0
gen time = _stop - _start
stset time, failure(event==1)  
strate gender if state == 1, per(1000)
strate gender if state == 2, per(1000)
sum _t, d

forval j = 1/2 {
	qui stpm3 i.state##i.ageupg if gender == `j', scale(lncumhazard) df(4) tvc(i.state i.ageupg) dftvc(4) eform vce(cluster patid)
	qui predict pd0 pd1 pd2 pd3 d0 d1 d2 d3, hazard per(1000) timevar(0 10, step(0.1)) frame(irate`j', replace)	///
			at1(state 1 ageupg 0) 	///
			at2(state 1 ageupg 1) 	///
			at3(state 1 ageupg 2) 	///
			at4(state 1 ageupg 3) 	///
			at5(state 2 ageupg 0) 	///
			at6(state 2 ageupg 1) 	///
			at7(state 2 ageupg 2) 	///
			at8(state 2 ageupg 3)
}	
frame change irate1
gen sex = "Men"
frame change irate2
gen sex = "Women"
fframeappend, using(irate1)
tempfile irates
save `irates', replace

forval j = 1/2 {
	frame change default
	qui stpm3 i.state##i.ageupg if gender == `j', scale(lncumhazard) df(4) tvc(i.state i.ageupg) dftvc(4) eform vce(cluster patid)
	forval a = 0/3 {
		qui frame change default
		qui predict pd`a' d`a', hazard ci per(1000) timevar(0 10, step(0.1)) ///
			at1(state 1 ageupg `a') at2(state 2 ageupg `a') contrast(difference) contrastvars(ddp) frame(framed, replace)
		qui frame change framed	
		qui keep tt ddp*
		qui gen sex = `j'
		qui gen ageg = `a'
		qui tempfile s`j'_a`a'
		qui save `s`j'_a`a'', replace
	}
}
clear
forval j = 1/2 {
	forval a = 0/3 {
		append using `s`j'_a`a''
	}
}
reshape wide ddp ddp_lci ddp_uci, i(tt sex) j(ageg)
sort sex tt
tostring sex, replace
replace sex = "Men" if sex == "1"
replace sex = "Women" if sex == "2"
tempfile drates
save `drates', replace

use `irates', replace
merge 1:1 tt sex using `drates', update
order sex tt, first
gsort -sex tt
drop _merge
save "$off/Results/Rates_PreD_T2DM.dta", replace


******************************************************************
**#MULTIPLE IMPUTATION FOR GLM - PROBS AND LOS********************
******************************************************************

*--------------------------------------------*
***Probs - RELATIVE [DEATH COMBINED] [HPC 1]-*
*--------------------------------------------*
cls
clear all
global off ".../Analysis"
use "$off/Datasets/ps_probs2.dta", clear

rename ps_p?__25 ps_p?_0_25 
rename ps_p?__5  ps_p?_0_5 
rename ps_p?__75 ps_p?_0_75 

egen float nmiss = rowmiss(ps_p1_0-ps_p8_10)
tab nmiss, m
drop if nmiss != 0
drop nmiss
sencode ethnic, replace
gen ps_p99_10 = ps_p4_10 + ps_p8_10 + ps_p7_10 + ps_p3_10			/*State all death combined*/
keep patid ageg gender bmig imd2010_5 ethnic smk ps_p*_10
sdecode smk, replace
sencode smk, gsort(-smk) gen(smkn)
drop smk
order smkn, after(ethnic)
egen float cmiss = rowmiss(bmi smk ethnic imd2010_5)
tab cmiss gender, m
mdesc age bmi smk ethnic imd2010_5

stset, clear

preserve
clear
tempfile mprobs
save `mprobs', emptyok replace
restore

set seed 14795
local nimp = 10

forvalues g = 1/2 {
	foreach s in 1 2 5 6 99 {						/*Non-absorbing + death*/
		preserve
		qui keep if gender == `g'
		qui keep ps_p`s'_10 ageg bmig smkn imd2010_5 ethnic
		qui mi set mlong
		qui mi register imputed bmig smkn imd2010_5 ethnic		/*missing*/
		qui mi register regular ps_p`s'_10 ageg				/*outcome and non-missing*/
		qui mi impute chained (ologit) bmig imd2010_5   ///
				      (mlogit) smkn ethnic	///
				      = ps_p`s'_10 ageg, add(`nimp')
		qui mi estimate: glm ps_p`s'_10 i.ageg i.bmig i.smkn i.imd2010_5 i.ethnic, eform link(log) vce(robust) allbaselevels family(gaussian)
		qui parmest, fast
		qui gen sex = `g'
		qui gen state = `s'
		qui append using `mprobs'
		qui save `mprobs', replace
		di "RelProbs -- Sex = `g' | State = `s' -- $S_TIME"
		restore
	}
}
qui use `mprobs', clear
qui gen time = "10 years"
qui gen out = "States_Relative_Probs_MI"
qui drop t p dof
qui drop if parm == "_cons"
qui sencode parm, gen(sparm)
qui sort sex state sparm
compress
save "$off/Results/States_Relative_prob2_MI.dta", replace


*------------------------------------------*
***LOS - RELATIVE [DEATH COMBINED] [HPC 2]-*
*------------------------------------------*
cls
clear all
global off ".../Analysis"
use "$off/Datasets/ps_los2.dta", clear

rename ps_los?__25 ps_los?_0_25 
rename ps_los?__5  ps_los?_0_5 
rename ps_los?__75 ps_los?_0_75 

egen float nmiss = rowmiss(ps_los1_0-ps_los8_10)
tab nmiss, m
drop if nmiss != 0
drop nmiss
sencode ethnic, replace
gen ps_los99_10 = ps_los4_10 + ps_los8_10 + ps_los7_10 + ps_los3_10		/*State all death combined*/
keep patid ageg gender bmig imd2010_5 ethnic smk ps_los*_10
sdecode smk, replace
sencode smk, gsort(-smk) gen(smkn)
drop smk
order smkn, after(ethnic)
egen float cmiss = rowmiss(bmi smk ethnic imd2010_5)
tab cmiss gender, m
mdesc age bmi smk ethnic imd2010_5

stset, clear

preserve
clear
tempfile mlos
save `mlos', emptyok replace
restore

set seed 74579
local nimp = 10

forvalues g = 1/2 {
	foreach s in 1 2 5 6 99 {						/*Non-absorbing + death*/
		preserve
		qui keep if gender == `g'
		qui keep ps_los`s'_10 ageg bmig smkn imd2010_5 ethnic
		qui mi set mlong
		qui mi register imputed bmig smkn imd2010_5 ethnic		/*missing*/
		qui mi register regular ps_los`s'_10 ageg			/*outcome and non-missing*/
		qui mi impute chained (ologit) bmig imd2010_5   ///
				      (mlogit) smkn ethnic	///
				      = ps_los`s'_10 ageg, add(`nimp')
		qui mi estimate: glm ps_los`s'_10 i.ageg i.bmig i.smkn i.imd2010_5 i.ethnic, link(identity) vce(robust) allbaselevels family(gaussian)
		qui parmest, fast
		qui gen sex = `g'
		qui gen state = `s'
		qui append using `mlos'
		qui save `mlos', replace
		di "RelLOS -- Sex = `g' | State = `s' -- $S_TIME"
		restore
	}
}
qui use `mlos', clear
qui gen time = "10 years"
qui gen out = "States_Relative_LOS_MI"
qui drop t p dof
qui drop if parm == "_cons"
qui sencode parm, gen(sparm)
qui sort sex state sparm
compress
save "$off/Results/States_Relative_los2_MI.dta", replace



*####################################################################################################################################*
**************************************************************************************************************************************
**#OUTPUT -- GRAPHS -- TABLES -- HPC**************************************************************************************************
**************************************************************************************************************************************
*####################################################################################################################################*

cls
clear all

global off ".../Analysis"

*------------------------*
**#Table-1 BASELINE------*
*------------------------*
use "$off/Datasets/db2.dta", replace
decode gender, gen(gendern)
sencode gendern, gsort(-gendern) replace
sdecode smk, replace
sencode smk, gsort(-smk) gen(smkn)
drop smk
 
baselinetable   								/*
*/  agein(cts tab("p50 (p25, p75)")) 						/*
*/  ageg (cat)									/*
*/  bmi(cts tab("p50 (p25, p75)"))						/*
*/  bmig(cat) smkn(cat) imd2010_5(cat) ethnic(cat)				/*
*/  , by(gendern, total missing) reportmissing countformat(%10.0fc) notable 	/*
*/  exportexcel("$off/Results/Table1_baseline", replace)

sum surv, d
total surv
egen float cmiss = rowmiss(bmi smk ethnic imd2010_5)
tab cmiss gender, m


*------------------------*
**#Figure-1 TRANSITIONS--*
*------------------------*
use "$off/Datasets/multistate_structure.dta", replace 
matrix tmat = (.,1,.,2,3,.,.,.\.,.,4,.,.,.,.,.\.,.,.,.,.,.,.,.\.,.,.,.,.,.,.,.\.,.,.,.,.,5,.,6\.,.,.,.,.,.,7,.\.,.,.,.,.,.,.,.\.,.,.,.,.,.,.,.)
adopath ++ ".../Analysis_msm"

preserve
keep if gender == 1
msboxes, transmatrix(tmat) id(patid)                       					///
	xvalues(0.00 0.33 0.66 0.33 0.33 0.66 0.99 0.66)            				///
	yvalues(0.75 1.00 1.00 0.75 0.25 0.50 0.50 0.00)            				///
	statenames(pT2D_[S1] Kp_[S2] DKp_[S3] Dp_[S4] T2D_[S5] Kt_[S6] DKt_[S7] Dt_[S8]) 	///
	transnames(T1 T2 T3 T4 T5 T6 T7)							///
	boxheight(0.10) ysize(4) xsize(8) 	
graph save "Graph" "$off/Results/Figure1_TraM.gph", replace
restore

preserve
keep if gender == 2
msboxes, transmatrix(tmat) id(patid)                       					///
	xvalues(0.00 0.33 0.66 0.33 0.33 0.66 0.99 0.66)            				///
	yvalues(0.75 1.00 1.00 0.75 0.25 0.50 0.50 0.00)            				///
	statenames(pT2D_[S1] Kp_[S2] DKp_[S3] Dp_[S4] T2D_[S5] Kt_[S6] DKt_[S7] Dt_[S8]) 	///
	transnames(T1 T2 T3 T4 T5 T6 T7)							///
	boxheight(0.10) ysize(4) xsize(8) 	
graph save "Graph" "$off/Results/Figure1_TraW.gph", replace
restore

graph combine "$off/Results/Figure1_TraW.gph" 	///
	      "$off/Results/Figure1_TraM.gph"	///
	      , scale(0.7) xsize(9) ysize(5) imargin(small) nocopies
graph save "Graph" "$off/Results/Figure1_Transitions.gph", replace
graph close _all
qui erase "$off/Results/Figure1_TraM.gph"
qui erase "$off/Results/Figure1_TraW.gph"


*---------------------------------------*
**#Figure-2 PROBS CUMULATIVE------------*
*---------------------------------------*
use "$off/Results/States_Overall_probs2.dta", clear
drop se_prob filen metric
order sex ageg time state 
gsort sex ageg time state
bys sex ageg time (state): egen float totprob = total(prob)
sum totprob
drop totprob _ci*
label define ageg 0 "<55 years", modify
label define ageg 1 "≥55 to <65 years", modify
label define ageg 2 "≥65 to <75 years", modify
label define ageg 3 "≥75 years", modify

rename prob est_
reshape wide est_, i(sex ageg time) j(state)
gen c1 =      est_1
gen c2 = c1 + est_2
gen c3 = c2 + est_5
gen c4 = c3 + est_6
gen c5 = c4 + est_4
gen c6 = c5 + est_8
gen c7 = c6 + est_3
gen c8 = c7 + est_7

foreach var of varlist est_1-c8 {
	replace `var' = `var'*100
}

twoway (area  c1        time, sort fcolor(%50) fintensity(35) lwidth(none) yaxis(2)) 		///
	   (rarea c2 c1 time, sort fcolor(%50) fintensity(35) lwidth(none))  			///
	   (rarea c3 c2 time, sort fcolor(%50) fintensity(35) lwidth(none)) 			///
	   (rarea c4 c3 time, sort fcolor(%50) fintensity(35) lwidth(none)) 			///
	   (rarea c5 c4 time, sort fcolor(blue%80)   fintensity(50) lwidth(none)) 		///
	   (rarea c6 c5 time, sort fcolor(orange%80) fintensity(50) lwidth(none)) 		///
	   (rarea c7 c6 time, sort fcolor(purple%80) fintensity(50) lwidth(none)) 		///
	   (rarea c8 c7 time, sort fcolor(black%80)  fintensity(50) lwidth(none)) 		///
	   (line c1 time, sort lcolor(black) lwidth(vthin))					///
	   (line c2 time, sort lcolor(black) lwidth(vthin))					///
	   (line c3 time, sort lcolor(black) lwidth(vthin))					///
	   (line c4 time, sort lcolor(black) lwidth(vthin) lpattern(shortdash))			///
	   (line c5 time, sort lcolor(black) lwidth(vthin) lpattern(shortdash))			///
	   (line c6 time, sort lcolor(black) lwidth(vthin) lpattern(shortdash))			///
	   (line c7 time, sort lcolor(black) lwidth(vthin) lpattern(shortdash))			///
	   (line c8 time, sort lcolor(black) lwidth(vthin) lpattern(shortdash))			///
	   , ytitle(`"Probability (%)"') ytitle(, size(medsmall)) ylabel(0(10)100, format(%3.0f) axis(1)) ymtick(##2) 		///
	   ytitle(, size(medsmall)) ylabel(0(10)100, format(%3.0f) axis(2)) ymtick(##2, axis(2))	 						///
	   xtitle(`"Time since prediabetes diagnosis (years)"') xlabel(0(1)10) xmtick(##2) 									///
	   by(sex ageg, iyaxes note("") scale(0.6) rows(2)) subtitle(, nobox size(medlarge)) by(, legend(position(6)))  	///
	   legend(order(						  	  						///
	   1 "{bf:[S1]} pT2D" 						  							///		/*c1 - est_1*/
	   2 "{bf:[S2]} Cancer (pT2D)" 				  								///		/*c2 - est_2*/
	   3 "{bf:[S5]} T2D" 						  							///		/*c3 - est_5*/
	   4 "{bf:[S6]} Cancer (T2D)" 				  								///		/*c4 - est_6*/
	   5 "{bf:[S4]} Death (pT2D)"				  								///		/*c5 - est_4*/
	   6 "{bf:[S8]} Death (T2D)" 				  								///		/*c6 - est_8*/
	   7 "{bf:[S3]} Death (pT2D + Cancer)"		  									///		/*c7 - est_3*/
	   8 "{bf:[S7]} Death (T2D + Cancer)")		  									///		/*c8 - est_7*/
	   rows(2) colgap(minuscule) rowgap(minuscule) keygap(minuscule) symxsize(5) size(small) nobox region(fcolor(none))) 	///
	   name("Fig2", replace) xsize(6.5) ysize(3) by(, graphregion(margin(tiny))) nodraw
graph save "Fig2" "$off/Results/Figure2_CumPROBS.gph", replace
rename est_? state_?
export excel using "$off/Results/Aggregate_results.xls", sheet("Figure_2", modify) firstrow(variables) keepcellfmt
  

*------------------------------------*
**#Figure-3 LOS CUMULATIVE-----------*
*------------------------------------*
use "$off/Results/States_Overall_los2.dta", clear
drop se_los filen metric
order sex ageg time state 
gsort sex ageg time state
bys sex ageg time (state): egen float totlos = total(los)
sum totlos
drop totlos
label define ageg 0 "<55 years", modify
label define ageg 1 "≥55 to <65 years", modify
label define ageg 2 "≥65 to <75 years", modify
label define ageg 3 "≥75 years", modify

renames los _ci_lb _ci_ub \ est_ min95_ max95_ 
reshape wide est_ min95_ max95_, i(sex ageg time) j(state)
gen c1 =      est_1
gen c2 = c1 + est_2
gen c3 = c2 + est_5
gen c4 = c3 + est_6
gen c5 = c4 + est_4
gen c6 = c5 + est_8
gen c7 = c6 + est_3
gen c8 = c7 + est_7

twoway 	(area  c1    time, sort fcolor(%50) fintensity(35) lwidth(none) yaxis(2)) 	///
	(rarea c2 c1 time, sort fcolor(%50) fintensity(35) lwidth(none))  		///
	(rarea c3 c2 time, sort fcolor(%50) fintensity(35) lwidth(none)) 		///
	(rarea c4 c3 time, sort fcolor(%50) fintensity(35) lwidth(none)) 		///
	(rarea c5 c4 time, sort fcolor(blue%80)   fintensity(50) lwidth(none)) 		///
	(rarea c6 c5 time, sort fcolor(orange%80) fintensity(50) lwidth(none)) 		///
	(rarea c7 c6 time, sort fcolor(purple%80) fintensity(50) lwidth(none)) 		///
	(rarea c8 c7 time, sort fcolor(black%80)  fintensity(50) lwidth(none)) 		///
	(line c1 time, sort lcolor(black) lwidth(vthin))				///
	(line c2 time, sort lcolor(black) lwidth(vthin))				///
	(line c3 time, sort lcolor(black) lwidth(vthin))				///
	(line c4 time, sort lcolor(black) lwidth(vthin) lpattern(shortdash))		///
	(line c5 time, sort lcolor(black) lwidth(vthin) lpattern(shortdash))		///
	(line c6 time, sort lcolor(black) lwidth(vthin) lpattern(shortdash))		///
	(line c7 time, sort lcolor(black) lwidth(vthin) lpattern(shortdash))		///
	(line c8 time, sort lcolor(black) lwidth(vthin) lpattern(shortdash))		///
	, ytitle(`"Length of stay (years)"') ytitle(, size(medsmall)) ylabel(0(1)10, format(%3.0f) axis(1)) ymtick(##2) 	///
	ytitle(, size(medsmall)) ylabel(0(1)10, format(%3.0f) axis(2)) ymtick(##2, axis(2))	 				///
	xtitle(`"Time since prediabetes diagnosis (years)"') xlabel(0(1)10) xmtick(##2) 									///
	by(sex ageg, iyaxes note("") scale(0.6) rows(2)) subtitle(, nobox size(medlarge)) by(, legend(position(6))) 		///
	legend(order(						  	  							///
	1 "{bf:[S1]} pT2D" 						  							///		/*c1 - est_1*/
	2 "{bf:[S2]} Cancer (pT2D)" 				  								///		/*c2 - est_2*/
	3 "{bf:[S5]} T2D" 						  							///		/*c3 - est_5*/
	4 "{bf:[S6]} Cancer (T2D)" 				  								///		/*c4 - est_6*/
	5 "{bf:[S4]} Death (pT2D)"				  								///		/*c5 - est_4*/
	6 "{bf:[S8]} Death (T2D)" 				  								///		/*c6 - est_8*/
	7 "{bf:[S3]} Death (pT2D + Cancer)"		  									///		/*c7 - est_3*/
	8 "{bf:[S7]} Death (T2D + Cancer)")		  									///		/*c8 - est_7*/
	rows(2) colgap(minuscule) rowgap(minuscule) keygap(minuscule) symxsize(5) size(small) nobox region(fcolor(none))) 	///
	name("Fig3", replace) xsize(6.5) ysize(3) by(, graphregion(margin(tiny))) nodraw
graph save "Fig3" "$off/Results/Figure3_CumLOS.gph", replace
drop min95_* max95_*
rename est_? state_?   
export excel using "$off/Results/Aggregate_results.xls", sheet("Figure_3", modify) firstrow(variables) keepcellfmt


*--------------------------------------------*
**#Figure-4 (Individual Probs with 95%CI)----*
*--------------------------------------------*
use "$off/Results/States_Overall_probs2.dta", clear
drop se_prob filen metric
order sex ageg state time  
gsort sex ageg state time 
renames prob _ci_lb _ci_ub \ est min95 max95
keep if time == 10
drop time
foreach var of varlist est-max95 {
	replace `var' = `var'*100
}
renames est-max95, prefix(p_)
gsort sex state ageg
tempfile f4p
save `f4p', replace

*--------------------------------------------*
*#Figure-4 (Individual LOS with 95%CI)-------*
*--------------------------------------------*
use "$off/Results/States_Overall_los2.dta", clear
drop se_los filen metric
order sex ageg state time  
gsort sex ageg state time 
renames los _ci_lb _ci_ub \ est min95 max95
keep if time == 10
drop time
renames est-max95, prefix(los_)
gsort sex state ageg
tempfile f4l
save `f4l', replace

*--------------------------------------------------*
*#Figure-4 (Combined Probs + LOS with 95%CI)-------*
*--------------------------------------------------*
use `f4p', replace
merge 1:1 sex ageg state using `f4l', update
drop _merge
label define ageg 0 "<55 years", modify
label define ageg 1 "≥55 to <65 years", modify
label define ageg 2 "≥65 to <75 years", modify
label define ageg 3 "≥75 years", modify
tostring state, replace
replace state = "Prediabetes [S1]" 		if state == "1"
replace state = "Cancer (pT2D) [S2]" 		if state == "2"
replace state = "Death (pT2D + Cancer) [S3]" 	if state == "3"
replace state = "Death (pT2D) [S4]" 		if state == "4"
replace state = "Type 2 Diabetes [S5]" 		if state == "5"
replace state = "Cancer (T2D) [S6]" 		if state == "6"
replace state = "Death (T2D + Cancer) [S7]"	if state == "7"
replace state = "Death (T2D) [S8]" 		if state == "8"
sencode state, replace
gsort sex state ageg

twoway  (scatter los_est p_est 			 if ageg == 0, msize(tiny) msymbol(square) mcolor(blue))  				///
	(pcspike los_est p_min95 los_est p_max95 if ageg == 0, lwidth(vthin) lcolor(blue))  						///
	(pcspike los_min95 p_est los_max95 p_est if ageg == 0, lwidth(vthin) lcolor(blue))  						///
	(scatter los_est p_est 					if ageg == 1, msize(tiny) msymbol(square) mcolor(gold))  		///
	(pcspike los_est p_min95 los_est p_max95 if ageg == 1, lwidth(vthin) lcolor(gold))  						///
	(pcspike los_min95 p_est los_max95 p_est if ageg == 1, lwidth(vthin) lcolor(gold))  						///
	(scatter los_est p_est 					if ageg == 2, msize(tiny) msymbol(square) mcolor(red))   		///
	(pcspike los_est p_min95 los_est p_max95 if ageg == 2, lwidth(vthin) lcolor(red))   						///
	(pcspike los_min95 p_est los_max95 p_est if ageg == 2, lwidth(vthin) lcolor(red))   						///
	(scatter los_est p_est 					if ageg == 3, msize(tiny) msymbol(square) mcolor(green)) 		///
	(pcspike los_est p_min95 los_est p_max95 if ageg == 3, lwidth(vthin) lcolor(green)) 						///
	(pcspike los_min95 p_est los_max95 p_est if ageg == 3, lwidth(vthin) lcolor(green)) 						///
	, by(sex state, yrescale xrescale cols(4) note("")) subtitle(, size(vsmall) nobox)						///
	ytitle(`"Length of stay (years)"') ytitle(, size(vsmall)) ylabel(#5, labsize(vsmall) format(%3.2f) glpattern(dot)) ymtick(##2) 	///
	xtitle(`"Probability (%)"') xtitle(, size(vsmall)) xlabel(#5, labsize(vsmall) format(%03.1f) glpattern(dot)) xmtick(##2) 	///
	by(, legend(position(6))) 													///
	legend(order(															///
	1  "<55 years" 															///		
	4  "≥55 to <65 years" 														///		
	7  "≥65 to <75 years" 														///		
	10 "≥75 years") 														///		
	rows(1) colgap(minuscule) rowgap(minuscule) keygap(minuscule) size(vsmall) nobox region(fcolor(none)))				///
	name("Fig4", replace) xsize(6.5) ysize(5) by(, graphregion(margin(tiny))) nodraw
graph save "Fig4" "$off/Results/Figure4_ProbsLOS.gph", replace
export excel using "$off/Results/Aggregate_results.xls", sheet("Figure_4", modify) firstrow(variables) keepcellfmt


*-----------------------------------*
**#Table-S1 BASELINE BY AGE GROUP---*
*-----------------------------------*
use "$off/Datasets/db2.dta", replace
decode gender, gen(gendern)
gen agegn = ageg
tostring agegn, replace
replace gendern = "0Women" if gendern == "Women"
replace gendern = "1Men" if gendern == "Men"
gen gt = gendern + "_Ageg" + agegn
sdecode smk, replace
sencode smk, gsort(-smk) gen(smkn)
drop smk
 
baselinetable   							/*
*/  agein(cts tab("p50 (p25, p75)")) 					/*
*/  bmi(cts tab("p50 (p25, p75)"))					/*
*/  bmig(cat) smkn(cat) imd2010_5(cat) ethnic(cat)			/*
*/  , by(gt, total missing) reportmissing countformat(%10.0fc) notable 	/*
*/  exportexcel("$off/Results/TableS1_basAgeg", replace)


*------------------------*
**#Table-S2 RATES--------*
*------------------------*
use "$off/Datasets/multistate_structure.dta", replace 
stset _stop, enter(_start) failure(_status==1)

preserve
collapse (count) patid, by(gender _trans)
rename gender Sex
tempfile denoms
save `denoms', replace
restore

forvalues k = 1/7 {
	qui strate gender if _trans`k' == 1, per(1000) output("$off/Results/rates`k'_sex", replace)
}

clear
tempfile rates
save `rates', emptyok replace
forvalues k = 1/7 {
	qui use "$off/Results/rates`k'_sex", clear
	qui gen transition = `k'
	qui append using `rates'
	qui save `rates', replace
}
use `rates', replace
gsort -gender transition
gen 	tn = "Prediabetes [S1] -> Cancer [S2]" 		if transition == 1
replace tn = "Prediabetes [S1] -> Death [S4]"		if transition == 2
replace tn = "Prediabetes [S1] -> Type 2 Diabetes [S5]"	if transition == 3
replace tn = "Cancer [S2] -> Death [S3]"		if transition == 4
replace tn = "Type 2 Diabetes [S5] -> Cancer [S6]"	if transition == 5
replace tn = "Type 2 Diabetes [S5] -> Death [S8]"	if transition == 6
replace tn = "Cancer [S6] -> Death [S7]"		if transition == 7
renames gender _D _Y _Rate _Lower _Upper transition \ Sex Events PYRs1000 Rate1000 LB95 UB95 _trans

foreach var of varlist PYRs1000-UB95 {
	tostring `var', format(%5.1f) force replace 
} 
merge m:1 _trans Sex using `denoms', update
gen prop100 = Events*100/patid
tostring prop100, format(%5.1f) force replace 
tostring Events, replace
drop PYRs1000
gen Rate = Rate1000 + " (" + LB95 + ", " + UB95 + ")"
drop Rate1000-UB95
order Events, before(prop100)
renames patid _trans \ Npeople Transition
tostring Npeople, replace
tostring Transition, replace
replace Transition = "T" + Transition
order Sex Transition tn Npeople Events prop100 Rate
gsort -Sex Transition
drop _merge
split tn, parse(" -> ")
renames tn1 tn2 Events \ From To Ntrans
drop tn
order From To, after(Transition)
export excel using "$off/Results/TableS2_rates.xls", firstrow(variables) replace

forvalues k = 1/7 {
	qui erase "$off/Results/rates`k'_sex.dta"
}


*----------------------------------------------*
**#Table-S3 NUMBER OF SPECIFIC CANCER SITES----*
*----------------------------------------------*
use "$off/Datasets/multistate_structure.dta", replace 
stset _stop, enter(_start) failure(_status==1)
keep if _trans1 == 1 | _trans5 == 1
keep if _d == 1
keep patid gender cancer_site _trans1 _trans5
bys gender: tab _trans*
gen transition = 1 if _trans1 == 1
replace transition = 5 if _trans5 == 1
drop patid _trans*
contract gender transition cancer_site
gsort gender transition -_freq
bys gender transition: egen float totk = total(_freq)
gen percent = _freq*100/totk
bys gender transition: gen cump = sum(percent)
bys gender transition: gen seq = _n
bys gender transition: distinct cancer_site
bys gender: distinct cancer_site
keep if seq <=20
drop totk seq
renames _freq cump cancer_site \ frequency cum_perc icd10
tostring percent cum_perc, replace format(%5.1f) force
tostring frequency, replace format(%5.0fc) force
icd10 generate cancer_site = icd10, description version(2019)
order gender transition icd10 cancer_site, first
replace cancer = subinstr(cancer_site, "Malignant neoplasm: ", "", .)
replace cancer = subinstr(cancer_site, "Malignant neoplasm of ", "", .)
replace cancer = proper(cancer_site)
compress
export excel using "$off/Results/TableS3_cancersites.xls", firstrow(variables) keepcellfmt replace


*----------------------------------------------*
**#Table-S4 PREDICTED PROPORTIONS OF DEATHS----*
*----------------------------------------------*
use "$off/Results/States_Overall_probs2.dta", clear
drop se_prob filen metric
order sex ageg state time  
gsort sex ageg state time 
renames prob _ci_lb _ci_ub \ est min95 max95
keep if time == 10
drop time
keep if state == 3 | state == 4 | state == 7 | state == 8
tostring state, replace
replace state = "_S" + state
foreach var of varlist est-max95 {
	replace `var' = `var'*100
}
renames est-max95, prefix(p_)
gsort sex state ageg
reshape wide p_*, i(sex ageg) j(state) string
drop p_min* p_max*
rename p_est_* prob_*
egen float prob_death = rowtotal(prob_S3- prob_S8)
foreach j in 3 4 7 8 {
	gen rp_S`j' = prob_S`j'*100/prob_death
}
foreach var of varlist prob_death rp_* {
	tostring `var', replace force format(%03.1f)
}
drop prob_S*
order rp_S4 rp_S3 rp_S8 rp_S7, after(prob_death)
export excel using "$off/Results/TableS4_propdeath.xls", firstrow(variables) replace


*----------------------------------------------*
**#Table-S5 COMPARISON WITH VS W/OUT MISSING---*
*----------------------------------------------*
use "$off/Datasets/db2.dta", replace
decode gender, gen(gendern)
gen agegn = ageg
tostring agegn, replace
sdecode smk, replace
sencode smk, gsort(-smk) gen(smkn)
drop smk
replace gendern = "0Women" if gendern == "Women"
replace gendern = "1Men" if gendern == "Men"
egen float cmiss = rowmiss(bmi smk ethnic imd2010_5)
tab cmiss, m
replace cmiss = 1 if cmiss >=1
tostring cmiss, replace
gen gt = gendern + "_allmiss" + cmiss

baselinetable   						/*
*/  agein(cts tab("p50 (p25, p75)")) 				/*
*/  ageg (cat)							/*
*/  bmi(cts tab("p50 (p25, p75)"))				/*
*/  bmig(cat) smkn(cat) imd2010_5(cat) ethnic(cat)		/*
*/  , by(gt) countformat(%10.0fc)  				/*
*/  exportexcel("$off/Results/TableS5_missvars", replace)


*-----------------------------------------------------------*
**#REVIEWR ONLY -- COMPARISON WITH VS W/OUT SMOKING---------*
*-----------------------------------------------------------*
use "$off/Datasets/db2.dta", replace
decode gender, gen(gendern)
gen agegn = ageg
tostring agegn, replace
gen missmok = "1" if smk == .
replace missmok = "0" if smk != .
tab missmok smk, m
replace gendern = "0Women" if gendern == "Women"
replace gendern = "1Men" if gendern == "Men"
gen gt = gendern + "_missmok" + missmok

baselinetable   						/*
*/  agein(cts tab("p50 (p25, p75)")) 				/*
*/  ageg (cat)							/*
*/  bmi(cts tab("p50 (p25, p75)"))				/*
*/  bmig(cat) imd2010_5(cat) ethnic(cat)			/*
*/  , by(gt) countformat(%10.0fc)  				/*
*/  exportexcel("$off/Results/REVIEWER_missmok", replace)


**********************************************************
**#FOR DISCUSSION - EVER VISITING VS PROB OF STAY/REMAIN**
**********************************************************
*You start from the disease of interest, and sum all following states
*For example, for death these are only death states, as absorbing
use "$off/Results/States_Overall_probs2.dta", clear
drop se_prob filen metric
order sex ageg state time  
gsort sex ageg state time 
renames prob _ci_lb _ci_ub \ est min95 max95
keep if time == 10
drop time min95 max95
reshape wide est, i(sex ageg) j(state)
gen s5678_everT2D = est5 + est6 + est7 + est8
gen s2367_everCAN = est2 + est3 + est6 + est7
gen s3_s32 = est3 / (est3 + est2)
foreach var of varlist est1-s3_s32 {
	replace `var' = `var'*1000
	tostring `var', force format(%5.0f) replace
}
export excel using "$off/Results/Ever_Stay.xls", firstrow(variables) keepcellfmt replace


*------------------------------------*
**#Figure-S1 FLOWDIAGRAM-------------*
*------------------------------------*


*------------------------------------*
**#Figure-S2 TRANSITIONS FLOW SANKEY-*
*------------------------------------*
use "$off/Datasets/pseudocounts_structure.dta", replace 
net install sankey, from("https://raw.githubusercontent.com/asjadnaqvi/stata-sankey/main/installation/") replace
which sankey /*v 1.9 -- 25 Jun 2025 -- https://github.com/asjadnaqvi/stata-sankey*/

forvalues k = 1/2 {
	preserve
	keep if gender == `k'
	count if fromstate == 1
	local tot = r(N)
	keep patid fromstate tostate
	collapse (count) patid, by(fromstate tostate)
	drop if fromstate == tostate
	set obs `=_N+1'
	replace fromstate = 1 if fromstate == .
	replace tostate = 1 if tostate == .
	replace patid = `tot' if patid == .
	gen gender = `k'
	tempfile res`k'
	save `res`k'', replace 
	restore
}

clear
forvalues k = 1/2 {
	append using `res`k''
}
gen layer     = 1 if fromstate == 1 & tostate == 1				/*layers defined for non-absorbing states*/
replace layer = 2 if fromstate == 1 & tostate != 1
replace layer = 3 if fromstate == 2 | fromstate == 5
replace layer = 4 if fromstate == 6
sort gender layer fromstate
gen genderl = "Men" if gender == 1
replace genderl = "Women" if gender == 2
label drop statelbl
foreach var of varlist fromstate tostate {
	tostring `var', replace force
	replace `var' = "S" + `var'
}	

foreach nm in Men Women {
	sankey patid if genderl == "`nm'", from(fromstate) to(tostate) by(layer) smooth(5) recenter(top) 	///
	gap(20) boxwid(5) colorby(level) labangle(90) labg(-10) alpha(50) lwidth(vthin) lcolor(black%50) 	///
	stock format(%12.2gc) name("Sankey_`nm'", replace) title("`nm'", size(medsmall)) novalues
}
graph combine Sankey_Men Sankey_Women, cols(1) nocopies xsize(5) ysize(7) name("Sankey", replace) nodraw
graph save "Sankey" "$off/Results/FigureS2_Sankey.gph", replace
graph close _all


*-----------------------------------*
**#Figure-S3 RATES PRED VS T2DM-----*
*-----------------------------------*
use "$off/Results/Rates_PreD_T2DM.dta", clear
sencode sex, replace

twoway (line pd0 pd1 pd2 pd3 tt,     sort lcolor(blue gold red green)) 										///
	   (line d0   d1  d2  d3 tt, sort lcolor(blue gold red green) lpattern(shortdash ..))							///
	   , by(sex, note("") scale(0.8)) subtitle(, nobox)											///
	   yscale(log range(3 60)) ytitle("Cancer incidence rate" "(per 1000 person-years)") ylabel(3 5 10 20 40 60) 				///
	   xtitle("Time since diagnosis of prediabetes or type 2 diabetes (years)") xlabel(0(1)10) xsize(6) ysize(3.5) by(, legend(pos(6)))	///
	   legend(order(															///
	   1 "pT2DM <55 years" 															///		
	   2 "pT2DM ≥55 to <65 years" 														///		
	   3 "pT2DM ≥65 to <75 years" 														///		
	   4 "pT2DM ≥75 years" 															///
	   5 "  T2DM <55 years" 														///		
	   6 "  T2DM ≥55 to <65 years" 														///		
	   7 "  T2DM ≥65 to <75 years" 														///		
	   8 "  T2DM ≥75 years")														///	   
	   rows(2) colgap(small) rowgap(minuscule) keygap(minuscule) size(small) nobox region(fcolor(none)))					///
	   name("FigS3r", replace) nodraw
	   
twoway (line ddp0 ddp1 ddp2 ddp3 tt,     sort lcolor(blue gold red green)) 						///
	   (rarea ddp_lci0 ddp_uci0 tt,  sort fcolor(blue%20)  lwidth(none))						///
	   (rarea ddp_lci1 ddp_uci1 tt,  sort fcolor(gold%20)  lwidth(none))						///
	   (rarea ddp_lci2 ddp_uci2 tt,  sort fcolor(red%20)   lwidth(none))						///
	   (rarea ddp_lci3 ddp_uci3 tt,  sort fcolor(green%20) lwidth(none))						///	   
	   , by(sex, note("") scale(0.8)) subtitle(, nobox)	by(, legend(off))					///
	   ytitle("Cancer incidence rate difference" "(per 1000 person-years)") ylabel(-4	0 4 8 12 16 20 24 28)	///
	   xtitle("Time since diagnosis of prediabetes or type 2 diabetes (years)") xlabel(0(1)10) xsize(6) ysize(3.5) 	///
	   name("FigS3d", replace) nodraw
	   
graph combine FigS3r FigS3d, xcommon cols(1) xsize(7) ysize(6) scale(0.9) name("FigS3", replace) nocopies
graph save "FigS3" "$off/Results/FigureS3_RatesPreDT2DM.gph", replace
export excel using "$off/Results/Aggregate_results.xls", sheet("FigS3", modify) firstrow(variables) keepcellfmt
graph close _all


*-----------------------------------*
**#Figure-S4 STRATIFIED EVENTS------*
*-----------------------------------*
clear
cd "$off/Results/"									
foreach m in ageg bmig smkn imd2010_5 ethnic {
	qui save "Events_`m'", emptyok replace					/*This is necessary for saving, replace below*/
}

cd "$off/Results/"	
use "$off/Datasets/multistate_structure.dta", replace
egen float cmiss = rowmiss(bmi smk ethnic imd2010_5)
drop if cmiss>0
distinct patid if gender == 1
distinct patid if gender == 2
sdecode smk, replace
sencode smk, gsort(-smk) gen(smkn)
drop smk

qui stset _stop, enter(_start) failure(_status==1)
foreach m in ageg bmig smkn imd2010_5 ethnic {
	qui groups gender _trans `m' if _status == 1, sepby(_trans) saving(Events_`m', replace)
}

clear
tempfile sevents
save `sevents', emptyok replace
foreach m in ageg bmig smkn imd2010_5 ethnic {
	qui cap use "$off/Results/Events_`m'", clear
	qui cap sdecode `m', replace
	qui cap tostring `m', replace
	qui cap rename `m' level
	qui gen groupn = "`m'"
	qui append using `sevents'
	qui save `sevents', replace
}
use `sevents', clear
gen variablen 	  = 1 if groupn == "ageg"
replace variablen = 2 if groupn == "bmig"
replace variablen = 3 if groupn == "smkn"
replace variablen = 4 if groupn == "imd2010_5"
replace variablen = 5 if groupn == "ethnic"
gsort -gender variablen level _trans
drop _percent variablen
renames gender _trans _freq \ sex trans events_tra
replace groupn = "Age (years)" 				if groupn == "ageg"
replace groupn = "Body mass index (kg/m2)"		if groupn == "bmig"
replace groupn = "Smoking status"			if groupn == "smkn"
replace groupn = "Index of multiple deprivation 2010" 	if groupn == "imd2010_5"
replace groupn = "Ethnicity" 				if groupn == "ethnic"
replace level  = "1st fifth (least deprived)" 		if level == "1"
replace level  = "2nd" 					if level == "2"
replace level  = "3rd" 					if level == "3"
replace level  = "4th" 					if level == "4"
replace level  = "5th fifth (most deprived)"  		if level == "5"
sencode level, replace
sencode groupn, replace
order sex groupn level trans events_tra
reshape wide events_, i(sex groupn level) j(trans)
egen float tot_row = rowtotal(events*)
sdecode sex, replace
sencode sex, gsort(-sex) replace
gsort sex groupn level

graph hbar (asis) events_tra*, 											///
	  over(level, label(labsize(vsmall))) over(sex, label(labsize(small))) stack scale(0.8)			///
	  bar(1, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(2, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(3, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(4, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(5, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(6, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(7, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  ylabel(0(10000)70000, labsize(vsmall) glpattern(solid)) ymtick(##5) ytitle("No. of individuals")	///
	  legend(position(2) ring(0)) legend(order(1 "T1" 2 "T2" 3 "T3" 4 "T4" 5 "T5" 6 "T6" 7 "T7")		///
	  cols(1) colgap(minuscule) rowgap(minuscule) keygap(minuscule) symxsize(2) symysize(2) size(vsmall) nobox region(fcolor(none))) ///
	  name("FigS4a", replace) xsize(5.5) ysize(6) nodraw

graph hbar (asis) events_tra*, percentages 									///
	  over(level, label(labsize(vsmall))) over(sex, label(labsize(small))) stack scale(0.8)			///
	  bar(1, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(2, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(3, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(4, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(5, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(6, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  bar(7, fintensity(35) lcolor(black) lwidth(vvthin)) 							///
	  ylabel(0(10)100, labsize(vsmall) glpattern(solid)) ymtick(##5) ytitle("Percent") legend(off) 		///
 	  name("FigS4b", replace) xsize(5.5) ysize(6) nodraw
	 
graph combine FigS4a FigS4b, cols(2) nocopies name("FigS4", replace) nodraw	  
graph save "FigS4" "$off/Results/FigureS4_StratifiedEvents.gph", replace
export excel using "$off/Results/Aggregate_results.xls", sheet("FigS4", modify) firstrow(variables) keepcellfmt
graph close _all

foreach m in ageg bmig smkn imd2010_5 ethnic {
	qui cap erase "$off/Results/Events_`m'.dta"
}


*-----------------------------------------------*
**#Figure-S5/6 MI - Risk Ratio Forest Probs-----*
*-----------------------------------------------*
use "$off/Results/States_Relative_prob2.dta", clear
gen mod = "CC"
append using "$off/Results/States_Relative_prob2_MI.dta"
replace mod = "MI" if mod == ""
replace min95 = 0 if min95 == .
replace max95 = 0 if max95 == .
drop eq npart time out sparm
order sex state, first
split parm, p(".")
drop parm1
rename parm2 plotid
tostring sex, replace
replace plotid = "Age (years)"  if plotid == "ageg"
replace plotid = "BMI (kg/m2)"  if plotid == "bmig"
replace plotid = "Smoking"  	if plotid == "smkn"
replace plotid = "IMD (fifths)" if plotid == "imd2010_5"
replace plotid = "Ethnicity" 	if plotid == "ethnic"
gen plotidn = plotid
replace plotidn = ""			if estimate != 0 | mod == "MI"
replace sex = "Men" 			if sex == "1"
replace sex = "Women" 			if sex == "2"
replace parm = "<55"			if parm == "0b.ageg"
replace parm = "≥55 to <65"		if parm == "1.ageg"
replace parm = "≥65 to <75"		if parm == "2.ageg"
replace parm = "≥75"			if parm == "3.ageg"
replace parm = "<25"			if parm == "0b.bmig"
replace parm = "≥25 to <30"		if parm == "1.bmig"
replace parm = "≥30 to <35"		if parm == "2.bmig"
replace parm = "≥35 to <40"		if parm == "3.bmig"
replace parm = "≥40"			if parm == "4.bmig"
replace parm = "Non"			if parm == "1b.smkn"
replace parm = "Former"			if parm == "2.smkn"
replace parm = "Current"		if parm == "3.smkn"
replace parm = "1st - least"		if parm == "1b.imd2010_5"
replace parm = "2nd"			if parm == "2.imd2010_5"
replace parm = "3rd"			if parm == "3.imd2010_5"
replace parm = "4th"			if parm == "4.imd2010_5"
replace parm = "5th - most"		if parm == "5.imd2010_5"
replace parm = "White"			if parm == "1b.ethnic"
replace parm = "Non-White"		if parm == "2.ethnic"
gen parmcol = parm
sencode parm, replace
gsort -sex state parm
tostring state, replace
replace state = "Prediabetes (pT2D) [S1]" if state == "1"
replace state = "Cancer (pT2D) [S2]" 	  if state == "2"
replace state = "Type 2 Diabetes [S5]" 	  if state == "5"
replace state = "Cancer (T2D) [S6]" 	  if state == "6"
replace state = "Death [S3+S4+S7+S8]" 	  if state == "99"
sencode state, replace
sort sex state parm mod
drop if estimate == 0 & mod == "MI"
replace plotid = plotid + "_" + mod
replace parmcol = "" if mod == "MI"

local xl = ln(1)
foreach n in Women Men {
	forvalues s = 1/5 {
		qui local mylabel: label (state) `s'
		qui forestplot estimate min95 max95 if sex == "`n'" & state == `s'				///
			, eform effect("") lcols(plotidn parmcol)                        	 		///
			nonull nonames noov nosu nowt dp(2) classic boxscale(80) astext(50) textsize(200) 	///
			xlabel(, labsize(9pt) format(%3.1f) nogrid)  						///
			xline(`xl', lpattern(solid) lwidth(vthin)) 						///
			spacing(4.5) yline(3.5 12.5 17.5 26.5, lwidth(vthin) lpattern(vshortdash)) 		///
			xtitle("Risk Ratio", size(9pt))	leftjustify ciopts(lwidth(vthin)) plotid(plotid)    	///
			box1opts(mcolor(red)      ms(square))        ci1opts(lcolor(red))    			///
			box2opts(mcolor(red)      ms(circle_hollow)) ci2opts(lcolor(red))   			///
			box3opts(mcolor(blue)     ms(square))        ci3opts(lcolor(blue))    			///
			box4opts(mcolor(blue)     ms(circle_hollow)) ci4opts(lcolor(blue))   			///
			box5opts(mcolor(black)    ms(square))        ci5opts(lcolor(black))    			///
			box6opts(mcolor(black)    ms(circle_hollow)) ci6opts(lcolor(black))   			///
			box7opts(mcolor(orange)   ms(square))        ci7opts(lcolor(orange))    		///
			box8opts(mcolor(orange)   ms(circle_hollow)) ci8opts(lcolor(orange))   			///
			box9opts(mcolor(dkgreen)  ms(square))        ci9opts(lcolor(dkgreen))    		///
			box10opts(mcolor(dkgreen) ms(circle_hollow)) ci10opts(lcolor(dkgreen))   		///
			title("`n', `mylabel'", size(medsmall)) name("P`n'_`s'", replace) xsize(6) ysize(9) scale(0.75) nodraw
		}
}
graph combine PWomen_1 PWomen_2 PWomen_3 PWomen_4 PWomen_5, nocopies ycommon rows(1) xsize(11) ysize(5) scale(1.3) name("FigS5", replace)
graph save "FigS5" "$off/Results/FigureS5_RelPROBS_W.gph", replace
graph combine PMen_1 PMen_2 PMen_3 PMen_4 PMen_5, nocopies ycommon rows(1) xsize(11) ysize(5) scale(1.3) name("FigS6", replace)
graph save "FigS6" "$off/Results/FigureS6_RelPROBS_M.gph", replace
drop _EFFECT stderr parmcol plotidn
order sex state plotid parm mod
replace plotid = substr(plotid, 1, length(plotid) - 3) 
gsort -sex state parm mod
replace mod = "Complete-case" if mod == "CC"
replace mod = "Multiple imputation (10)" if mod == "MI"
replace mod = "(Reference)" if estimate == 0
foreach var of varlist estimate-max95 {
	replace `var' = exp(`var')
}
compress
export excel using "$off/Results/Aggregate_results.xls" if sex == "Women", sheet("FigS5", modify) firstrow(variables) keepcellfmt
export excel using "$off/Results/Aggregate_results.xls" if sex == "Men",   sheet("FigS6", modify) firstrow(variables) keepcellfmt


*------------------------------------------*
**#Figure-S7/8 MI Forest LOS difference----*
*------------------------------------------*
use "$off/Results/States_Relative_los2.dta", clear
gen mod = "CC"
append using "$off/Results/States_Relative_los2_MI.dta"
replace mod = "MI" if mod == ""
drop eq npart time out sparm
order sex state, first
split parm, p(".")
drop parm1
rename parm2 plotid
tostring sex, replace
replace plotid = "Age (years)"  if plotid == "ageg"
replace plotid = "BMI (kg/m2)"  if plotid == "bmig"
replace plotid = "Smoking"  	if plotid == "smkn"
replace plotid = "IMD (fifths)" if plotid == "imd2010_5"
replace plotid = "Ethnicity" 	if plotid == "ethnic"
gen plotidn = plotid
replace plotidn = ""			if estimate != 0
replace sex = "Men" 			if sex == "1"
replace sex = "Women" 			if sex == "2"
replace parm = "<55"			if parm == "0b.ageg"
replace parm = "≥55 to <65"		if parm == "1.ageg"
replace parm = "≥65 to <75"		if parm == "2.ageg"
replace parm = "≥75"			if parm == "3.ageg"
replace parm = "<25"			if parm == "0b.bmig"
replace parm = "≥25 to <30"		if parm == "1.bmig"
replace parm = "≥30 to <35"		if parm == "2.bmig"
replace parm = "≥35 to <40"		if parm == "3.bmig"
replace parm = "≥40"			if parm == "4.bmig"
replace parm = "Non"		    	if parm == "1b.smkn"
replace parm = "Former"			if parm == "2.smkn"
replace parm = "Current"		if parm == "3.smkn"
replace parm = "1st - least"		if parm == "1b.imd2010_5"
replace parm = "2nd"			if parm == "2.imd2010_5"
replace parm = "3rd"			if parm == "3.imd2010_5"
replace parm = "4th"			if parm == "4.imd2010_5"
replace parm = "5th - most"		if parm == "5.imd2010_5"
replace parm = "White"			if parm == "1b.ethnic"
replace parm = "Non-White"		if parm == "2.ethnic"
gen parmcol = parm
sencode parm, replace
gsort -sex state parm
tostring state, replace
replace state = "Prediabetes (pT2D) [S1]" if state == "1"
replace state = "Cancer (pT2D) [S2]" 	  if state == "2"
replace state = "Type 2 Diabetes [S5]" 	  if state == "5"
replace state = "Cancer (T2D) [S6]" 	  if state == "6"
replace state = "Death [S3+S4+S7+S8]" 	  if state == "99"
sencode state, replace
sort sex state parm mod
drop if estimate == 0 & mod == "MI"
replace plotid = plotid + "_" + mod
replace parmcol = "" if mod == "MI"

foreach n in Women Men {
	forvalues s = 1/5 {
		qui local mylabel: label (state) `s'
		qui forestplot estimate min95 max95 if sex == "`n'" & state == `s'				///
			, effect("") lcols(plotidn parmcol)                        	 			///
			nonull nonames noov nosu nowt dp(2) classic boxscale(80) astext(50) textsize(200) 	///
			xlabel(, labsize(9pt) format(%3.2f) nogrid)  						///
			xline(0, lpattern(solid) lwidth(vthin)) 						///
			spacing(4.5) yline(3.5 12.5 17.5 26.5, lwidth(vthin) lpattern(vshortdash)) 		///
			xtitle("Length of stay difference (years)", size(9pt))					///
			leftjustify ciopts(lwidth(vthin)) plotid(plotid)    					///
			box1opts(mcolor(red)      ms(square))        ci1opts(lcolor(red))    			///
			box2opts(mcolor(red)      ms(circle_hollow)) ci2opts(lcolor(red))   			///
			box3opts(mcolor(blue)     ms(square))        ci3opts(lcolor(blue))    			///
			box4opts(mcolor(blue)     ms(circle_hollow)) ci4opts(lcolor(blue))   			///
			box5opts(mcolor(black)    ms(square))        ci5opts(lcolor(black))    			///
			box6opts(mcolor(black)    ms(circle_hollow)) ci6opts(lcolor(black))   			///
			box7opts(mcolor(orange)   ms(square))        ci7opts(lcolor(orange))    		///
			box8opts(mcolor(orange)   ms(circle_hollow)) ci8opts(lcolor(orange))   			///
			box9opts(mcolor(dkgreen)  ms(square))        ci9opts(lcolor(dkgreen))    		///
			box10opts(mcolor(dkgreen) ms(circle_hollow)) ci10opts(lcolor(dkgreen))   		///
			title("`n', `mylabel'", size(medsmall)) name("L`n'_`s'", replace) xsize(6) ysize(9) scale(0.75) nodraw
		}
}
graph combine LWomen_1 LWomen_2 LWomen_3 LWomen_4 LWomen_5, nocopies ycommon rows(1) xsize(11) ysize(5) scale(1.3) name("FigS7", replace)
graph save "FigS7" "$off/Results/FigureS7_RelLOS_W.gph", replace
graph combine LMen_1 LMen_2 LMen_3 LMen_4 LMen_5, nocopies ycommon rows(1) xsize(11) ysize(5) scale(1.3) name("FigS8", replace)
graph save "FigS8" "$off/Results/FigureS8_RelLOS_M.gph", replace
drop _EFFECT stderr parmcol plotidn
order sex state plotid parm mod
replace plotid = substr(plotid, 1, length(plotid) - 3) 
gsort -sex state parm mod
replace mod = "Complete-case" if mod == "CC"
replace mod = "Multiple imputation (10)" if mod == "MI"
replace mod = "(Reference)" if estimate == 0
compress
export excel using "$off/Results/Aggregate_results.xls" if sex == "Women", sheet("FigS7", modify) firstrow(variables) keepcellfmt
export excel using "$off/Results/Aggregate_results.xls" if sex == "Men",   sheet("FigS8", modify) firstrow(variables) keepcellfmt
graph close _all
clear all

