clear all 

* specify a working directory (cd) 
* path has to be specified for this file to work (delete //)

local path "/Users/robert/Desktop/JCA"

mkdir cross_site_imputation

cd "`path'/cross_site_imputation"

*** install mi impute from ***
ssc install mi_impute_from, replace

********************************************************************************

*** generate data ***

********************************************************************************

cap program drop singlestudy
program define singlestudy, rclass

  syntax /// 
  	[, nobs(real 1000) ///
  	study(real 1) ///
  	filename(string)]
  
  drop _all 
  local nobs = runiformint(500,10000)
  qui set obs `nobs'
  
  local pc =  runiform(0.1, 0.4)
  gen c = rbinomial(1, `pc')
  
  	
  local a0 = runiform(0.05, 0.1)
  local a1 = rnormal(ln(4), 0.1)	
  gen x = rbinomial(1, invlogit(`a0'+`a1'*c))

  local b0 = runiform(logit(0.05), logit(0.3))
  local b1 = rnormal(ln(1.5), 0.1)
  local b2 = rnormal(ln(4), 0.1)
  gen y = rbinomial(1, invlogit(`b0'+`b1'*x+`b2'*c))
  
  gen study = `study'
  if "`filename'" != "" qui save "`filename'", replace 
end

*** generate studies ****

forv i = 1/5 {
	singlestudy, study(`i') filename(study_`i')
}

********************************************************************************

*** control approach ****

********************************************************************************

capture frame drop metadata
frame create metadata 
qui frame metadata:{ 
	set obs 5 // number of studies
	gen effect = .
	gen se = .
	gen study = .
	gen size = .
}  


forv i = 1/5{
	use study_`i', replace
	logit y x c
	local obs = _N
	frame metadata:{
		replace effect = _b[x] if _n == `i'
		replace se = _se[x] if _n == `i'
		replace study = `i' if _n == `i'
		replace size = `obs'  if _n == `i'
	}
}


frame metadata: list 
frame metadata: meta set effect se , studysize(size)
frame metadata: meta summarize,  eform fixed

*** generate missing data ***

forv i = 4/5{
	qui use study_`i', clear 
	qui count 
	di in red "Study `i': N = " r(N)
	replace c = . 
	save study_`i', replace
}


*** fit imputation model in all studies with available data ***

forv i = 1/3 {
	qui use study_`i', replace
	qui logit c y x
	mat ib = e(b)
	mat iV = e(V)
	svmat ib
	qui export delimited ib* using b_study`i'.txt in 1 , replace 
	svmat iV 
	qui export delimited iV* using v_study`i'.txt if iV1 != . , replace 
}


local b_file "b_study1 b_study2 b_study3"
local v_file "v_study1 v_study2 v_study3"

********************************************************************************

*** unadjusted pooled estimate ***

********************************************************************************

forv i = 1/5{
	use study_`i', replace
	logit y x
	local obs = _N
	frame metadata:{
		replace effect = _b[x] if _n == `i'
		replace se = _se[x] if _n == `i'
		replace study = `i' if _n == `i'
		replace size = `obs'  if _n == `i'
	}
}

frame metadata: list 
frame metadata: meta set effect se , studysize(size)
frame metadata: meta summarize,  eform fixed

********************************************************************************

*** adjusted complete case analysis ***

********************************************************************************

capture frame drop metadata
frame create metadata 
qui frame metadata:{ 
	set obs 3
	gen effect = .
	gen se = .
	gen study = .
	gen size = .
}  

forv i = 1/3{
	use study_`i', replace
	logit y x c
	local obs = _N
	frame metadata:{
		replace effect = _b[x] if _n == `i'
		replace se = _se[x] if _n == `i'
		replace study = `i' if _n == `i'
		replace size = `obs'  if _n == `i'
	}

}

frame metadata: meta set effect se , studysize(size)
frame metadata: meta summarize,  eform fixed

********************************************************************************

*** adjusted for available data estimate ***

********************************************************************************

capture frame drop metadata
frame create metadata 
qui frame metadata:{ 
	set obs 5
	gen effect = .
	gen se = .
	gen study = .
	gen size = .
}  

forv i = 1/3{
	use study_`i', replace
	logit y x c
	local obs = _N
	frame metadata: replace effect = _b[x] if _n == `i'
	frame metadata: replace se = _se[x] if _n == `i'
	frame metadata: replace study = `i' if _n == `i'
	frame metadata: replace size = `obs'  if _n == `i'
}


forv i = 4/5{
	use study_`i', replace
	logit y x
	local obs = _N
		frame metadata:{
		replace effect = _b[x] if _n == `i'
		replace se = _se[x] if _n == `i'
		replace study = `i' if _n == `i'
		replace size = `obs' if _n == `i'
	}
}


frame metadata: meta set effect se , studysize(size)
frame metadata: meta summarize, eform fixed

********************************************************************************

*** single study imputation adjusted estimate ***

********************************************************************************

local pick = runiformint(1,3) 
use study_`pick', clear
	logit c y x
 
	// export matrices 
	mat ib = e(b)
	mat iV = e(V)
	svmat ib
	qui export delimited ib* using b_study`pick'.txt in 1 , replace 
	svmat iV 
	qui export delimited iV* using v_study`pick'.txt if iV1 != . , replace 


forv i = 4/5 {

	use study_`i', clear 
	mi set wide
	mi register imputed c
	
	mi_impute_from_get, b(b_study`pick') v(v_study`pick') colnames(y x _cons) imodel(logit) 
	mat ib = r(get_ib)
	mat iV = r(get_iV)
	
	mi impute from c , add(10) b(ib) v(iV) imodel(logit)

	mi estimate, post noi: logit y x c

	local obs = _N
	frame metadata:{
		replace effect = _b[x] if _n == `i'
		replace se = _se[x] if _n == `i'
		replace study = `i' if _n == `i'
		replace size = `obs'  if _n == `i'
	}
}


forv i = 1/3{
	use study_`i', replace
	logit y x c
	local obs = _N
	frame metadata{
		replace effect = _b[x] if _n == `i'
		replace se = _se[x] if _n == `i'
		replace study = `i' if _n == `i'
		replace size = `obs'  if _n == `i'
	}
}

frame metadata: meta set effect se , studysize(size)
frame metadata: meta summarize, eform fixed

********************************************************************************

*** multiple study imputation adjusted pooled estimate ***

********************************************************************************

forv i = 1/3 {
	qui use study_`i', replace
	qui logit c y x
	mat ib = e(b)
	mat iV = e(V)
	svmat ib
	qui export delimited ib* using b_study`i'.txt in 1 , replace 
	svmat iV 
	qui export delimited iV* using v_study`i'.txt if iV1 != . , replace 
}

local b_file "b_study1 b_study2 b_study3"
local v_file "v_study1 v_study2 v_study3"

// use the list of files in mi_impute_from_get

forv k = 4/5{
	
	// impute in study 4 & 5 
	use study_`k', clear 
		
	mi set wide
	mi register imputed c
		
	mat drop _all
	mi_impute_from_get, b(`b_file') v(`v_file') colnames(y x _cons) imodel(logit) // weighted average is automatically taken
	mat ib = r(get_ib)
	mat iV = r(get_iV)
	
	mi impute from c , b(ib) v(iV) add(10) imodel(logit)
	mi estimate, post noi: logit y x c

	local obs = _N
	frame metadata:{
		replace effect = _b[x] if _n == `k'
		replace se = _se[x] if _n == `k'
		replace study = `k' if _n == `k'
		replace size = `obs' if _n == `k'
	}
}


forv i = 1/3{
	use study_`i', replace
	logit y x c
	local obs = _N
	frame metadata{
		replace effect = _b[x] if _n == `i'
		replace se = _se[x] if _n == `i'
		replace study = `i' if _n == `i'
		replace size = `obs'  if _n == `i'
	}
}

frame metadata: meta set effect se , studysize(size)
frame metadata: meta summarize,  eform fixed

********************************************************************************

*** Empirical heterogenity ***

********************************************************************************

*** fit imputation model in all studies with available data *** (already been done previously)
forv i = 1/3 {
	qui use study_`i', replace
	qui logit c y x
	mat ib = e(b)
	mat iV = e(V)
	svmat ib
	qui export delimited ib* using b_study`i'.txt in 1 , replace 
	svmat iV 
	qui export delimited iV* using v_study`i'.txt if iV1 != . , replace 
}


local b_file "b_study1 b_study2 b_study3"
local v_file "v_study1 v_study2 v_study3"


** import 3 files with beta-coefficients from prediction model

forv i = 1/3 {
	tempname get_b_`i'
	qui import delimited using "b_study`i'.txt", clear 
	qui mkmat * , matrix(`get_b_`i'')	
}

** import the 3 files with the variance/covariance matrix from the prediction model

forv i = 1/3 {
	tempname get_v_`i'
	qui import delimited using "v_study`i'.txt", clear 
	qui mkmat *, matrix(`get_v_`i'')
}


** All coefficients from the 3 prediction models are now stored in form of matrices
** put all matrices into a frame (create a dataframe)

cap frame drop random_impmodel  
frame create random_impmodel study y1 y2 y3 v11 v12 v13 v22 v23 v33
*frame random_impmodel: desc 

forv t = 1/3 {
 
	local post_est_b_`t' ""
	local post_est_V_`t' ""
	
	forv i = 1/3 {
		local post_est_b_`t' "`post_est_b_`t'' (`=`get_b_`t''[1,`i']')"
	}
	
	forv i = 1/3 {
		forv j = `i'/3 {   // triangular and diagonal elements of the matrix for the V matrix
            local post_est_V_`t' "`post_est_V_`t'' (`=`get_v_`t''[`i',`j']')"
        }
    }
	
	qui frame post random_impmodel (`t') `post_est_b_`t'' `post_est_V_`t''
	
}

*frame random_impmodel: list

** fit meta regrssion model with REML and identity var/cov matrix
frame random_impmodel: /// 
	meta mvregress y* , wcovvariables(v*) random(reml, covariance(identity))
	

** clean up to be used with mi impute from

mat ib = e(b)[1, 1..3]
mat colnames ib = y x _cons

mat iV = e(V)[1..3, 1..3]
mat colnames iV = y x _cons
mat rownames iV = y x _cons

mat li ib 
mat li iV

** use the matrices in the studies with incomplete data

forv k = 4/5{
	
	// impute in study 4 & 5 
	use study_`k', clear 
		
	mi set wide
	mi register imputed c
		
	mi impute from c , b(ib) v(iV) add(10) imodel(logit)
	mi estimate, post noi: logit y x c

	local obs = _N
	frame metadata:{
		replace effect = _b[x] if _n == `k'
		replace se = _se[x] if _n == `k'
		replace study = `k' if _n == `k'
		replace size = `obs' if _n == `k'
	}
}

** do the analysis for all studies with complete data
forv i = 1/3{
	use study_`i', replace
	logit y x c
	local obs = _N
	frame metadata{
		replace effect = _b[x] if _n == `i'
		replace se = _se[x] if _n == `i'
		replace study = `i' if _n == `i'
		replace size = `obs'  if _n == `i'
	}
}

frame metadata: meta set effect se , studysize(size)
frame metadata: meta summarize,  eform random(reml)

********************************************************************************
exit 
