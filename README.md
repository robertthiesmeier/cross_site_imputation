# Cross-site imputation to recover covariates without sharing individual-level pooled data
In the file below, we work you through an example code of the analysis presneted in this paper (add link here) :page_facing_up:. All datasets are simulated.

## Content 
Missing data is a common challenge across scientific disciplines. Current imputation methods require the availability of individual data to impute missing values. Often, however, missingness requires using external data for the imputation. Therefore, we introduce a new Stata command, `mi impute from`, designed to impute missing values using linear predictors and their related covariance matrix from imputation models estimated in one or multiple external studies. This allows for the imputation of any missing values without sharing individual data between studies. 

## The command: `mi impute from` :computer:
To impute missing data across study sites, we use the Stata package `mi impute from`. The command can be downloaded from the SSC Archive in Stata:

```ruby
ssc install mi_impute_from
```
[!NOTE]
In this [preprint](https://arxiv.org/pdf/2410.02982v1), we describe the underlying method and present the syntax of `mi impute from` alongside practical examples of missing data in collaborative research projects.

# Cross-site imputation 

## How does it work? :pencil2:

For simplicity, we consider three binary variables: an outcome, Y, and exposure, X, and a confounder C. The relations between the three variables is as follows: 

![DAG_conf](https://github.com/user-attachments/assets/f5a3d195-dc43-4d09-b147-5d467fb01b04)

To generate the data according to the outlined mechanism above, apply the following program:

<details>
**Data Generating Mechanism (DGM)**

To run the code with the simulated data in this example, please check the DGM. 

  ```ruby

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

```

We generate five studies according to the specified mechanism above: 

```ruby

forv i = 1/5 {
	singlestudy, study(`i') filename(study_`i')
}

```

</details>

### Control approach: No missing data :passport_control:

This analysis serves as a control approach to evaluate the performance of the imputation approach. Can we successfully recover (i.e., impute) missing covariates in a single example? 
In the following examples, we work with `frames` in Stata. Even though we generate all data ourselves, we cannot pool the individual data to make the example as realistic as possible. 
We can initiate the dataframe fr our analysis in Stata as follows:

```ruby

capture frame drop metadata
frame create metadata 
qui frame metadata:{ 
	set obs 5 // number of studies
	gen effect = .
	gen se = .
	gen study = .
	gen size = .
}  

```

Now, for each study in our data network, we perform the same analysis (federated analysis) and save the estimates for each study in the frame we initiated in the previous step.

```ruby

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

```

Finally, the estimates can be meta-analysed with a fixed or random meta-analytical model.

```ruby

frame metadata: list 
frame metadata: meta set effect se , studysize(size)
frame metadata: meta summarize,  eform fixed

```

#### Generate missing data to apply `mi impute from` :exclamation::heavy_exclamation_mark:

After having conducted the control analysis, we can generate systematically missing data on the confounder C in Study 4 and 5. 

```ruby

forv i = 4/5{
	qui use study_`i', clear 
	qui count 
	di in red "Study `i': N = " r(N)
	replace c = . 
	save study_`i', replace
}

```

### Approach :one: Federated analysis using 5 studies without adjustment for C

In this analysis, we consider all studies, however, none of them incldues the adjustment for the confounder C.

```ruby

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

```
### Approach :two: Federated analysis using only studies with complete informtion on all variables

This approach takes into consideration 3/5 studies and disregards study 4 and 5 as a result of systematically missing data on C at these sites.

```ruby

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

```

### Approach :three: Federated analysis using all studies with complete and incomplete informtion on all variables

In this approach we aim to include all studies in the analysis, regardless whether or not we are able to adjust for C in some of the studies with missing information. 
First, we fit the fully-adjusted outcome model in study 1 to 3: 

```ruby

forv i = 1/3{
	use study_`i', replace
	logit y x c
	local obs = _N
	frame metadata: replace effect = _b[x] if _n == `i'
	frame metadata: replace se = _se[x] if _n == `i'
	frame metadata: replace study = `i' if _n == `i'
	frame metadata: replace size = `obs'  if _n == `i'
}

```

Second, we fit the model not adjusting for the confounding variable C in study 4 and 5. 

```ruby

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

```

Finally, we can again analyse all studies with a meta-analytical model.

```ruby

frame metadata: meta set effect se , studysize(size)
frame metadata: meta summarize, eform fixed

```

### Approach :four: Federated analysis using all studies and recovering the missing variable C.

The next two approaches consider all studies and applying cross-site imputation to recoevr the variables with missing informations at sites where values are missing. In this approach, we consider a randomly selected study that has complete information on C (one out of the three studies) and fir the imputation model in that study. The imputation regression coefficients are then exported to text files that can be easily shared with other studies.

```ruby

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

```

At the study sites with missing data (Study 4 and Study 5), we import the regression coefficients and impute the issing values on C 10-times. A detailed explanation on this step can be found here (add link to im impute from). Finally, we fit the outcome model to each imputed dataset and combine the estimates with Rubin's rules.

```ruby

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

```

We have now estimates the C-adjusted effect of the exposure X on the outcome Y after imputing C in Study 4 and 5. Next, we can derive the estimates for Study 1 to 3 as we have done the in the previous steps and appky a meta-analysis to derive a single pooled estimate.

```ruby

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

```

### Approach :five: Federated analysis using all studies and recovering the missing variable C.

In this approach, we also consider all five studies using cross-site imputation to recover missing values of C in Study 4 and 5. However, here we consider Study 1 to 3 to fit an imputaton model as opposed to only considering a single study as the basis for imputation. To do so, we first fit the imputation model in the studies with available data on the confounder C. 

```ruby

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

```

We save all files and transport them to the sites with missing data. Here, we proceed as outlined in Approach 4. The command `mi_impute_from_get` recognises multiple input files and takes a weighted average. 

```ruby

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

```

We add the estimates for the studies with complete data and apply a meta-analysis in the end. 

```ruby

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

```

## Wrap-up :white_check_mark:
All five approaches can be implemented without the need for any real data and you can test the package `mi impute from`. Again, please refer to this preprint (add link here) for a detailed description of the steps and assumptions made that are pivotal to understand the concept of cross-site imputation.

## Further reading
:label: The underlying imputation method and a simulation study are described in: [Thiesmeier, R., Bottai, M., & Orsini, N. (2024). Systematically missing data in distributed data networks: multiple imputation when data cannot be pooled. Journal of Statistical Computation and Simulation, 1–19](https://doi.org/10.1080/00949655.2024.2404220)

:label: The documentation of `mi impute from` is available here: [Thiesmeier, R., Bottai, M., & Orsini, N. (2024). Imputing Missing Values with External Data](https://arxiv.org/pdf/2410.02982) 

:label: The first version of `mi impute from` was presented at the [2024 UK Stata Conference in London](https://www.stata.com/meeting/uk24/slides/UK24_Orsini.pdf)

:warning: If you find any error please notfiy us: robert.thiesmeier@ki.se



