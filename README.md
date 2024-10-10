# Cross-site imputation to recover covariates without sharing individual-level pooled data
In the file below, we work you through an example code of the analysis presneted in this paper (add link here) :page_facing_up:. All datasets are simulated.

## Content 
Missing data is a common challenge across scientific disciplines. Current imputation methods require the availability of individual data to impute missing values. Often, however, missingness requires using external data for the imputation. Therefore, we introduce a new Stata command, `mi impute from`, designed to impute missing values using linear predictors and their related covariance matrix from imputation models estimated in one or multiple external studies. This allows for the imputation of any missing values without sharing individual data between studies. 

## The command: `mi impute from` :computer:
To impute missing data across study sites, we use the Stata package `mi impute from`. The command can be downloaded from the SSC Archive in Stata:

```ruby
ssc install mi_impute_from
```

In this [preprint](https://arxiv.org/pdf/2410.02982v1), we describe the underlying method and present the syntax of `mi impute from` alongside practical examples of missing data in collaborative research projects.

# Cross-site imputation 

## How does it work? :pencil2:

<details>
  The data used for this example can be generated with this program: 

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
	
	local pc =  runiform(0.1, 0.5)
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
	
	// no adjustment for c
	logit y x, or nolog
	ret scalar b1_mod1 = _b[x]
	ret scalar b1_se_mod1 = _se[x]
	
	// adjustment for c
	logit y x c, or nolog
	ret scalar b1_mod2 = _b[x]
	ret scalar b1_se_mod2 = _se[x]
	test x 
	ret scalar pvalue = r(p)
	
end

```
</details>
