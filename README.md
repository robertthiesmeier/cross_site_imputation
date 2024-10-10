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

For simplicity, we consider three binary variables: an outcome, Y, and exposure, X, and a confounder C. The relations between the three variables is as follows: 

![DAG_conf](https://github.com/user-attachments/assets/f5a3d195-dc43-4d09-b147-5d467fb01b04)

To generate the data according to the outlined mechanism above, apply the following program:

<details>
## Data Generating Mechanism (DGM)

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

### Analysis 1: No missing data (control)

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

#### Generate missing data to apply `mi impute from`

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



