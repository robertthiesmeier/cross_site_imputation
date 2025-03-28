# Cross-site imputation to recover covariates without in federated analyses without sharing individual-level pooled data
On this site, we illustrate computer code to apply cross-site imputation. A preprint is available [here](https://www.medrxiv.org/content/10.1101/2024.12.19.24319364v1.full.pdf). :page_facing_up:. All datasets are simulated, and you can try out the code and implementation of the software package yourself.

## What is the context? 
Missing data is a common challenge across scientific disciplines. Most of the currently implemnted imputation methods require the availability of individual data to impute missing values. Often, however, missingness requires using external data for the imputation. Therefore, propose a new imputation approach - cross-site multiple imputation - designed to impute missing values using linear predictors and their related covariance matrix from imputation models estimated in one or multiple external studies. This allows for the imputation of any missing values without sharing individual data between studies. The idea was previously discussed [here](https://www.tandfonline.com/doi/full/10.1080/00949655.2024.2404220). In this short tutorial on cross-site imputation, we will work with the newly developed Stata code `mi impute from` that facilitates the imputation of missing values. 

### Download `mi impute from` :computer:
To impute missing data across study sites, we use the Stata package `mi impute from`. The command can be downloaded from the SSC Archive in Stata:

```ruby
ssc install mi_impute_from
```
> [!NOTE]
> In this [preprint](https://arxiv.org/pdf/2410.02982v1), we describe the underlying method and present the syntax of `mi impute from` alongside practical examples of missing data in collaborative research projects.

# Cross-site imputation 

## How does it work? :pencil2:

For simplicity, we consider three binary variables: an outcome, Y, and exposure, X, and a confounder C. The relations between the three variables is as follows: 

![DAG_conf](https://github.com/user-attachments/assets/f5a3d195-dc43-4d09-b147-5d467fb01b04)

> [!TIP]
> To generate the data according to the outlined mechanism above and NOT use any real data, apply the DGM program below.

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
We can initiate the dataframe for our analysis in Stata as follows:

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

#### Generate missing data to apply `mi impute from` :heavy_exclamation_mark:

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

In this analysis, we consider all studies, however, none of them includes the adjustment for the confounder C.

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
### Approach :two: Federated analysis using only studies with complete information on all variables

This approach takes into consideration 3/5 studies and disregards study 4 and 5 as a result of systematically missing data on C at these sites.

```ruby

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

```

### Approach :three: Federated analysis using all studies with complete and incomplete information on all variables

In this approach we aim to include all studies in the analysis, regardless of whether or not we are able to adjust for C in some of the studies with missing information. 
First, we fit the fully adjusted outcome model in study 1 to 3: 

```ruby

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

The next two approaches consider all studies and applying cross-site imputation to recover the variables with missing information at sites where values are missing. In this approach, we consider a randomly selected study that has complete information on C (one out of the three studies) and fir the imputation model in that study. The imputation regression coefficients are then exported to text files that can be easily shared with other studies.

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

At the study sites with missing data (Study 4 and Study 5), we import the regression coefficients and impute the missing values on C 10-times. A detailed explanation on this step can be found [here](https://github.com/robertthiesmeier/mi_impute_from.git). Finally, we fit the outcome model to each imputed dataset and combine the estimates with Rubin's rules.

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

We have now estimated the C-adjusted effect of the exposure X on the outcome Y after imputing C in Study 4 and 5. Next, we can derive the estimates for Study 1 to 3 as we have done the in the previous steps and apply a meta-analysis to derive a single pooled estimate.

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

In this approach, we also consider all five studies using cross-site imputation to recover missing values of C in Study 4 and 5. However, here we consider Study 1 to 3 to fit an imputaton model as opposed to only considering a single study as the basis for imputation. To do so, we first fit the imputation model in the studies with available data on the confounder C. After the first step (fitting the prediction model in the studies with available data), we save all files and transport them to the sites with missing data. Here, we proceed as outlined in Approach 4. The command `mi_impute_from_get` recognises multiple input files and takes a weighted average. 


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

## :high_brightness: Extension: Empirical heterogenity between study sites
When multiple studies are used to fit the prediction model, it is often desirable to account for the statistical (i.e., empirical) heterogenity between the study sites. To account for these differences in the final imputation model, we can fit a meta-regression model with random effects combining the regression coefficients from multiple studies. The example below illustrates this approach. 

> [!TIP]
> In the following section, we walk you through the code. For smooth implemntation, refer to the do.file and run everything in one go to avoid any error messages.

First, as shown above, the imputation must be fit at all sites with complete data on the systematically missing confounder.

```ruby

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

```

📧 We can send the .txt files to the sites with missing data. 
At the receiving site, we then have to transform the files into matrices. 

```ruby

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

```

▶️ We now imported all regression coefficients and their variance/covariances from the three prediction models. 
Next, we would like to take a random meta regression model of the three set of coefficients to derive our final imputation regression coefficients. In the previous examples, when multiple files were input, `mi_impute_from_get` faciliated a weighted average using inverse variance method. Here, we show how to facililate a random meta regression model to respect the statistical heterogenity between sites. (🕵️ Expand to see code shown in details)

<details>
	
```ruby

capture frame drop random_impmodel  
frame create random_impmodel study y1 y2 y3 v11 v12 v13 v22 v23 v33

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

```
</details>

Once we created a dataframe and organised the beta-coefficients and their variance/covariances, we can fit the meta regression model with random effects. 

```ruby

frame random_impmodel: meta mvregress y* , wcovvariables(v*) random(reml, covariance(identity))

```

The output has to be cleaned up and organised a little bit to be used for `mi impute from`. (🕵️ Expand to see code shown in details)

<details>
	
```ruby
	
** clean up the results to be used with mi impute from
mat ib = e(b)[1, 1..3]
mat colnames ib = y x _cons

mat iV = e(V)[1..3, 1..3]
mat colnames iV = y x _cons
mat rownames iV = y x _cons

mat li ib 
mat li iV

```

</details>

We can then proceed as previously, and use the final set of imputation coefficients to be used with `mi impute from`.

```ruby

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

** do the analysis for all studies with complete data and fit a final random-effects meta-analysis

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

```

☑️ We now used the cross-site imputation approach allowing for heterogenity in the imputation model.

## Wrap-up :white_check_mark:
All five approaches can be implemented without the need for any real data and you can test the package `mi impute from`. In addition, we showed how to incorporate empirical heterogenity between the sites into the final imputation model by fitting a meta-regression model with random effects on the regression coefficients coming from multiple sites. This approach is explained in more detail in Resche-Rigon et al. (2018). 


> [!IMPORTANT]
> Again, please refer to this [preprint](https://www.medrxiv.org/content/10.1101/2024.12.19.24319364v1.full.pdf) for a more detailed description of the steps and assumptions made that are pivotal to understand the concept of cross-site imputation.

## Further reading
:label: The idea of cross-site imputation, an applied example, and a discussion around the assumptions of the method are presented in: [Thiesmeier, R. Madley-Dowd, P. Orsini, N., Ahlqvist, V. (2024). Cross-site imputation for recovering variables without individual pooled data.](https://www.medrxiv.org/content/10.1101/2024.12.19.24319364v1.full.pdf+html)

:label: The underlying imputation method and a simulation study are described in: [Thiesmeier, R., Bottai, M., & Orsini, N. (2024). Systematically missing data in distributed data networks: multiple imputation when data cannot be pooled. Journal of Computational Statistics and Simulation, 94(17), 3807–3825](https://doi.org/10.1080/00949655.2024.2404220)

:label: The documentation of `mi impute from` is available here: [Thiesmeier, R., Bottai, M., & Orsini, N. (2024). Imputing Missing Values with External Data.](https://arxiv.org/pdf/2410.02982) 

:label: The software package `mi impute from` in Stata is stored here: [Thiesmeier R, Bottai M, Orsini R. (2024). MI_IMPUTE_FROM: Stata module to impute using an external imputation model. Statistical Software Components S459378, Boston College Department of Economics](https://ideas.repec.org/c/boc/bocode/s459378.html)

:label: The first version of `mi impute from` was presented at the [2024 UK Stata Conference in London.](https://www.stata.com/meeting/uk24/slides/UK24_Orsini.pdf)

:label: Cross-site imputation was presented at the Royal Statistical Society International Conference in Brighton, UK, in September 2024, and at the International Biometric Society Conference in Atlanta, GA, USA, in December 2024.

:label: [Multiple imputation by chained equations for systematically and sporadically missing multilevel data](https://journals.sagepub.com/doi/10.1177/0962280216666564?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) by Resche-Rigon M. and White I. illustarting the theoretical foundations of a two-stage imputation process for multilevel data. 

:label: A first introduction and illustration of [MI algorithms in distributed health data networks](https://www.nature.com/articles/s41467-020-19270-2) by Chang et al. (2020). 

:warning: If you find any errors, please notfiy us: robert.thiesmeier@ki.se
