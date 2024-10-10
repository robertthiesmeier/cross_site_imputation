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
