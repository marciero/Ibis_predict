# A framework for Bayesian modeling of hospital visits 

## Overview

We combine real and simulated data for illustration of a modeling framework for  
the number of hospital and ER visits for patients, and subsequently cost and other statistics of interest, for users  of 
Senscio Systems Ibis tablet device for self care management. We use data such as 

* demographic - age, gender, geographic (city, town, zip code)
* socioeconomic- income or financial class strata, 
* health data - ICD 10  and procedure codes, medications, etc. 


 With many thousands of predictor variables, we can perform feature engineering and selection to reduce the size of the predictor set. Such techniques include

* regularization
* likelihood encoding
* principle component analysis


## Modeling approaches

We propose to model `visits` as a function of the other predictors using Poisson or Negative Binomial regression, which are both appropriate for count data. The latter is often used with over or under dispersed data because there is an additional shape parameter to model. 

We primarily use a Bayesian approach where possible.

* models likelihood of parameters given the data (which we know), rather the other way around as with null-hypotheis-significane-test (NHST) type approaches
* distributions rather than point estimates, for all parameters, not just outcome variable
* ability to model any function of parameters or outcome variables, such as cost for given number of visits.
* ability to easily examine changes in outcome with change in a single predictor; eg how the $distribution$ of `visits` varies between users and non users of Ibis for patients of a given gender, financial class, medical profile, etc. 
* computationally expensive

### Condition categories

We can use condition categories such as the CMS HCC (hierarchical condition categories) to map ICD 10 codes to a reduced set of predictors. HCC are based on clinical relationships of diagnoses, and are used to perform risk adjustment for Medicare Advantage plans. Potential abuses from "gaming" the coding are not an issue if we do the mapping ourselves, assuming "honest" ICD 10 diagnoses. Gaming is less an issue also because the HCCs will be predictors in a model, and not used to generate a cost directly. 

CMS has tables defining the ICD 10 $\rightarrow$ HCC assignments and I've done a partial mapping manually on a sample data set. We may also be able to obtain data from Chronic Conditions Warehouse, an entity contracted by CMS for the purpose of curating Medicare patient data for researchers. They can provide software which outputs and single risk score for a given patient. In fact, the CMS risk score actually a form of *likelihood encoding* of categorical variables. 




