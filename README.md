# `lmpermutation`

This package provides functions to compute permutation tests in linear models with nuisances variables. The package is focused on several goals :

* Providing to users the most recent methods to handle nuisance variables for permutation tests in the linear models.
* Giving to users tools to compute most common tests in linear model (t test, ANOVA and repeated measure ANOVA).
* Providing an extension for the multiple comparisons problems in linear models with a focus for EEG data.

## The `lmperm()` function
This function is constructed as an extension of the the `lm()` function for permutation test. It produce t statistics with univariate and bivariate pvalue by permutation.

## The `aovperm()` function
This function is constructed as an extension of the the `aov()` function for permutation test. It produces marginal F statistics (type III). Repeated measures anova are feasible using the same notations used in an `aov()` formula with `+Error(id/within)` to specify the random effects.

## The `clusterlm()` function
This function compute cluster-mass statistics for multiple comparisons. It is designed for ERP analysis of unichannel EEG data. The left part of formula object must be a matrix or dataframe which columns represents multiple responses tested on the same experimental design (specified by right part of the formula). This function provides several methods to handle nuisance variables, a F or t statistics, an extension for repeated measure anova and several methods for the multiple comparisons lit the threshold-free cluster enhancement. 

# Contact
If you need help to use the package or want to report errors, contact Jaromil Frossard at <jaromil.frossard@unige.ch>.



# References
For permutation tests with nuisance variables :

* Kherad-Pajouh, S., & Renaud, O. (2010). An exact permutation method for testing any effect in balanced and unbalanced fixed effect ANOVA. Computational Statistics & Data Analysis, 54(7), 1881-1893.
* Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M., & Nichols, T. E. (2014). Permutation inference for the general linear model. Neuroimage, 92, 381-397.

For permutation test in repeated measure ANOVA :

* Kherad-Pajouh, S., & Renaud, O. (2015). A general permutation approach for analyzing repeated measures ANOVA and mixed-model designs. Statistical Papers, 56(4), 947-967.

For cluster-mass statistics for the muliple comparison problems :

* Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG-and MEG-data. Journal of neuroscience methods, 164(1), 177-190.

For the Threshold-free cluster enhancement method :

* Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster enhancement: addressing problems of smoothing, threshold dependence and localisation in cluster inference. Neuroimage, 44(1), 83-98.




