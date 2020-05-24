[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/scoringRules)](https://cran.r-project.org/package=scoringRules) 
[![Downloads](https://cranlogs.r-pkg.org/badges/scoringRules)](https://cranlogs.r-pkg.org/badges/scoringRules)

# scoringRules 

An R package to compute scoring rules for fixed (parametric) and simulated forecast distributions. Authored by Alexander Jordan (University of Bern), Fabian Krüger (Karlsruhe Institute of Technology (KIT)) and Sebastian Lerch (KIT and Heidelberg Institute for Theoretical Studies), with contributions from Maximiliane Graeter (KIT). 

## Highlights
  - Coherent, dictionary-like reference for computing scoring rules in a wide range of situations
  - Previously unavailable closed-form expressions of the CRPS for many parametric distributions
  - Efficient implementation thanks to R/Rcpp 
  - Whenever more than one implementation variant exists, we offer statistically principled default choices
  
## Installation

CRAN version:
```r
install.packages("scoringRules")
```

Development version (GitHub):
```r
# install.packages("devtools")
library(devtools)
install_github("FK83/scoringRules")
```

## Background

Scoring rules are functions S(F, y) which evaluate the accuracy of a forecast distribution F, given that an outcome y was observed. The **scoringRules** package contains functions to compute scoring rules, for a variety of  distributions F that come up in applied work, and two popular choices of S. Two main classes of distributions are

  - Parametric distributions like normal, t, and gamma. For example, most weather forecasts (which apply statistical postprocessing to physical models) take such a form. 
  - Distributions that are not known analytically, but are indirectly described through a sample of simulaton draws. For example, Bayesian forecasts produced via Markov Chain Monte Carlo (MCMC) take this form. 

The scoring rules we cover are the continuous ranked probability score (CRPS; Matheson and Winkler, *Management Science*, 1976) and the logarithmic score (Good, *Journal of the Royal Statistical Society (Series B)*, 1952). The package further provides functions to compute the multivariate energy and variogram scores for forecast distributions given by discrete samples.

## History
  - August 2019: Version 1.0.0 on CRAN, vignette published in the *Journal of Statistical Software*
  - November 2017: Version 0.9.4 on CRAN, including a detailed vignette 
  - July 2017: Vignette for closed-form expressions of the CRPS (Alexander Jordan), and functions for CRPS-based fitting of truncated/censored distributions
  - July 7, 2016: Version 0.9 published on CRAN
  - September 15, 2014: First commit 
