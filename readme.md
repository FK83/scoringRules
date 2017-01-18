# scoringRules 

An R package to compute scoring rules for fixed (parametric) and simulated forecast distributions. Authored by Alexander Jordan (Heidelberg Institute for Theoretical Studies, HITS), Fabian Krüger (Heidelberg University) and Sebastian Lerch (HITS and Karlsruhe Institute of Technology, KIT), with contributions from Maximiliane Graeter (KIT). 

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

Scoring rules are functions S(F, y) which evaluate the accuracy of a forecast distribution F, given that an outcome y was observed. The **scoringRules** package contains functions to compute scoring rules, for a variety of distributions F that come up in applied work, and two popular choices of S. Two main classes of distributions are

  - Parametric distributions like normal, t, and gamma. For example, most weather forecasts (which apply statistical postprocessing to physical models) take such a form. 
  - Distributions that are not known analytically, but are indirectly described through a sample of simulaton draws. For example, Bayesian forecasts produced via Markov Chain Monte Carlo (MCMC) take this form. 

The scoring rules we cover are the continuous ranked probability score (CRPS; Matheson and Winkler, 1976) and the logarithmic score (Good, 1952). The package further provides functions to compute the multivariate energy and variogram scores for forecast distributions given by discrete samples.

## History
  - January 2017: Added multivariate energy and variogram score
  - January 11, 2017: Version 0.9.2 published on CRAN
  - December 2016: Small fixes for truncated/censored distributions
  - September 6, 2016: Version 0.9.1 published on CRAN 
  - August 2016: Added replication materials for our paper on MCMC based forecasting 
    (see <http://arxiv.org/abs/1608.06802>)
  - July 7, 2016: Version 0.9 published on CRAN
  - April -- June 2016: Various small changes (improved consistency and better documentation)
  - August 3, 2015: Some name changes to header functions
  - February 6, 2015: Re-design of internal function structure
  - November 19, 2014: Split functions according to parametric versus simulated forecast distributions
  - September 15, 2014: First commit 
