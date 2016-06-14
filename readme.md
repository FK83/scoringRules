# scoringRules 

An R package to compute scoring rules for fixed (parametric) and simulated forecast distributions

## Highlights
  - Coherent, dictionary-like reference for computing scoring rules in a wide range of situations
  - Previously unavailable closed-form expressions of the CRPS for many parametric distributions
  - Efficient implementation thanks to R/Rcpp 
  - Whenever more than one implementation variant exists, we offer statistically principled default choices
  
## Installation

```r
# install.packages("devtools")
library(devtools)
install_github("FK83/scoringRules")
library(scoringRules)
```

## Background

Scoring rules are functions $S(F, y)$ which evaluate the accuracy of a forecast distribution $F$, given that an outcome $y$ was observed. The **scoringRules** package contains functions to compute scoring rules, for a variety of distributions F that come up in applied work, and three popular choices of S. Two main classes of distributions are

  - Parametric distributions like normal, t, and gamma. For example, most weather forecasts (which apply statistical postprocessing to physical models) take such a form. 
  - Distributions that are not known analytically, but are indirectly described through a sample of simulaton draws. For example, Bayesian forecasts produced via Markov Chain Monte Carlo (MCMC) take this form. 

The scoring rules we cover are the continuous ranked probability score (CRPS; Matheson and Winkler, 1976) and the logarithmic score (Good, 1952).

## History
  - April 2016: Various small changes (improved consistency and better documentation)
  - August 3, 2015: Some name changes to header functions
  - February 6, 2015: Re-design of internal function structure
  - November 19, 2014: Split functions according to parametric versus simulated forecast distributions
  - September 15, 2014: First commit 
