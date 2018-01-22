
<!-- README.md is generated from README.Rmd. Please edit that file -->
deconvolve
==========

The R package *deconvolve* provides tools for performing non-parametric deconvolution on measurement error problems. It contains functions for finding bandwidths, deconvolved densities and non-parametric regression estimates.

Installation
------------

You can install the **development** version from [Github](https://github.com/timothyhyndman/deconvolve).

``` r
# install.packages("devtools")
devtools::install_github("timothyhyndman/deconvolve")
```

Usage
-----

``` r
library(deconvolve)

n <- 200
W <- GenerateTestData(n, dist_type = "chi", error_type = "norm")
xx <- seq(min(W), max(W), length.out = 100)
d <- deconvolve(W, xx)
```

Available methods
-----------------

This is a list of every method that this package contains

##### Deconvolution

-   Homoscedastic Errors
-   Heteroscedastic Errors
-   Unknown Errors
-   *Unknown Errors with Replicates (coming soon)*

##### Bandwidth selection

-   Plug-in Estimator for Homoscedastic Errors
-   Plug-in Estimator for Heteroscedastic Errors
-   Plug-in Estimator for Unknown Errors
-   *Plug-in Estimator for Unknown Errors with Replicates (coming soon)*
-   CV
-   SIMEX

##### Regression

-   Homoscedastic Errors

License
-------

This package is free and open source software, licensed under GPL (&gt;=2).
