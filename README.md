
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

License
-------

This package is free and open source software, licensed under GPL (&gt;=2).
