deconvolve
=======================

The R package *deconvolve* provides tools for performing non-parametric deconvolution on measurement error problems.

## Installation

You can install the **development** version from [Github](https://github.com/timothyhyndman/deconvolve).

```s
install.packages("devtools")
devtools::install_github("timothyhyndman/deconvolve")
```

## Usage

```s
library(deconvolve)
library(ggplot2)

# Example 1

# Example 2
```

## List of functions

### Exported

#### Deconvolution

* fdecUknown
* fdecUknownhet
* DeconErrSymPmf
* DeconErrSymPmfToPdf

#### Bandwidth Selection

* CVdeconv
* PI_deconvUknownth4
* PI_deconvUKnownth4het

#### Nonparametric Regression With Errors in Variables

* NWDecUKnown
* hSIMEXUknown (bandwidth)

#### Other

* GenerateTestData (convenience function for trying out tests)
* Plot (all functions for quickly plotting results)
* rlap (used in examples)

### Not Exported

* ComputePhi (used in Symmetric Error)
* deconvolve.R (documentation only)
* KernelWeight (used in Symmetric Error)
* outerop (used in various functions - I'm not convinced it's necessary)
* phiK2 (Default kernel characteristic function for when none is supplied)
* PhiUSpline (used in Symmetric Error PmfToPdf)
* PI_DeconvUEstTh4 (used in Symmetric Error PmfToPdf)
* NWDecridgeL10CUknown (used in hSIMEXUknown I think?)

## License

This package is free and open source software, licensed under GPL (>=2).