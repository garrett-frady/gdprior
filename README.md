
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gdprior

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/garrett-frady/gdprior.svg?branch=master)](https://travis-ci.com/garrett-frady/gdprior)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/garrett-frady/gdprior?branch=master&svg=true)](https://ci.appveyor.com/project/garrett-frady/gdprior)
<!-- badges: end -->

R package for “Frady, G., Dey, D.K. and Mohammed, S., 2023. Gaussian and
Diffused-gamma Feature Extraction Applied to Sparse, High Dimensional
Spatio-Temporal Data by Local Modeling.”

Code to perform estimation, feature extraction, and prediction under the
GD prior structure.

## Installation

You can install the development version of gdprior like so:

``` r
# install the package (devtools required)
if(!require(devtools)) install.packages("devtools")
devtools::install_github('garrett-frady/gdprior')
library(gdprior)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(gdprior)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
