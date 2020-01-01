
<!-- README.md is generated from README.Rmd. Please edit that file -->

# reducedscan

<!-- badges: start -->

[![Travis-CI Build
Status](https://travis-ci.org/fboehm/reducedscan.svg?branch=master)](https://travis-ci.org/fboehm/reducedscan)

<!-- badges: end -->

Please note that the ‘reducedscan’ project is released with a
[Contributor Code of Conduct](.github/CODE_OF_CONDUCT.md). By
contributing to this project, you agree to abide by its terms.

The goal of `reducedscan` is to determine allelic series at a QTL and
scan based on that series.

## Installation

You can install the released version of reducedscan from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("reducedscan")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fboehm/reducedscan")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(reducedscan)
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
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
