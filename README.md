
# 'ncontRanalyser'
<!-- badges: start -->
<!-- badges: end -->

The goal of 'ncontRanalyser' is to ...

UNDER CONSTRUCTION!! UNDER CONSTRUCTION!!

Don't mind me :)


## Installation

You can install the development version of 'ncontRanalyser' like so:

``` r
# requires devtools

# then use
install_github("sproc2104/ncontRanalyser")

```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library('ncontRanalyser')
## basic example code
```

This is the second case, in which a reference structure is not used to generate reference quartiles, but these are provided manually.

Quartiles must be provided in the vector form, where
``` r
# The values within the quartile vectors MUST be in the order below (decreasing).
ab_quartiles<- c(max_val, Q1, Q2, Q3)
ag_quartiles<- c(max_val, Q1, Q2, Q3)
```
These can then be passed to the mutant_analysis control, which will assign antigen and antibody residue quartiles based on these values.
``` r
mut<- mutant_analysis(h_file, k_file, ab_quartiles, ag_quartiles)
```

The sidnificance of the quartiles is arbitrary, but in previous analysis, I have selected those in the top quartile and above as composing an AlphaFold-predicted epitope (AFPE).
