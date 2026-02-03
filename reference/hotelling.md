# Hottelling T^2 test for multivariate regression

Hottelling T^2 test compares estimated regression coefficients to
specified values under the null. This tests a global hypothesis for all
specified coefficients. It uses the F-distribution as the null for the
test statistic which is exact under finite sample size.

## Usage

``` r
hotelling(beta, Sigma, n, mu_null = rep(0, length(beta)))
```

## Arguments

- beta:

  regressioin coefficients

- Sigma:

  covariance matrix of regression coefficients

- n:

  sample size used for estimation

- mu_null:

  the values of the regression coefficients under the null hypothesis.
  Defaults to all zeros

## Details

The Hotelling T2 test is not defined when n - p \< 1. Returns
`data.frame` with `stat = pvalue = NA`.

## Examples

``` r
library(clusterGeneration)
library(mvtnorm)

# sample size
n = 30

# number of response variables
m = 2

# Error covariance
Sigma = genPositiveDefMat(m)$Sigma

# regression parameters
beta = matrix(.6, 1, m)

# covariates
X = matrix(rnorm(n), ncol=1)

# Simulate response variables
Y = X %*% beta + rmvnorm(n, sigma = Sigma)

# Multivariate regression
fit = lm(Y ~ X)

# extract coefficients and covariance
# corresponding to the x variable
beta = coef(fit)['X',]
S = vcov(fit)[c(2,4), c(2,4)]

# perform Hotelling test
hotelling(beta, S, n)
#>    n p      tsq     stat    p.value
#> 1 30 2 6.502292 3.139038 0.05888199
```
