# Fixed effect meta-analysis for correlated test statistics

Fixed effect meta-analysis for correlated test statistics using the
Lin-Sullivan method.

## Usage

``` r
LS(beta, stders, cor = diag(1, length(beta)))
```

## Arguments

- beta:

  regression coefficients from each analysis

- stders:

  standard errors corresponding to betas

- cor:

  correlation matrix between of test statistics. Default considers
  uncorrelated test statistics

## Value

- beta::

  effect size

- se::

  effect size standard error

- p::

  p-value

## Details

Perform fixed effect meta-analysis for correlated test statistics using
method of Lin and Sullivan (2009). By default, correlation is set to
identity matrix to for independent test statistics.

This method requires the correlation matrix to be symmatric positive
definite (SPD). If this condition is not satisfied, results will be NA.
If the matrix is not SPD, there is likely an issue with how it was
generated.

However, evaluating the correlation between observations that are not
pairwise complete can give correlation matricies that are not SPD. In
this case, consider running `Matrix::nearPD( x, corr=TRUE)` to produce
the nearest SPD matrix to the input.

## References

Lin D, Sullivan PF (2009). “Meta-analysis of genome-wide association
studies with overlapping subjects.” *The American Journal of Human
Genetics*, **85**(6), 862–872.
<https://doi.org/10.1016/j.ajhg.2009.11.001>.

## Examples

``` r
library(clusterGeneration)
#> Loading required package: MASS
library(mvtnorm)

# sample size
n = 30

# number of response variables
m = 6

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

# Correlation between residuals
C = cor(residuals(fit))

# Extract effect sizes and standard errors from model fit
df = lapply(coef(summary(fit)), function(a) 
 data.frame(beta = a["X", 1], se = a["X", 2]))
df = do.call(rbind, df)

# Run fixed effects meta-analysis, 
# assume identity correlation  
LS( df$beta, df$se)
#>        beta        se            p
#> 1 0.6854324 0.1459224 2.637027e-06
 
# Run fixed effects meta-analysis, 
# account for correlation 
LS( df$beta, df$se, C)
#>        beta        se            p
#> 1 0.7280211 0.1596967 5.145312e-06
```
