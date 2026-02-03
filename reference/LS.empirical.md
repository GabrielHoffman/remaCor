# Fixed effect meta-analysis for correlated test statistics

Fixed effect meta-analysis for correlated test statistics using the
Lin-Sullivan method using Monte Carlo draws from the null distribution
to compute the p-value.

## Usage

``` r
LS.empirical(
  beta,
  stders,
  cor = diag(1, length(beta)),
  nu,
  n.mc.samples = 10000,
  seed = 1,
  useGamma = TRUE
)
```

## Arguments

- beta:

  regression coefficients from each analysis

- stders:

  standard errors corresponding to betas

- cor:

  correlation matrix between of test statistics. Default considers
  uncorrelated test statistics

- nu:

  residual degrees of freedom

- n.mc.samples:

  number of Monte Carlo samples

- seed:

  random seed so results are reproducable

- useGamma:

  if `TRUE`, use gamma approximation to fit empirical distribution of
  test statistics and compute p-value. Ff `FALSE`, report p-value as
  `(1 + sum(obs > stat)) / (length(stat)+1)`. Using `TRUE` allows
  computation of small p-values with fewer Monte Carlo samples.

## Details

The theoretical null for the Lin-Sullivan statistic for fixed effects
meta-analysis is chisq when the regression coefficients are estimated
from a large sample size. But for finite sample size, this null
distribution is not well characterized. In this case, we are not aware
of a closed from cumulative distribution function. Instead we draw
covariance matrices from a Wishart distribution, sample coefficients
from a multivariate normal with this covariance, and then compute the
Lin-Sullivan statistic. A gamma distribution is then fit to these draws
from the null distribution and a p-value is computed from the cumulative
distribution function of this gamma.

## References

Hoffman GE, Roussos P (2025). “Fast, flexible analysis of differences in
cellular composition with crumblr.” *biorxiv*.
<https://doi.org/10.1101/2025.01.29.635498>.

## See also

[`LS()`](http://gabrielhoffman.github.io/remaCor/reference/LS.md)

## Examples

``` r
library(clusterGeneration)
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

# Meta-analysis assuming infinite sample size
# but the p-value is anti-conservative
LS(df$beta, df$se, C)
#>        beta         se            p
#> 1 0.5974489 0.08610924 3.969383e-12

# Meta-analysis explicitly modeling the finite sample size
# Gives properly calibrated p-values
# nu is the residual degrees of freedom from the model fit
LS.empirical(df$beta, df$se, C, nu=n-2)
#>        beta         se            p
#> 1 0.5974489 0.08610924 9.033952e-08
```
