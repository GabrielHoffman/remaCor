# Random effect meta-analysis for correlated test statistics

Random effect meta-analysis for correlated test statistics using RE2C

## Usage

``` r
RE2C(beta, stders, cor = diag(1, length(beta)), twoStep = FALSE)
```

## Arguments

- beta:

  regression coefficients from each analysis

- stders:

  standard errors corresponding to betas

- cor:

  correlation matrix between of test statistics. Default considers
  uncorrelated test statistics

- twoStep:

  Apply two step version of RE2C that is designed to be applied only
  after the fixed effect model.

## Value

- stat1::

  statistic testing effect mean

- stat2::

  statistic testing effect heterogeneity

- RE2Cp::

  RE2 p-value accounting for correlelation between tests

- RE2Cp.twoStep::

  two step RE2C test after fixed effect test. Only evaluated if
  `twoStep==TRUE`

- QE::

  test statistic for the test of (residual) heterogeneity

- QEp::

  p-value for the test of (residual) heterogeneity

- Isq::

  I^2 statistic

## Details

Perform random effect meta-analysis for correlated test statistics using
RE2 method of Han and Eskin (2011) , or RE2 for correlated test
statistics from Han et al. (2016) . Also uses RE2C method of Lee et al.
(2017) to further test for heterogenity in effect size. By default,
correlation is set to identity matrix to for independent test
statistics.

This method requires the correlation matrix to be symmatric positive
definite (SPD). If this condition is not satisfied, results will be NA.
If the matrix is not SPD, there is likely an issue with how it was
generated.

For computing `Q`, `QEp`, and `Isq`, the effective number of independent
studies is determined by `sum(solve(cor))` (Liu and Liang 1997) . When
the correlation is diagonal, the studies are independent and this value
is simply the number of studies.

## References

Han B, Duong D, Sul JH, de Bakker PI, Eskin E, Raychaudhuri S (2016). “A
general framework for meta-analyzing dependent studies with overlapping
subjects in association mapping.” *Human Molecular Genetics*, **25**(9),
1857–1866. <https://doi.org/10.1093/hmg/ddw049>.  
  
Han B, Eskin E (2011). “Random-effects model aimed at discovering
associations in meta-analysis of genome-wide association studies.” *The
American Journal of Human Genetics*, **88**(5), 586–598.
<https://doi.org/10.1016/j.ajhg.2011.04.014>.  
  
Lee CH, Eskin E, Han B (2017). “Increasing the power of meta-analysis of
genome-wide association studies to detect heterogeneous effects.”
*Bioinformatics*, **33**(14), i379–i388.
<https://doi.org/10.1093/bioinformatics/btx242>.  
  
Liu G, Liang K (1997). “Sample size calculations for studies with
correlated observations.” *Biometrics*, 937–947.
<https://doi.org/10.2307/2533554>.

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

# Run fixed effects meta-analysis, 
# assume identity correlation  
LS( df$beta, df$se)
#>        beta       se           p
#> 1 0.4004915 0.147962 0.006795195

# Run random effects meta-analysis,
# assume identity correlation  
RE2C( df$beta, df$se)
#>      stat1 stat2      RE2Cp RE2Cp.twoStep       QE       QEp Isq
#> 1 7.326322     0 0.01255951            NA 3.585407 0.6105057   0

# Run fixed effects meta-analysis, 
# account for correlation 
LS( df$beta, df$se, C)
#>        beta        se          p
#> 1 0.4657671 0.1684852 0.00570209

# Run random effects meta-analysis,
# account for correlation 
RE2C( df$beta, df$se, C)
#>      stat1     stat2       RE2Cp RE2Cp.twoStep       QE       QEp      Isq
#> 1 7.642119 0.3393766 0.009417392            NA 3.780033 0.3747031 4.994907
```
