# Random effects meta-analysis for correlated test statistics

Standard approaches to meta-analysis assumes that effect sizes are
statistically independent. Here we provide methods for fixed and random
effects meta-analysis when the correlation between effect sizes are
known.

#### Fixed effects meta-analysis

[`LS()`](http://gabrielhoffman.github.io/remaCor/reference/LS.md)
implements fixed effect meta-analysis for correlated test statistics
using method of Lin and Sullivan (2009). By default, correlation is set
to identity matrix to for independent test statistics.

#### Random effects meta-analysis

[`RE2C()`](http://gabrielhoffman.github.io/remaCor/reference/RE2C.md)
implements random effect meta-analysis for correlated test statistics
that jointly tests deviation of the mean from zero as well as effect
size heterogenity. This method uses the RE2 method of Han and Eskin
(2011), or RE2 for correlated test statistics from Han et al. (2016). By
default, correlation is set to identity matrix to for independent test
statistics. (In addition, this function computes the two step RE2C
method of Lee, Eskin, and Han (2017) to further test for heterogenity in
effect size after applying a fixed effect test.)

- `stat1`: statistic testing effect mean

- `stat2`: statistic testing effect heterogeneity

- `RE2Cp`: RE2 p-value accounting for correlelation between tests. (This
  is the p-value appropriate for most questions)

- `RE2Cp.twoStep`: two step RE2C test after fixed effect test. Only
  evaluated if `twoStep==TRUE`. (not typically used)

- `QE`: test statistic for the test of (residual) heterogeneity

- `QEp`: p-value for the test of (residual) heterogeneity

- `Isq`: I^2 statistic

  `QE`, `QEp` and `Isq` are only evaluted if correlation is diagonal

#### Examples

``` r
library(remaCor)
library(metafor)
library(mvtnorm)
library(clusterGeneration )

# sample size
n = 30

# number of response variables
m = 2

# Error covariance
Sigma = genPositiveDefMat(m)$Sigma

# regression parameters
beta = matrix(0, 1, m)

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

# Standard fixed effects meta-analysis
# of independent effects with metafor pacakge
rma( df$beta, sei=df$se, method="FE")
```

    ## 
    ## Fixed-Effects Model (k = 2)
    ## 
    ## I^2 (total heterogeneity / total variability):   0.00%
    ## H^2 (total variability / sampling variability):  0.16
    ## 
    ## Test for Heterogeneity:
    ## Q(df = 1) = 0.1647, p-val = 0.6849
    ## 
    ## Model Results:
    ## 
    ## estimate      se     zval    pval    ci.lb   ci.ub    
    ##  -0.0834  0.2450  -0.3402  0.7337  -0.5636  0.3969    
    ## 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Standard random effects meta-analysis
# of independent effects with metafor pacakge
rma( df$beta, sei=df$se, method="REML")
```

    ## 
    ## Random-Effects Model (k = 2; tau^2 estimator: REML)
    ## 
    ## tau^2 (estimated amount of total heterogeneity): 0 (SE = 0.2444)
    ## tau (square root of estimated tau^2 value):      0
    ## I^2 (total heterogeneity / total variability):   0.00%
    ## H^2 (total variability / sampling variability):  1.00
    ## 
    ## Test for Heterogeneity:
    ## Q(df = 1) = 0.1647, p-val = 0.6849
    ## 
    ## Model Results:
    ## 
    ## estimate      se     zval    pval    ci.lb   ci.ub    
    ##  -0.0834  0.2450  -0.3402  0.7337  -0.5636  0.3969    
    ## 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Run fixed effects meta-analysis, assume identity correlation  
# Use Lin-Sullivan method
LS( df$beta, df$se)
```

    ##          beta        se         p
    ## 1 -0.08335047 0.2450277 0.7337303

``` r
# Run fixed effects meta-analysis, accounting for correlation  
# Use Lin-Sullivan method
LS( df$beta, df$se, C)
```

    ##         beta        se         p
    ## 1 -0.1014787 0.1857379 0.5848224

``` r
# Run random effects meta-analysis, assume identity correlation  
RE2C( df$beta, df$se)
```

    ##      stat1 stat2     RE2Cp RE2Cp.twoStep        QE      QEp Isq
    ## 1 0.115714     0 0.7667089            NA 0.1646886 0.684876   0

``` r
# Run random effects meta-analysis, accounting for correlation 
RE2C( df$beta, df$se, C)
```

    ##       stat1 stat2    RE2Cp RE2Cp.twoStep        QE       QEp Isq
    ## 1 0.2985031     0 0.630749            NA 0.1701623 0.9702919   0

    ##       stat1 stat2    RE2Cp RE2Cp.twoStep        QE       QEp Isq
    ## 1 0.2985031     0 0.630749            NA 0.1701623 0.9702919   0

## Session info

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] clusterGeneration_1.3.8 MASS_7.3-65             mvtnorm_1.3-3          
    ## [4] metafor_4.8-0           numDeriv_2016.8-1.1     metadat_1.4-0          
    ## [7] Matrix_1.7-4            remaCor_0.0.20          ggplot2_4.0.1          
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.6       jsonlite_2.0.0     compiler_4.5.2     EnvStats_3.1.0    
    ##  [5] Rcpp_1.1.1         stringr_1.6.0      jquerylib_0.1.4    systemfonts_1.3.1 
    ##  [9] scales_1.4.0       textshaping_1.0.4  yaml_2.3.12        fastmap_1.2.0     
    ## [13] lattice_0.22-7     R6_2.6.1           plyr_1.8.9         knitr_1.51        
    ## [17] rbibutils_2.4.1    desc_1.4.3         bslib_0.10.0       RColorBrewer_1.1-3
    ## [21] rlang_1.1.7        cachem_1.1.0       stringi_1.8.7      mathjaxr_2.0-0    
    ## [25] xfun_0.56          fs_1.6.6           sass_0.4.10        S7_0.2.1          
    ## [29] cli_3.6.5          pkgdown_2.2.0      withr_3.0.2        magrittr_2.0.4    
    ## [33] Rdpack_2.6.5       digest_0.6.39      grid_4.5.2         nlme_3.1-168      
    ## [37] lifecycle_1.0.5    vctrs_0.7.1        evaluate_1.0.5     glue_1.8.0        
    ## [41] farver_2.1.2       codetools_0.2-20   ragg_1.5.0         reshape2_1.4.5    
    ## [45] rmarkdown_2.30     tools_4.5.2        htmltools_0.5.9

### References

Han, Buhm, Dat Duong, Jae Hoon Sul, Paul IW de Bakker, Eleazar Eskin,
and Soumya Raychaudhuri. 2016. “A General Framework for Meta-Analyzing
Dependent Studies with Overlapping Subjects in Association Mapping.”
*Human Molecular Genetics* 25 (9): 1857–66.
<https://doi.org/10.1093/hmg/ddw049>.

Han, Buhm, and Eleazar Eskin. 2011. “Random-Effects Model Aimed at
Discovering Associations in Meta-Analysis of Genome-Wide Association
Studies.” *The American Journal of Human Genetics* 88 (5): 586–98.
<https://doi.org/10.1016/j.ajhg.2011.04.014>.

Lee, CH, Eleazar Eskin, and Buhm Han. 2017. “Increasing the Power of
Meta-Analysis of Genome-Wide Association Studies to Detect Heterogeneous
Effects.” *Bioinformatics* 33 (14): i379–88.
<https://doi.org/10.1093/bioinformatics/btx242>.

Lin, Dan-Yu, and Patrick F Sullivan. 2009. “Meta-Analysis of Genome-Wide
Association Studies with Overlapping Subjects.” *The American Journal of
Human Genetics* 85 (6): 862–72.
<https://doi.org/10.1016/j.ajhg.2009.11.001>.
