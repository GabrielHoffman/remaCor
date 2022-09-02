## Random effects meta-analysis for correlated test statistics

<p align="center">
<img src=man/figures/image.png width="400">
</p>

Standard approaches to meta-analysis assumes that effect sizes are statistically independent. Here we provide methods for fixed and random effects meta-analysis when the correlation between effect sizes are known.

## Usage
```r
# Run fixed effects meta-analysis, accounting for correlation 
LS( beta, stders, Sigma)

# Run random effects meta-analysis, accounting for correlation 
RE2C( beta, stders, Sigma)
```


## Install from GitHub
```r
devtools::install_github("GabrielHoffman/remaCor")
```
