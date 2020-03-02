# remaCor: Random effects meta-analysis for correlated test statistics

<p align="center">
<img src=https://users.hpc.mssm.edu/~hoffmg01/software/remaCor/logo.png width="400">
</p>

Standard approaches to meta-analysis assumes that effect sizes are statistically independent. Here we provide methods for fixed and random effects meta-analysis when the correlation between effect sizes are known.

## Install from GitHub

```r
library(devtools)

install_github("GabrielHoffman/remaCor")
```

## Usage
Basic usage:
```r
# Run random effects meta-analysis, accounting for correlation 
RE2C( beta, stders, Sigma)
```

## [Vignette](https://users.hpc.mssm.edu/~hoffmg01/software/remaCor/remaCor.html)
## [Manual](https://users.hpc.mssm.edu/~hoffmg01/software/remaCor/remaCor-manual.pdf)

This is a developmental version.


