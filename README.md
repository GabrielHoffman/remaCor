# remaCor: Random effects meta-analysis for correlated test statistics

<p align="center">
<img src=https://hoffmg01.u.hpc.mssm.edu/software/remaCor/logo.png width="400">
</p>

Standard approaches to meta-analysis assumes that effect sizes are statistically independent. Here we provide methods for fixed and random effects meta-analysis when the correlation between effect sizes are known.

## Install from GitHub

```
library(devtools)

install_github("GabrielHoffman/remaCor")
```

## Usage
Basic usage:
```{r }
# Run random effects meta-analysis, accounting for correlation 
RE2C( beta, stders, Sigma)
```

## [See Vignette](https://hoffmg01.u.hpc.mssm.edu/software/remaCor/remaCor.html)

This is a developmental version.


