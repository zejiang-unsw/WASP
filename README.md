# WASP

An open-source wavelet tool for improving prediction accuracy for natural system models

## Requirements
<pre>
Dependencies:
  waveslim, stats, tidyr, ggplot2, sp

Suggest:
    zoo, readr,
    cowplot, SPEI, FNN, 
    NPRED, synthesis, fitdistrplus
</pre>

## Installation

You can install the package via devtools from [GitHub](https://github.com/) with:

``` r
devtools::install_github("zejiang-unsw/WASP", dependencies = TRUE)
```

or via CRAN with: 

``` r
install.packages("WASP")
```

## Citation
Jiang, Z., Sharma, A., & Johnson, F. (2021). Variable transformations in the spectral domain – Implications for hydrologic forecasting. Journal of Hydrology, 126816. [doi](https://doi.org/10.1016/J.JHYDROL.2021.126816)

Jiang, Z., Rashid, M. M., Johnson, F., & Sharma, A. (2020). A wavelet-based tool to modulate variance in predictors: an application to predicting drought anomalies. Environmental Modelling & Software, 135, 104907. [doi](https://doi.org/10.1016/j.envsoft.2020.104907)

Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3), e2019WR026962. [doi](https://doi.org/10.1029/2019WR026962)
