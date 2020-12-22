# BHRM Package

Bowen Lei, Se Yoon Lee, Bani Mallick

# BHRM: Bayesian Hierarchical Richard Model

This R package is designed to help users to train the Bayesian Hierarchical Richard Model (BHRM) which is based on the Richards curve:

<img src="https://latex.codecogs.com/gif.latex?f(t&space;;&space;\theta_1,&space;\theta_2,&space;\theta_3,&space;\xi)=\theta_1&space;\cdot[&space;1&space;&plus;&space;\xi&space;\cdot&space;\exp&space;\{-\theta_2&space;\cdot&space;(&space;t&space;-&space;\theta_3)&space;\}&space;]^{-1/\xi}" title="f(t ; \theta_1, \theta_2, \theta_3, \xi)=\theta_1 \cdot[ 1 + \xi \cdot \exp \{-\theta_2 \cdot ( t - \theta_3) \} ]^{-1/\xi}" />

The model can uncover a hidden pattern from growth curves. At the same time, users can choose covariate version and identify important predictors that largely affect on the shape of the curve f in terms of the three curve parameters.

The details of the model is as below:
![](https://github.com/StevenBoys/BHRM/blob/main/Image/BHRM_formula.png?raw=true =200x100)

## Installation

```
require(devtools)
devtools::install_github("StevenBoys/BHRM", build_vignettes = T)
```

## Usage

```
library(BHRM)
```


## References

[1] S. Y. Lee, B. Lei, and B. K. Mallick, Estimation of COVID-19 spread curves integrating global data and borrowing information, PLoS ONE 15 (7), (2020).
