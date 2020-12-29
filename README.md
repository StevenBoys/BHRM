# BHRM: An R package implementing Bayesian Hierarchical Richard Model (BHRM) to Extrapolate Infection Trajectories and Identify Risk Factors for the COVID-19 Outbreak

## Contents
* [Overview](#overview)
* [Installation](#installation)
* [Examples](#examples)
* [References](#References)

## Overview
## Installation
## Examples
## References

# BHRM: Bayesian Hierarchical Richard Model

This R package is designed to help users to train the Bayesian Hierarchical Richard Model (BHRM) which is based on the Richards curve:

<img src="https://latex.codecogs.com/gif.latex?f(t&space;;&space;\theta_1,&space;\theta_2,&space;\theta_3,&space;\xi)=\theta_1&space;\cdot[&space;1&space;&plus;&space;\xi&space;\cdot&space;\exp&space;\{-\theta_2&space;\cdot&space;(&space;t&space;-&space;\theta_3)&space;\}&space;]^{-1/\xi}" title="f(t ; \theta_1, \theta_2, \theta_3, \xi)=\theta_1 \cdot[ 1 + \xi \cdot \exp \{-\theta_2 \cdot ( t - \theta_3) \} ]^{-1/\xi}" />

![](https://github.com/StevenBoys/BHRM/blob/main/Image/Richard_f.png?raw=true)

The model can uncover a hidden pattern from growth curves. At the same time, users can choose covariate version and identify important predictors that largely affect on the shape of the curve in terms of the three curve parameters.

The details of the model is as below:
![](https://github.com/StevenBoys/BHRM/blob/main/Image/BHRM_formula.png?raw=true)

## Installation

```
require(devtools)
devtools::install_github("StevenBoys/BHRM", build_vignettes = T)
```

## Usage

```
library(BHRM)
```

Take the COVID-19 data as an example. 

![](https://github.com/StevenBoys/BHRM/blob/main/Image/infect_COVID-19.png?raw=true)

The time_series_data include infection trajectories for several global countries and the design_matrix include potential covariates. 
```
# load the data
data("design_matrix")
data("time_series_data")
Y = time_series_data[, -c(1:2)]; X = design_matrix[, -c(1:2)]
# standardize the design matrix
norm_vec = function( x , y = rep(0,length(x)) ) {
  norm_vector = sqrt(sum((x-y)^2))
  
  if (norm_vector == Inf){
    norm_vector = 1e+308
  } else {
    norm_vector = norm_vector
  } 
  return(norm_vector)
}
for (j in 1:ncol(X)){
  X[,j] = X[,j] - mean(X[,j])
  X[,j] = X[,j]/norm_vec(x = X[,j], y = rep(0,nrow(X)))
}
X = as.matrix(X)
```

We choose the Bayesian hierarchical Richard model with covariates to analyse the data.
```
# set the hyperparameters
seed.no = 1 ; burn = 20000 ; nmc = 20000 ; thin = 30; varrho = 0
pro.var.theta.2 = 0.0002 ; pro.var.theta.3 = 0.05; mu = 0 ; rho.sq = 1
t.values = list(); num_days = 14
for(i in 1:nrow(Y)){
  t.values[[i]] = c(1:(ncol(Y) - num_days))
}
Y = Y[, c(1:(ncol(Y) - num_days))]
# run the model
res_cov = BHRM_cov(Y = Y, X = X, t.values = t.values, seed.no = seed.no, burn = burn,   
                   nmc = nmc, thin = thin, varrho = varrho, pro.var.theta.2 = pro.var.theta.2, 
                   pro.var.theta.3 = pro.var.theta.3, mu = mu, rho.sq = rho.sq)  
```

We can use `var_sele` function to check the variable selection results.
```
# check the important factors for theta1
var_selection = var_sele(beta.vec = res_cov$thinned.theta.1.vec)

# check the names of the top covariates selected
var_selection$id_sele
# [1] 10 30  7 15 19 23 28 27  2 24

# plot the figure for 95% credible interval of each covariates
var_selection$figure
```

![](https://github.com/StevenBoys/BHRM/blob/main/Image/var_sele.png?raw=true)

We can also use `predict` to make predictions and use `plot_RM` to make a plot and compare the real trajectory and predictions.
```
# make predictions
predict_list = predict(res_cov, Y, 1)
# make a plot to see the performance of the predictions
plot_RM(predict_list$prediction, Y[1, ])
```

![](https://github.com/StevenBoys/BHRM/blob/main/Image/prediction.png?raw=true)

The flat time point is defined as the the time point whereat only (1−gamma)theta1 cases can maximally take place to reach the final epidemic size.

![](https://github.com/StevenBoys/BHRM/blob/main/Image/flat_time_point.png?raw=true)

We can also compute the flat point of the estimated Richard curve using function `flat_time_point`.
```
theta.1 = mean(res_cov$thinned.theta.1.vec[1, ])
theta.2 = mean(res_cov$thinned.theta.2.vec[1, ])
theta.3 = mean(res_cov$thinned.theta.3.vec[1, ])
xi = mean(res_cov$thinned.xi.vec[1, ])
flat_time_point(theta.1,theta.2,theta.3,xi)
# 156.0234
```

## References

[1] S. Y. Lee, B. Lei, and B. K. Mallick, Estimation of COVID-19 spread curves integrating global data and borrowing information, PLoS ONE 15 (7), (2020).
