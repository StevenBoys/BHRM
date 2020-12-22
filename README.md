# BHRM Package

Bowen Lei, Se Yoon Lee, Bani Mallick

# BHRM: Bayesian Hierarchical Richard Model

This R package is designed to help users to train the Bayesian Hierarchical Richard Model (BHRM) which is based on the Richards curve:

<img src="https://latex.codecogs.com/gif.latex?f(t&space;;&space;\theta_1,&space;\theta_2,&space;\theta_3,&space;\xi)=\theta_1&space;\cdot[&space;1&space;&plus;&space;\xi&space;\cdot&space;\exp&space;\{-\theta_2&space;\cdot&space;(&space;t&space;-&space;\theta_3)&space;\}&space;]^{-1/\xi}" title="f(t ; \theta_1, \theta_2, \theta_3, \xi)=\theta_1 \cdot[ 1 + \xi \cdot \exp \{-\theta_2 \cdot ( t - \theta_3) \} ]^{-1/\xi}" />

The model can uncover a hidden pattern from growth curves. At the same time, users can choose covariate version and identify important predictors that largely affect on the shape of the curve f in terms of the three curve parameters.

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

Take the COVID-19 data as an example. The time_series_data include infection trajectories for several global countries and the design_matrix include potential covariates. 
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
```

We choose the Bayesian hierarchical Richard model with covariates to analyse the data.
```
# set the hyperparameters
seed.no = 1 ; burn = 40000 ; nmc = 20000 ; thin = 30; varrho = 0
pro.var.theta.2 = 0.0002 ; pro.var.theta.3 = 0.05; mu = 0 ; rho.sq = 1
# run the model
res_cov = BHRM_cov(Y = Y, X = X, seed.no = seed.no, burn = burn, nmc = nmc,  
                   thin = thin, varrho = varrho, pro.var.theta.2 = pro.var.theta.2, 
                   pro.var.theta.3 = pro.var.theta.3, mu = mu, rho.sq = rho.sq)  
```

We can use `var_sele` function to check the variable selection results.
```
# check the important factors for theta1
var_selection = var_sele(beta.vec = res_cov$thinned.theta.1.vec, names = names(X))
# check the names of the top covariates selected
var_selection$names_sele
# plot the figure for 95% credible interval of each covariates
var_selection$figure
```




## References

[1] S. Y. Lee, B. Lei, and B. K. Mallick, Estimation of COVID-19 spread curves integrating global data and borrowing information, PLoS ONE 15 (7), (2020).
