# BHRM: An R package implementing Bayesian Hierarchical Richards Model to Estimate Infection Trajectories and Identify Risk Factors for the COVID-19 Outbreak

![Figure 1: Extrapolated infectory for grand average over 40 countries](https://github.com/StevenBoys/BHRM/blob/main/Image/Global_average2.png)

## Contents
* [Overview](#overview)
* [Installation](#installation)
* [BHRM](#bhrm)
* [Examples](#examples)
* [References](#References)

## Overview
An R package `BHRM` of the paper titled  **["Estimation of COVID-19 spread curves integrating global data and borrowing information", PLOS ONE, (2020)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0236860)** are available here. This is a joint project of Ph.D. students, Bowen Lei (bowenlei@stat.tamu.edu) and [Se Yoon Lee](https://sites.google.com/view/seyoonlee) (seyoonlee@stat.tamu.edu), and a University Distinguished Professor [Bani K. Mallick](https://www.stat.tamu.edu/~bmallick/) (bmallick@stat.tamu.edu) at Texas A&M University. R package `BHRM` contains relevant R codes to implement Bayesian Hierarchical Richards Model (BHRM) applied to the COVID-19 dataset obtained from multiple countries. 

The sources of the datasets are: 
1. [Center for Systems Science and Engineering at Johns Hopkins University](https://github.com/CSSEGISandData/COVID-19)
2. [World Bank](https://data.worldbank.org/)
3. [World Health Organization](https://apps.who.int/gho/data/node.main)
4. [National Oceanic and Atmospheric Administration](https://www.noaa.gov/)

## Installation

```
require(devtools)
devtools::install_github("StevenBoys/BHRM")
library(BHRM)
```

## BHRM
Richards growth curve has been widely used to describe epidemiology for real-time prediction of outbreak of diseases. We propose a Bayesian hierarchical model based on the Richards curve (BHRM) to accommodate the global COVID-19 data. We aim to uncover a hidden pattern from the infection trajectory for each country and then extrapolate the curve. At the same time, we want to identify important predictors that largely affect on the shape the curve. The details of the hierarchy of the model is shown in the figure below.

<div align=center><img src="https://github.com/StevenBoys/BHRM/blob/main/Image/BHRM_formula.png?raw=true" width="700" height="500" alt=" "/></div>

## Examples
We take the `time_series_data` and `design_matrix` in this package as an example. `time_series_data` include infection growth curve for 30 global countries and `design_matrix` include 20 potential predictors.

![](https://github.com/StevenBoys/BHRM/blob/main/Image/infect_COVID-19.png?raw=true)

Firstly, we load the data.
```
library(BHRM)
# load the data
data("design_matrix")
data("time_series_data")
Y = time_series_data[, -c(1:2)]; X = design_matrix[, -c(1:2)]
```

Then, we choose the [`BHRM_cov`](https://github.com/StevenBoys/BHRM/blob/main/R/BHRM_cov.R) which refers to Bayesian hierarchical Richards model with covariates to analyse the data.
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

We can use [`var_sele`](https://github.com/StevenBoys/BHRM/blob/main/R/var_sele.R) function to check the variable selection results.
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

We can also use [`extrapolate`](https://github.com/StevenBoys/BHRM/blob/main/R/extrapolate.R) to make extrapolations and use [`plot_RM`](https://github.com/StevenBoys/BHRM/blob/main/R/extrapolate.R) to make a plot and compare the real trajectory and extrapolated values.
```
# make extrapolations
extra_list = extrapolate(res_cov, Y, 1)
# make a plot to see the performance of the extrapolations
plot_RM(extra_list$mean, Y[1, ])
```

![](https://github.com/StevenBoys/BHRM/blob/main/Image/prediction.png?raw=true)

We can also compute the flat point of the estimated Richards curve using function [`flat_time_point`](https://github.com/StevenBoys/BHRM/blob/main/R/flat_time_point.R).
```
theta.1 = mean(res_cov$thinned.theta.1.vec[1, ])
theta.2 = mean(res_cov$thinned.theta.2.vec[1, ])
theta.3 = mean(res_cov$thinned.theta.3.vec[1, ])
xi = mean(res_cov$thinned.xi.vec[1, ])
flat_time_point(theta.1,theta.2,theta.3,xi)
# 156.0234
```

## References

[1] [Se Yoon Lee, Bowen Lei, and Bani K. Mallick. (2020) “Estimation of COVID19 spread curves integrating global data and borrowing information,” PLOS ONE](https://journals.plos.org/plosone/article/authors?id=10.1371/journal.pone.0236860)

[2] [Davidian, M., and Giltinan, D. M. (1995). Nonlinear models for repeated measurement data (Vol. 62). CRC press.](https://books.google.com/books?hl=en&lr=&id=0eSIBPAL4qsC&oi=fnd&pg=IA7&dq=nonlinear+mixed+effect+model+giltnan&ots=9frDPH3F4J&sig=L5Wz91waGu447OdyYHQ8Vp5ckQc#v=onepage&q=nonlinear%20mixed%20effect%20model%20giltnan&f=false)

