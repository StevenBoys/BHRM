# BHRM Package

Bowen Lei, Se Yoon Lee, Bani Mallick

# BHRM: Bayesian Hierarchical Richard Model

This R package is designed to help users to train the Bayesian Hierarchical Richard Model (BHRM) which is based on the Richards curve:
$$
f(t ; \theta_1, \theta_2, \theta_3, \xi) &= \theta_1 \cdot [ 1 + \xi \cdot \exp \{-\theta_2 \cdot ( t - \theta_3) \}   ]^{-1/\xi}
$$

The model can uncover a hidden pattern from growth curves. At the same time, users can choose covariate version and identify important predictors that largely affect on the shape of the curve $f(t ; \theta_1, \theta_2, \theta_3, \xi)$ in terms of the three curve parameters.

The details of the model is as below:
$$
y_{it}&=f(t ; \theta_{1i}, \theta_{2i}, \theta_{3i}, \xi_i) + \epsilon_{it}, \quad \,\epsilon_{it}\sim \mathcal{N}(0,\sigma^2),\,\,\,\,\,\,\,\,\,\,\, (i=1, \cdots, N, \, t = 1, \cdots, T),\\
\label{eq:linear_regression}
\theta_{li} &= \alpha_l + \textbf{x}_i^{\top} \bm{\beta}_l + \varepsilon_{li}, \quad\quad \varepsilon_{li}\sim \mathcal{N}(0,\sigma_l^2),\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,\,(i=1, \cdots, N, \, l = 1, 2,3),\\
\label{eq:Horseshoe}
\beta_{lj}&|\lambda_{lj},\tau_{lj}, \sigma_l^2  \sim \mathcal{N}(0, \sigma_l^2 \tau_l^2 \lambda_{lj}^2)
,\quad
\lambda_{lj},\tau_{lj} \sim \mathcal{C}^{+}(0,1),\,\,\,\,\,\,\,\,
(l = 1, 2,3,\, j = 1,\cdots, p),\\
\label{eq:xi_log_normal}
\xi_i &\sim \log\ \mathcal{N}(0,1),\,\,\,\,\,\,\,\,\,\,\,\,\,\quad\quad\quad\quad\quad\quad\quad\quad\quad\quad\quad\quad\quad  (i=1, \cdots, N),\\
\label{eq:improper_priors}
\alpha_l &\sim \pi(\alpha) \propto 1,\,\quad\quad
\sigma^2 ,\sigma_l^2 \sim \pi(\sigma^2) \propto 1/\sigma^2,
\quad\quad\quad\quad (l = 1, 2,3),
$$

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
