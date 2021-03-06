% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flat_time_point.R
\name{flat_time_point}
\alias{flat_time_point}
\title{Function that calculates the flat time point of a Richards curve given the parameters}
\usage{
flat_time_point(
  model,
  Y,
  i,
  gamma = c(0.9999, 0.999, 0.99),
  y_name = NULL,
  title = NULL,
  size_axis_text = 20,
  size_axis_title = 25,
  size_plot_title = 35
)
}
\arguments{
\item{model}{- a list which is the output of function BHRM_cov, BHRM, or BRM}

\item{Y}{- N-by-T time sereis training data for N countries for T days which are used to get the model}

\item{i}{- index that indicates which subject will be extrapolated}

\item{gamma}{- a progression constant between 0 and 1 which determines when it can be called flat time point.}

\item{y_name}{- the name of the y axis. If it's not specified, the name for the y axis will be "Values".}

\item{title}{- the title of the plot. If it's not specified, the plot won't include title.}

\item{size_axis_text}{- the text size of y and x axis and the default value is 20}

\item{size_axis_title}{- the text size of y and x title and the default value is 25}

\item{size_plot_title}{- the text size of plot title and the default value is 35}
}
\value{
A list of
       \item{figure}{ - figure that shows the flat time points in the trajectory}
       \item{flat_time_points}{ - a vector that stores the flat time points}
       \item{epi_size}{ - the estimated final epidemic size}
}
\description{
Function that calculates the flat time point of a Richards curve given the parameters
}
\examples{
data("time_series_data")
Y = time_series_data[, -c(1:2)]
seed.no = 1 ; burn = 20000 ; nmc = 20000 ; thin = 30; varrho = 0
pro.var.theta.2 = 0.0002 ; pro.var.theta.3 = 0.05; mu = 0 ; rho.sq = 1
num_days = 14; t.values = c(1:(ncol(Y) - num_days))
Y = Y[, c(1:(ncol(Y) - num_days))]
res = BHRM(Y = Y[1, ], t.values = t.values, seed.no = seed.no, burn = burn,
nmc = nmc, thin = thin, varrho = varrho, pro.var.theta.2 = pro.var.theta.2,
pro.var.theta.3 = pro.var.theta.3, mu = mu, rho.sq = rho.sq)
out = flat_time_point(res, Y, 1)
}
