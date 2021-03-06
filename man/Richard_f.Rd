% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extrapolate.R
\name{Richard_f}
\alias{Richard_f}
\title{Function that evaluate the Richards function}
\usage{
Richard_f(t, theta.1, theta.2, theta.3, xi)
}
\arguments{
\item{t}{- the time point to evaluate on the Richards curve}

\item{theta.1}{- K : final epidemic size : real number}

\item{theta.2}{- r : intrinsic growth rate : real number}

\item{theta.3}{- tau - disease turning point : real number}

\item{xi:}{shape parameter (measuring deviation from symmetric, xi = 1 => SYMMETRIC)}
}
\value{
The evaluation of Richards function at time point t given the parameters.
}
\description{
Function that evaluate the Richards function
}
\examples{
t = 1
theta.1 = 1000; theta.2 = 0.1; theta.3 = 10; xi = 1
res = Richard_f(t, theta.1,theta.2, theta.3, xi)
}
