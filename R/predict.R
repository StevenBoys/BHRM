#' Function that evaluate the Richard function
#'
#' @param t - the time point to evaluate on the Richard curve
#' @param theta.1 - K : final epidemic size : real number
#' @param theta.2 - r : intrinsic growth rate : real number
#' @param theta.3 - tau - disease turning point : real number
#' @param xi: shape parameter (measuring deviation from symmetric, xi = 1 => SYMMETRIC)
#'
#' @return The evaluation of Richard function at time point t given the parameters.
#' @export
#'
#' @examples
#' t = 1
#' theta.1 = 1000; theta.2 = 0.1; theta.3 = 10; xi = 1
#' res = Richard_f(t, theta.1,theta.2, theta.3, xi)
Richard_f = function(t, theta.1,theta.2, theta.3, xi){
  res = theta.1/( 1 + xi*exp(-theta.2*(t - theta.3)))^(1/xi)
  return(res)
}

#' Function that make predictions based on the trained model
#'
#' @param model - a list which is the output of function BHRM_cov, BHRM, or BRM
#' @param Y - N-by-T time sereis training data for N countries for T days which are used to get the model
#' @param i - index that indicates which subject will be predicted
#' @param num_days - integer that indicates the number of days to predict and the default value is 14
#' @param conf - a numerical value between 0 and 1 which indicates the confidence level of the interval
#' @param seed - random seed set in the function
#'
#' @return A list of
#'        \item{prediction}{ - the predicted values}
#'        \item{upper}{ - 95% upper bound for the predicted values}
#'        \item{lower}{ - 95% lower bound for the predicted values}
#' @export
#'
#' @examples
#' data("time_series_data")
#' Y = time_series_data[, -c(1:2)]
#' seed.no = 1 ; burn = 20000 ; nmc = 20000 ; thin = 30; varrho = 0
#' pro.var.theta.2 = 0.0002 ; pro.var.theta.3 = 0.05; mu = 0 ; rho.sq = 1
#' res = BRM(Y = Y[1, ], seed.no = seed.no, burn = burn, nmc = nmc,
#' thin = thin, varrho = varrho, pro.var.theta.2 = pro.var.theta.2,
#' pro.var.theta.3 = pro.var.theta.3, mu = mu, rho.sq = rho.sq)
#' prediction = predict(model, Y, 1)
predict = function(model, Y, i, num_days = 14, conf = 0.95, seed = 5){
  index_show = list()
  set.seed(seed)
  # get the prediction mean
  t.values = c(1:(ncol(Y)+num_days))
  y.temp = matrix(0, nrow = ncol(model$thinned.theta.1.vec), ncol = length(t.values))
  for (t in 1:length(t.values)){
    for (s in 1:ncol(res$thinned.theta.1.vec)){
      mu = Richard_f(t = t.values[t],
             theta.1 = res$thinned.theta.1.vec[i, s],
             theta.2 = res$thinned.theta.2.vec[i, s],
             theta.3 = res$thinned.theta.3.vec[i, s],
             xi = res$thinned.xi.vec[i, s])
      sigma.sq = res$thinned.sigma.sq[s]
      y.temp[s,t] = rnorm(nsample, mu, sqrt(sigma.sq))
    }
  }
  pred = round(colMeans(y.temp))
  # get the prediction interval
  up = 1 - (1 - conf)/2
  low = (1 - conf)/2
  upper = apply(y.temp, 2, function(t)quantile(t, up))
  lower = apply(y.temp, 2, function(t)quantile(t, low))

  list(prediction = pred, upper = upper, lower = lower)
}


