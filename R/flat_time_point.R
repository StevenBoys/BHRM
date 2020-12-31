#' Function that calculates the flat time point of a Richards curve given the parameters
#'
#' @param theta.1 - parameter theta1 of Richards curve.
#' @param theta.2 - parameter theta2 of Richards curve.
#' @param theta.3 - parameter theta3 of Richards curve.
#' @param xi - parameter xi of Richards curve.
#' @param gamma - a progression constant between 0 and 1 which determines when it can be called flat time point.
#'
#' @return The flat time point of a Richards curve given the parameters.
#' @export
#'
#' @examples
#' data("time_series_data")
#' Y = time_series_data[, -c(1:2)]
#' seed.no = 1 ; burn = 20000 ; nmc = 20000 ; thin = 30; varrho = 0
#' pro.var.theta.2 = 0.0002 ; pro.var.theta.3 = 0.05; mu = 0 ; rho.sq = 1
#' num_days = 14; t.values = c(1:(ncol(Y) - num_days))
#' Y = Y[, c(1:(ncol(Y) - num_days))]
#' res = BRM(Y = Y[1, ], t.values = t.values, seed.no = seed.no, burn = burn,
#' nmc = nmc, thin = thin, varrho = varrho, pro.var.theta.2 = pro.var.theta.2,
#' pro.var.theta.3 = pro.var.theta.3, mu = mu, rho.sq = rho.sq)
#' theta.1 = mean(res$thinned.theta.1.vec)
#' theta.2 = mean(res$thinned.theta.2.vec)
#' theta.3 = mean(res$thinned.theta.3.vec)
#' xi = mean(res$thinned.xi.vec)
#' flat_time_point(theta.1,theta.2,theta.3,xi)
flat_time_point = function(theta.1,theta.2,theta.3,xi,gamma = 0.99){
  nomi = log((1/xi) * ( (1/(gamma))^{xi} -1 ) )
  denom = theta.2
  res = theta.3 - nomi/denom
  return(res)
}
