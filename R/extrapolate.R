#' Function that evaluate the Richards function
#'
#' @param t - the time point to evaluate on the Richards curve
#' @param theta.1 - K : final epidemic size : real number
#' @param theta.2 - r : intrinsic growth rate : real number
#' @param theta.3 - tau - disease turning point : real number
#' @param xi: shape parameter (measuring deviation from symmetric, xi = 1 => SYMMETRIC)
#'
#' @return The evaluation of Richards function at time point t given the parameters.
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

#' Function that make extrapolation based on the trained model
#'
#' @param model - a list which is the output of function BHRM_cov, BHRM, or BRM
#' @param Y - N-by-T time sereis training data for N countries for T days which are used to get the model
#' @param i - index that indicates which subject will be extrapolated
#' @param num_days - integer that indicates the number of days to extrapolate and the default value is 14
#' @param conf - a numerical value between 0 and 1 which indicates the confidence level of the interval
#' @param seed - random seed set in the function
#'
#' @return A list of
#'        \item{mean}{ - the extrapolated values}
#'        \item{upper}{ - 95% upper bound for the extrapolated values}
#'        \item{lower}{ - 95% lower bound for the extrapolated values}
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
#' extra_list = extrapolate(model, Y, 1)
extrapolate = function(model, Y, i, num_days = 14, conf = 0.95, seed = 5){
  index_show = list()
  set.seed(seed)
  # get the extrapolated mean
  t.values = c(1:(ncol(Y)+num_days))
  y.temp = matrix(0, nrow = ncol(model$thinned.theta.1.vec), ncol = length(t.values))
  for (t in 1:length(t.values)){
    for (s in 1:ncol(model$thinned.theta.1.vec)){
      mu = Richard_f(t = t.values[t],
             theta.1 = model$thinned.theta.1.vec[i, s],
             theta.2 = model$thinned.theta.2.vec[i, s],
             theta.3 = model$thinned.theta.3.vec[i, s],
             xi = model$thinned.xi.vec[i, s])
      y.temp[s,t] = mu
    }
  }
  mean = round(colMeans(y.temp))
  # get the confidence interval
  up = 1 - (1 - conf)/2
  low = (1 - conf)/2
  upper = apply(y.temp, 2, function(t)quantile(t, up))
  lower = apply(y.temp, 2, function(t)quantile(t, low))

  list(mean = mean, upper = upper, lower = lower)
}

#' Function that make a plot and can see a comparison between the observations and extrapolations
#'
#' @param mean - a vector of the mean of extrapolated values.
#' @param y - a vector of the observations
#' @param upper - the upper bound of the extrapolated values. If it's not specified, the plot won't draw the interval.
#' @param lower - the lower bound of the extrapolated values. If it's not specified, the plot won't draw the interval.
#' @param num_days - integer that indicates the number of days to extrapolate and the default value is 14.
#' @param title - the title of the plot. If it's not specified, the plot won't include title.
#' @param y_name - the name of the y axis. If it's not specified, the name for the y axis will be "Values".
#' @param size_axis_text - the text size of y and x axis and the default value is 20
#' @param size_axis_title - the text size of y and x title and the default value is 25
#' @param size_plot_title - the text size of plot title and the default value is 35
#'
#' @return A plot that can see a comparison between the observations and extrapolations.
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
#' extra_list = extrapolate(model, Y, 1)
#' plot_RM(extra_list$mean, Y[1, ], extra_list$upper, extra_list$lower)
plot_RM = function(mean, y, upper = NULL, lower = NULL, num_days = 14, title = NULL, y_name = NULL,
                   size_axis_text = 20, size_axis_title = 25, size_plot_title = 35){
  library(ggplot2)
  library(scales)
  library(ggthemes)

  # define data frame for the mean values
  t.values = c(1:length(mean))
  df_mean = data.frame(cbind(t.values, mean))
  names(df_mean) <- c("t.values", "post_mean")

  # define the data frame for the observations
  y = as.vector(as.matrix(y))
  df_obs = as.data.frame(cbind(c(1:length(y)), y))
  names(df_obs) = c("t.values","obs")

  if(!is.null(upper[0])){
    # define the data frame for the pointwise confidence interval
    t.values = c((length(y)+1):length(mean))
    df_cov = data.frame(cbind(t.values, lower[t.values], upper[t.values]))
    names(df_cov) <- c("t.values", "lower_bound", "upper_bound")
    df_cov_polygon = as.data.frame(cbind(c(df_cov$t.values,rev(df_cov$t.values)),
                                         c(df_cov$lower_bound,rev(df_cov$upper_bound)))  )
    names(df_cov_polygon) = c("t.values","bounds")

    # calculate the range of the y axis
    y_min = min(c(y, mean, lower[t.values], upper[t.values]))
    y_max = max(c(y, mean, lower[t.values], upper[t.values]))
  }else{
    # calculate the range of the y axis
    y_min = min(c(y, mean))
    y_max = max(c(y, mean))
  }

  # set the name of y axis
  if(is.null(y_name)){
    y_name = "Values"
  }

  # make the plot based on the defined data frame
  g.res = ggplot() +
    geom_line(data = df_mean, mapping = aes(x = t.values, y = post_mean), col = "blue", size = 1.5) +
    geom_point(data = df_obs, mapping = aes(x = t.values, y = obs),size = 2.5, col = "#000000") +
    xlab("Days") + ylab(y_name) +
    theme(
      axis.text=element_text(size=size_axis_text),
      axis.title=element_text(size=size_axis_title),
      plot.title = element_text(size = size_plot_title, face = "bold"))+
    theme_hc() + scale_colour_hc() +
    ylim(y_min, y_max) +
    scale_y_continuous(labels=comma)

  if(!is.null(upper[0])){
    g.res = g.res +
      geom_polygon(data = df_cov_polygon, mapping = aes(x = t.values, y = bounds), fill = "#87CEEB", alpha = 0.5)
  }
  if(!is.null(title)){
    g.res = g.res + ggtitle(title)
  }

  return(g.res)
}
