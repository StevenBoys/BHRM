
flat_time_point_cal = function(theta.1,theta.2,theta.3,xi,gamma = 0.99){
  nomi = log((1/xi) * ( (1/(gamma))^{xi} -1 ) )
  denom = theta.2
  res = theta.3 - nomi/denom
  return(res)
}

#' Function that calculates the flat time point of a Richards curve given the parameters
#'
#' @param model - a list which is the output of function BHRM_cov, BHRM, or BRM
#' @param Y - N-by-T time sereis training data for N countries for T days which are used to get the model
#' @param i - index that indicates which subject will be extrapolated
#' @param gamma - a progression constant between 0 and 1 which determines when it can be called flat time point.
#' @param y_name - the name of the y axis. If it's not specified, the name for the y axis will be "Values".
#' @param title - the title of the plot. If it's not specified, the plot won't include title.
#' @param size_axis_text - the text size of y and x axis and the default value is 20
#' @param size_axis_title - the text size of y and x title and the default value is 25
#' @param size_plot_title - the text size of plot title and the default value is 35
#'
#' @return A list of
#'        \item{figure}{ - figure that shows the flat time points in the trajectory}
#'        \item{flat_time_points}{ - a vector that stores the flat time points}
#'        \item{epi_size}{ - the estimated final epidemic size}
#' @export
#'
#' @examples
#' data("time_series_data")
#' Y = time_series_data[, -c(1:2)]
#' seed.no = 1 ; burn = 20000 ; nmc = 20000 ; thin = 30; varrho = 0
#' pro.var.theta.2 = 0.0002 ; pro.var.theta.3 = 0.05; mu = 0 ; rho.sq = 1
#' num_days = 14; t.values = c(1:(ncol(Y) - num_days))
#' Y = Y[, c(1:(ncol(Y) - num_days))]
#' res = BHRM(Y = Y[1, ], t.values = t.values, seed.no = seed.no, burn = burn,
#' nmc = nmc, thin = thin, varrho = varrho, pro.var.theta.2 = pro.var.theta.2,
#' pro.var.theta.3 = pro.var.theta.3, mu = mu, rho.sq = rho.sq)
#' out = flat_time_point(res, Y, 1)
flat_time_point = function(model, Y, i, gamma = c(0.9999, 0.999, 0.99), y_name = NULL, title = NULL,
                           size_axis_text = 20, size_axis_title = 25, size_plot_title = 35){
  thinned.theta.1.vec_cov = model$thinned.theta.1.vec
  thinned.theta.2.vec_cov = model$thinned.theta.2.vec
  thinned.theta.3.vec_cov = model$thinned.theta.3.vec
  thinned.xi.vec = model$thinned.xi.vec

  # calculate the flat time points
  flat_time_points = rep(0, length(gamma))
  S = ncol(thinned.theta.1.vec_cov)
  for (e in 1:length(gamma)){
    temp = c()
    for (s in 1:S){
      temp[s] = flat_time_point_cal(theta.1 = thinned.theta.1.vec_cov[i,s],
                                theta.2 = thinned.theta.2.vec_cov[i,s],
                                theta.3 = thinned.theta.3.vec_cov[i,s],
                                xi = thinned.xi.vec[i,s],
                                gamma = gamma[e])
    }
    flat_time_points[e] = mean(temp)
  }

  # begin to make the plot
  max_flat_time = max(flat_time_points)
  extra = extrapolate(model, Y, i, num_days = max_flat_time + 10 - ncol(Y))

  library(ggplot2)
  library(scales)
  library(ggthemes)
  # set the relevant parameters
  mean = extra$mean; y = Y[i, ]; upper = extra$upper
  lower = extra$lower; num_days = max_flat_time + 10 - ncol(Y)

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
    geom_line(data = df_mean, mapping = aes(x = t.values, y = post_mean), col = "red", size = 1.2) +
    geom_point(data = df_obs, mapping = aes(x = t.values, y = obs), size = 2) +
    xlab("Dates") + ylab(y_name) +
    theme(
      axis.text=element_text(size=size_axis_text),
      axis.title=element_text(size=size_axis_title),
      plot.title = element_text(size = size_plot_title, face = "bold"))+
    theme_hc() + scale_colour_hc() +
    ylim(y_min, y_max) +
    scale_y_continuous(labels=comma) +
    geom_vline(xintercept = c(floor(flat_time_points)), linetype="dashed",
               color = "violet", size=1) +
    geom_vline(xintercept = c(length(y)), linetype="dotted",
               color = "grey", size=1)+
    geom_hline(yintercept = c(mean(thinned.theta.1.vec_cov[i,])), linetype="dashed",
               color = "blue", size=1)

  if(!is.null(upper[0])){
    g.res = g.res +
      geom_polygon(data = df_cov_polygon, mapping = aes(x = t.values, y = bounds), fill = "#87CEEB", alpha = 0.5)
  }
  if(!is.null(title)){
    g.res = g.res + ggtitle(title)
  }

  list(figure = g.res, flat_time_points = flat_time_points, epi_size = mean(thinned.theta.1.vec_cov[i,]))
}



