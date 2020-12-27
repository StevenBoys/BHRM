#' Function that do variable selection based on the mcmc samples of the coefficients
#'
#' @param beta.vec - p-by-S matrix for S mcmc samples of p potential factors
#' @param max.rank - The number of top covariates selected based on the posterior mean
#'
#' @return A list of
#'        \item{figure}{ - 95% credible interval of each covariates}
#'        \item{id_sele}{ - the column id for the top covariates selected}
#'        \item{names_sele}{ - if names are provided, names_sele are the names of the top covariates; otherwise names_sel are null}
#' @export
#'
#' @examples
#' data("design_matrix")
#' data("time_series_data")
#' Y = time_series_data[, -c(1:2)]; X = design_matrix[, -c(1:2)]
#' X = scale(X)
#' seed.no = 1 ; burn = 40000 ; nmc = 20000 ; thin = 30; varrho = 0
#' pro.var.theta.2 = 0.0002 ; pro.var.theta.3 = 0.05; mu = 0 ; rho.sq = 1
#' res_cov = BHRM_cov(Y = Y, X = X, seed.no = seed.no, burn = burn, nmc = nmc,
#' thin = thin, varrho = varrho, pro.var.theta.2 = pro.var.theta.2,
#' pro.var.theta.3 = pro.var.theta.3, mu = mu, rho.sq = rho.sq)
#' var_selection = var_sele(beta.vec = res_cov$thinned.theta.1.vec, names = names(X))
var_sele = function(beta.vec, max.rank = 10, names = NULL){
  library("ggplot2")
  # Input should be posterior realizations of p dimensional vector of beta
  # Each column of the vector should be each iteration
  p = nrow(beta.vec)
  # y = X beta + epsilon
  cred.prob = function(x) quantile(x = x, probs = c(0.05, 0.95))
  res1 = apply(X = beta.vec, MARGIN = 1, FUN = cred.prob )
  res2 = rowMeans(x = beta.vec)

  # get the top max.rank covariates
  max.rank = min(max.rank, nrow(X))
  gene.rank = c(1:p)[order(abs(rowMeans(beta.vec)),decreasing = TRUE)]
  temp1 <- c(1:max.rank)[rowMeans(beta.vec)[gene.rank[1:max.rank]] > 0]
  temp2 <- c(1:max.rank)[rowMeans(beta.vec)[gene.rank[1:max.rank]] < 0]
  is_top <- rep(0, p)
  is_top[gene.rank[temp1]] <- 1
  is_top[gene.rank[temp2]] <- -1

  df = data.frame( post.mean =  res2,
                   post.lower.bound = res1[1,],
                   post.upper.bound = res1[2,],

                   covariates.index=c(1:p),
                   is.top = is_top)
  df$is.top <- as.factor(df$is.top)

  get_ylim <- function(data, times = 1){
    u <- c(); l <- c()
    for(i in 1:nrow(data)){
      u <- c(u, quantile(data[i, ], 0.975))
      l <- c(l, quantile(data[i, ], 0.025))
    }
    c(max(u)*times, min(l)*times)
  }
  y_lim <- get_ylim(beta.vec)

  figure = ggplot(df, aes(factor(covariates.index),
                          y = post.mean,
                          ymin = post.lower.bound,
                          ymax = post.upper.bound,
                          col = is.top))+
    geom_pointrange() +
    #xlab("covariate index j") +
    xlab("") +
    ylab(expression(paste("95% credible interval of ", beta[j]))) +
    geom_hline(yintercept=0, linetype = "dashed", color = "red")+
    scale_color_manual(values = c('1' = '#1E90FF', '-1' = '#FF4500', '0' = '#191970')) +
    theme(text = element_text(size=20),
          axis.ticks.y = element_blank()) +
    ylim(y_lim)

  id_sele = gene.rank[1:max.rank]
  if(!is.null(names)){
    names_sele = names[id_sele]
  }else{
    names_sele = NULL
  }

  list(figure = figure, id_sele = id_sele, names_sele = names_sele)

}
