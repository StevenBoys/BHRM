#' Function that implement Bayesian Hierarchical Richard model without covariates
#'
#' @param Y - N-by-T time sereis data for N countries for T days
#' @param t.values - N list for time perids for N countries
#' @param seed.no - scalar for random seed
#' @param burn - no of burn
#' @param nmc - no of iterations after burn
#' @param thin - no of thining for the nmc
#' @param varrho - hyper-parameter for the diffuse prior error variances ~ IG(varrho,varrho)
#' @param pro.var.theta.2 - proposal variances for the MH algorithms to sample from theta.2.i
#' @param pro.var.theta.3 - proposal variances for the MH algorithms to sample from theta.3.i
#' @param mu - mu
#' @param rho.sq - rho.sq
#'
#' @return A list of
#'        \item{thinned.theta.1.vec}{ - thinned.theta.1.vec}
#'        \item{thinned.theta.2.vec}{ - thinned.theta.2.vec}
#'        \item{thinned.theta.3.vec}{ - thinned.theta.3.vec}
#'        \item{thinned.xi.vec}{ - thinned.xi.vec}
#'        \item{thinned.alpha.1}{ - thinned.alpha.1}
#'        \item{thinned.alpha.2}{ - thinned.alpha.2}
#'        \item{thinned.alpha.3}{ - thinned.alpha.3}
#'        \item{thinned.sigma.sq}{ - thinned.sigma.sq}
#'        \item{thinned.sigma.1.sq}{ - thinned.sigma.1.sq}
#'        \item{thinned.sigma.2.sq}{ - thinned.sigma.2.sq}
#'        \item{thinned.sigma.3.sq}{ - thinned.sigma.3.sq}
#'        \item{mu}{ - mu}
#'        \item{rho.sq}{ - rho.sq}
#' @export
#'
#' @examples
#' data("time_series_data")
#' Y = time_series_data[, -c(1:2)]
#' seed.no = 1 ; burn = 40000 ; nmc = 20000 ; thin = 30; varrho = 0
#' pro.var.theta.2 = 0.0002 ; pro.var.theta.3 = 0.05; mu = 0 ; rho.sq = 1
#' res_noncov = BHRM(Y = Y, seed.no = seed.no, burn = burn, nmc = nmc,
#' thin = thin, varrho = varrho, pro.var.theta.2 = pro.var.theta.2,
#' pro.var.theta.3 = pro.var.theta.3, mu = mu, rho.sq = rho.sq)
BHRM = function(Y,t.values,seed.no=1,burn=2000,nmc=2000,thin=10, varrho = 0, pro.var.theta.2 = 0.0002, pro.var.theta.3 = 0.05, mu = 0, rho.sq = 1){
  # Calling packages
  {
    library(mvtnorm)
    library(truncdist)
    library(nleqslv)
  }

  # About Data
  {
    N = nrow(Y)
    y = function(i){
      y = as.numeric(Y[i,c(t.values[[i]])])
      return(y)
    }
  }

  # MCMC setting
  {
    burn = burn
    nmc = nmc
    thin = thin
    S = burn + nmc
  }

  # Make rooms for parameters
  {
    # Model parameters & measurement error
    theta.1.vec = matrix(rep(0, N*S), nrow = N, ncol = S)
    theta.2.vec = matrix(rep(0, N*S), nrow = N, ncol = S)
    theta.3.vec = matrix(rep(0, N*S), nrow = N, ncol = S)
    xi.vec = matrix(rep(0, N*S), nrow = N, ncol = S)
    sigma.sq = rep(0,S)

    # Prior
    alpha.1 = rep(0,S)
    alpha.2 = rep(0,S)
    alpha.3 = rep(0,S)
    sigma.1.sq = rep(0,S)
    sigma.2.sq = rep(0,S)
    sigma.3.sq = rep(0,S)

  }

  # Decide initial values
  {
    y_temp = function(i){
      y = as.numeric(Y[i,1:max(t.values[[i]])])
      return(y)
    }
    for (i in 1:N){

      #temp.d = data.frame(y=y(i), t=c(1:T.days[i]))
      temp.d = data.frame(y=y_temp(i), t=c(1:max(t.values[[i]])))

      theta.1.vec[i,1] = max(temp.d$y)
      theta.3.vec[i,1] = which.max(diff(temp.d$y))

      xi.vec[i,1] = 0.1

      f = function(t, theta.1,theta.2, theta.3, xi){
        # (Original) Richard model
        # theta.1: K : final epidemic size : real number
        # theta.2: r : intrinsic growth rate : real number
        # theta.3: tau : disease turning point : real number
        # xi: shape parameter (measuring deviation from symmetric, xi = 1 => SYMMETRIC)
        res = theta.1/( 1 + xi*exp(-theta.2*(t - theta.3)))^(1/xi)
        return(res)
      }

      eps = 0.00000001
      temp.ft = function(x){
        nomi = log(f(t = theta.3.vec[i,1] + eps, theta.1 = theta.1.vec[i,1], theta.2 = x, theta.3 = theta.3.vec[i,1], xi = xi.vec[i,1])/
                     f(t = theta.3.vec[i,1], theta.1 = theta.1.vec[i,1], theta.2 = x, theta.3 = theta.3.vec[i,1], xi = xi.vec[i,1]))
        denom = eps
        # nomi/denom corresponds to the derivative of the logarithm of the Richard curve
        res = x/2 - nomi/denom
        return(res)
      }

      theta.2.vec[i,1] = nleqslv(x = 0.1, fn = temp.ft)$x

    }

    sigma.sq[1] = 1
    alpha.1[1] = 1
    alpha.2[1] = 1
    alpha.3[1] = 1
    sigma.1.sq[1] = 1
    sigma.2.sq[1] = 1
    sigma.3.sq[1] = 1

  }

  # Seed no of Gibbs sampler (for replication purpose)
  {
    set.seed(seed.no)
  }

  # Gibbs sampler
  for (s in 1:(S-1)){

    # Define necessary functions
    if (s == 1){

      # norm function
      norm_vec = function( x , y = rep(0,length(x)) ) {
        norm_vector = sqrt(sum((x-y)^2))

        if (norm_vector == Inf){
          norm_vector = 1e+308
        } else {
          norm_vector = norm_vector
        }
        return(norm_vector)
      }

      # One vector
      one.vec.N = rep(1,N)
      # Identity matrix
      I.N = diag(one.vec.N)

      # Richard model
      f = function(t, theta.1,theta.2, theta.3, xi){
        # (Original) Richard model
        # theta.1: K : final epidemic size : real number
        # theta.2: r : intrinsic growth rate : real number
        # theta.3: tau : disease turning point : real number
        # xi: shape parameter (measuring deviation from symmetric, xi = 1 => SYMMETRIC)
        res = theta.1/( 1 + xi*exp(-theta.2*(t - theta.3)))^(1/xi)
        return(res)
      }

      # T_i-dimensional vector for Richard model
      f.i.vec = function(i,theta.1.i,theta.2.i,theta.3.i, xi.i){
        res = matrix(data = f(c(t.values[[i]]), theta.1 = theta.1.i, theta.2 = theta.2.i, theta.3 = theta.3.i, xi = xi.i), ncol = 1)
        return(res)
      }

      # Ingredients for bulk updater for theta.1.vec
      h = function(t, theta.2, theta.3, xi){
        res = 1/( 1 + xi*exp(-theta.2*(t - theta.3)))^(1/xi)
        return(res)
      }

      h.i.vec = function(i,theta.2.i,theta.3.i,xi.i){
        res = matrix(data = h(c(t.values[[i]]), theta.2 = theta.2.i, theta.3 = theta.3.i, xi = xi.i), ncol = 1)
        return(res)
      }

      r = function(i, theta.2.i,theta.3.i, xi.i){
        res = t(y(i))%*%h.i.vec(i=i,theta.2.i=theta.2.i,theta.3.i=theta.3.i, xi.i=xi.i)
        return(res)
      }

      r.vec = function(theta.2.vec, theta.3.vec, xi.vec){
        temp = c()
        for (i in 1:N) {
          temp[i] = r(i = i, theta.2.i = theta.2.vec[i],theta.3.i = theta.3.vec[i],xi.i = xi.vec[i])
        }
        res = matrix(data = temp, nrow = N, ncol = 1)
        return(res)
      }

      H.mat = function(theta.2.vec, theta.3.vec, xi.vec){
        temp = c()
        for (i in 1:N) {
          temp[i] = norm_vec(x = h.i.vec(i = i, theta.2.i = theta.2.vec[i], theta.3.i = theta.3.vec[i], xi.i = xi.vec[i]))^2
        }
        res = diag(temp)
        return(res)

      }

      # Ingredients to calculate criterion function for MH to update theta.2.vec
      f.2.i = function(i, theta.2.i, old.theta.2.i,theta.1.i,theta.3.i, xi.i){
        temp.1 = f.i.vec(i = i, theta.2.i = theta.2.i, theta.1.i = theta.1.i, theta.3.i = theta.3.i, xi.i = xi.i)
        temp.2 = f.i.vec(i = i, theta.2.i = old.theta.2.i, theta.1.i = theta.1.i, theta.3.i = theta.3.i, xi.i = xi.i)
        res = temp.1 - temp.2
        return(res)
      }
      g.2.i = function(i, theta.2.i, old.theta.2.i,theta.1.i,theta.3.i, xi.i){
        temp.1 = norm_vec(x = f.i.vec(i = i, theta.2.i = theta.2.i, theta.1.i = theta.1.i, theta.3.i = theta.3.i, xi.i = xi.i))^2
        temp.2 = norm_vec(x = f.i.vec(i = i, theta.2.i = old.theta.2.i, theta.1.i = theta.1.i, theta.3.i = theta.3.i, xi.i = xi.i))^2
        res = temp.1 - temp.2
        return(res)
      }
      log.Likelihood.ft.2 = function(i,theta.2.i, old.theta.2.i, theta.1.i, theta.3.i, xi.i, sigma.sq){
        temp.1 = sum(t(y(i))%*%f.2.i(i = i,theta.2.i=theta.2.i,old.theta.2.i=old.theta.2.i,theta.1.i=theta.1.i,theta.3.i=theta.3.i, xi.i = xi.i))
        temp.2 = g.2.i(i = i,theta.2.i=theta.2.i,old.theta.2.i=old.theta.2.i,theta.1.i=theta.1.i,theta.3.i=theta.3.i, xi.i = xi.i)
        res = (1/sigma.sq)*(temp.1 - (1/2)*temp.2)
        return(res)
      }
      alpha.ft.2 = function(i = i, theta.2.i, old.theta.2.i, theta.1.i, theta.3.i, xi.i, sigma.sq, alpha.2, sigma.2.sq, beta.2.vec){
        temp = exp(log.Likelihood.ft.2(i = i, theta.2.i = theta.2.i, old.theta.2.i = old.theta.2.i, theta.1.i = theta.1.i, theta.3.i = theta.3.i, xi.i = xi.i, sigma.sq = sigma.sq) +
                     (1/sigma.2.sq)*((alpha.2)*(theta.2.i-old.theta.2.i) - (1/2)*(theta.2.i^2 - old.theta.2.i^2)) )
        res = max(min(temp,1),1e-10)
        return(res)
      }


      # Ingredients to calculate criterion function for MH to update theta.3.vec
      f.3.i = function(i, theta.3.i, old.theta.3.i,theta.1.i,theta.2.i, xi.i){
        temp.1 = f.i.vec(i = i, theta.3.i = theta.3.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i, xi.i = xi.i)
        temp.2 = f.i.vec(i = i, theta.3.i = old.theta.3.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i, xi.i = xi.i)
        res = temp.1 - temp.2
        return(res)
      }
      g.3.i = function(i, theta.3.i, old.theta.3.i,theta.1.i,theta.2.i, xi.i){
        temp.1 = norm_vec(x = f.i.vec(i = i, theta.3.i = theta.3.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i ,xi.i = xi.i))^2
        temp.2 = norm_vec(x = f.i.vec(i = i, theta.3.i = old.theta.3.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i ,xi.i = xi.i))^2
        res = temp.1 - temp.2
        return(res)
      }
      log.Likelihood.ft.3 = function(i, theta.3.i, old.theta.3.i, theta.1.i, theta.2.i, xi.i, sigma.sq){
        temp.1 = sum(t(y(i))%*%f.3.i(i = i,theta.3.i=theta.3.i,old.theta.3.i=old.theta.3.i,theta.1.i=theta.1.i,theta.2.i=theta.2.i, xi.i = xi.i))
        temp.2 = g.3.i(i = i,theta.3.i=theta.3.i,old.theta.3.i=old.theta.3.i,theta.1.i=theta.1.i,theta.2.i=theta.2.i, xi.i = xi.i)
        res = (1/sigma.sq)*(temp.1 - (1/2)*temp.2)
        return(res)
      }
      alpha.ft.3 = function(i = i, theta.3.i, old.theta.3.i, theta.1.i, theta.2.i, xi.i, sigma.sq, alpha.3, sigma.3.sq, beta.3.vec){
        temp = exp(log.Likelihood.ft.3(i = i, theta.3.i = theta.3.i, old.theta.3.i = old.theta.3.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i, xi.i = xi.i, sigma.sq = sigma.sq) +
                     (1/sigma.3.sq)*((alpha.3)*(theta.3.i-old.theta.3.i) - (1/2)*(theta.3.i^2 - old.theta.3.i^2)) )
        res = max(min(temp,1),1e-10)
        return(res)
      }

      # Ingredients for ESS to sample from eta : eta = log(xi)
      f.eta.i = function(i, eta.i, old.eta.i,theta.1.i,theta.2.i, theta.3.i){
        xi.i = exp(eta.i)
        old.xi.i = exp(old.eta.i)
        temp.1 = f.i.vec(i = i, xi.i = xi.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i, theta.3.i = theta.3.i)
        temp.2 = f.i.vec(i = i, xi.i = old.xi.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i, theta.3.i = theta.3.i)
        res = temp.1 - temp.2
        return(res)
      }
      g.eta.i = function(i, eta.i, old.eta.i,theta.1.i,theta.2.i, theta.3.i){
        xi.i = exp(eta.i)
        old.xi.i = exp(old.eta.i)
        temp.1 = norm_vec(x = f.i.vec(i = i, xi.i = xi.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i, theta.3.i = theta.3.i))^2
        temp.2 = norm_vec(x = f.i.vec(i = i, xi.i = old.xi.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i, theta.3.i = theta.3.i))^2
        res = temp.1 - temp.2
        return(res)
      }
      log.Likelihood.ft.eta = function(i, eta.i, old.eta.i,theta.1.i,theta.2.i, theta.3.i, sigma.sq){
        temp.1 = sum(t(y(i))%*% f.eta.i(i = i, eta.i = eta.i, old.eta.i = old.eta.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i, theta.3.i = theta.3.i))
        temp.2 = g.eta.i(i = i, eta.i = eta.i, old.eta.i = old.eta.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i, theta.3.i = theta.3.i)
        res = (1/sigma.sq)*(temp.1 - (1/2)*temp.2)
        return(res)
      }
      alpha.ft.eta = function(i, eta.i, old.eta.i,theta.1.i,theta.2.i, theta.3.i, sigma.sq){
        temp = exp(log.Likelihood.ft.eta(i = i,eta.i = eta.i, old.eta.i = old.eta.i, theta.1.i = theta.1.i, theta.2.i = theta.2.i, theta.3.i = theta.3.i, sigma.sq = sigma.sq))
        res = max(min(temp,1),1e-10)
        return(res)
      }
    }

    # Step 1 : updating theta.1.vec, theta.2.vec, and theta.3.vec
    {
      # 1-i : theta.1.vec (Closed-form update)
      {
        Sigma.theta.1.vec = solve((1/sigma.sq[s])*H.mat(theta.2.vec = theta.2.vec[,s], theta.3.vec = theta.3.vec[,s], xi.vec = xi.vec[,s]) +
                                    (1/sigma.1.sq[s])*I.N)
        mu.theta.1.vec = Sigma.theta.1.vec%*%((1/sigma.sq[s])*r.vec(theta.2.vec = theta.2.vec[,s],theta.3.vec = theta.3.vec[,s], xi.vec = xi.vec[,s]) +
                                                (1/sigma.1.sq[s])*(one.vec.N*alpha.1[s]))
        theta.1.vec[,(s+1)] = mvtnorm::rmvnorm(n=1, mean = mu.theta.1.vec, sigma = Sigma.theta.1.vec)
      }
      # 1-ii : theta.2.vec (MH)
      {
        for (i in 1:N){
          # Step 1: proposal step
          new.theta.2.i = rnorm(n = 1, mean = theta.2.vec[i,s], sd = sqrt(pro.var.theta.2))
          # Step 2: calculate acceptance ratio
          temp.ratio = alpha.ft.2(i = i, theta.2.i = new.theta.2.i,
                                  old.theta.2.i = theta.2.vec[i,s],
                                  theta.1.i = theta.1.vec[i,(s+1)],
                                  theta.3.i = theta.3.vec[i,s],
                                  xi.i = xi.vec[i,s],
                                  sigma.sq = sigma.sq[s],
                                  alpha.2 = alpha.2[s],
                                  sigma.2.sq = sigma.2.sq[s],
                                  beta.2.vec = beta.2.vec[,s])

          # Step 3: Criteria
          u = runif(n = 1, min = 0, max = 1)
          if (u <= temp.ratio){
            theta.2.vec[i,(s+1)] = new.theta.2.i
          }else{
            theta.2.vec[i,(s+1)] = theta.2.vec[i,s]
          }
        }
      }
      # 1-iii :theta.3.vec (MH)
      {

        for(i in 1:N){
          # Step 1: proposal step
          new.theta.3.i = rnorm(n = 1, mean = theta.3.vec[i,s], sd = sqrt(pro.var.theta.3))
          # Step 2: calculate acceptance ratio
          temp.ratio = alpha.ft.3(i = i,theta.3.i = new.theta.3.i,
                                  old.theta.3.i = theta.3.vec[i,s],
                                  theta.1.i = theta.1.vec[i,(s+1)],
                                  theta.2.i = theta.2.vec[i,(s+1)],
                                  xi.i = xi.vec[i,s],
                                  sigma.sq = sigma.sq[s],
                                  alpha.3 = alpha.3[s],
                                  sigma.3.sq = sigma.3.sq[s],
                                  beta.3.vec = beta.3.vec[,s])
          # Step 3: Criteria
          u = runif(n = 1, min = 0, max = 1)
          if (u <= temp.ratio){
            theta.3.vec[i,(s+1)] = new.theta.3.i
          }else{
            theta.3.vec[i,(s+1)] = theta.3.vec[i,s]
          }
        }
      }
      # 1-iv : xi.vec (ESS)
      {
        for(i in 1:N){
          # step 1: variable change
          old.eta.i = log(xi.vec[i,s])
          # step 2: ESS: eta.i ~ \pi(eta.i|-)
          {

            # old.eta.i   : previous realization for the eta.i
            # eta.i.star  : proposed eta.i in ESS
            # new.eta.i   : accepted eta.i

            # ESS:1 Choose an ellipse:
            nu = rnorm(n = 1, mean = mu, sd = sqrt(rho.sq))
            # ESS:2 Define a criterion function:
            # Defined above
            # ESS:3 Choose a threshold and fix:
            u = runif(n = 1, min = 0, max = 1)
            # ESS:4 Draw an initial proposal:
            theta = runif(n = 1, min = -pi, max = pi)
            eta.i.star = (old.eta.i-mu)*cos(theta) + (nu -mu)*sin(theta) + mu

            # ESS:5 ESS core step
            if (u < alpha.ft.eta(i = i, eta.i = eta.i.star, old.eta.i = old.eta.i, theta.1.i = theta.1.vec[i,(s+1)],theta.2.i = theta.2.vec[i,(s+1)],theta.3.i = theta.3.vec[i,(s+1)],sigma.sq = sigma.sq[s])){
              new.eta.i = eta.i.star
            }else{
              # Define the bracket:
              theta.min = -pi ; theta.max = pi
              while (u >= alpha.ft.eta(i = i, eta.i = eta.i.star, old.eta.i = old.eta.i, theta.1.i = theta.1.vec[i,(s+1)],theta.2.i = theta.2.vec[i,(s+1)],theta.3.i = theta.3.vec[i,(s+1)],sigma.sq = sigma.sq[s])){
                # Shrink the bracket and try a new point:
                if (theta > 0) { theta.max =  theta} else { theta.min = theta}
                theta = runif(n = 1, min = theta.min , max = theta.max)
                eta.i.star = (old.eta.i-mu)*cos(theta) + (nu - mu)*sin(theta) + mu
              }
              new.eta.i = eta.i.star
            }

          }
          # step 3: transform-back
          xi.vec[i,(s+1)] = exp(new.eta.i)
        }
      }
    }

    # Step 2 : updating sigma.sq
    {
      temp.1 = length(unlist(t.values))/2 + varrho
      temp.2 = c()
      for(i in 1:N){
        temp.2[i] = norm_vec(x = y(i), y = f.i.vec(i = i, theta.1.i = theta.1.vec[i,(s+1)],theta.2.i = theta.2.vec[i,(s+1)],theta.3.i = theta.3.vec[i,(s+1)], xi.i = xi.vec[i,(s+1)]))^2
      }
      temp.3 = (1/2)*sum(temp.2) + varrho
      sigma.sq[s+1] = 1/rgamma(n = 1, shape = temp.1, rate = temp.3)
    }

    # Step 3 : updating alpha.1, alpha.2, and alpha.3
    {
      # 3-i : alpha.1
      {
        temp.1 = (1/N)*t(one.vec.N)%*%(theta.1.vec[,(s+1)])
        temp.2 = sqrt(sigma.1.sq[s]/N)
        alpha.1[s+1] = rnorm(n = 1, mean = temp.1, sd = temp.2)
      }
      # 3-ii : alpha.2
      {
        temp.1 = (1/N)*t(one.vec.N)%*%(theta.2.vec[,(s+1)])
        temp.2 = sqrt(sigma.2.sq[s]/N)
        alpha.2[s+1] = rnorm(n = 1, mean = temp.1, sd = temp.2)
      }
      # 3-iii : alpha.3
      {
        temp.1 = (1/N)*t(one.vec.N)%*%(theta.3.vec[,(s+1)])
        temp.2 = sqrt(sigma.3.sq[s]/N)
        alpha.3[s+1] = rnorm(n = 1, mean = temp.1, sd = temp.2)
      }
    }

    # Step 4 : updating sigma.1.sq, sigma.2.sq, and sigma.3.sq
    {
      # 4-i : sigma.1.sq
      {
        temp.1 = (N)/2 + varrho
        temp.2 = (1/2)*norm_vec(x = theta.1.vec[,(s+1)], y = one.vec.N*alpha.1[s+1] )^2 + varrho
        sigma.1.sq[s+1] = 1/rgamma(n = 1, shape = temp.1, rate = temp.2)
      }
      # 4-ii : sigma.2.sq
      {
        temp.1 = (N)/2 + varrho
        temp.2 = (1/2)*norm_vec(x = theta.2.vec[,(s+1)], y = one.vec.N*alpha.2[s+1])^2 + varrho
        sigma.2.sq[s+1] = 1/rgamma(n = 1, shape = temp.1, rate = temp.2)
      }
      # 4-iii : sigma.3.sq
      {
        temp.1 = (N)/2 + varrho
        temp.2 = (1/2)*norm_vec(x = theta.3.vec[,(s+1)], y = one.vec.N*alpha.3[s+1])^2 + varrho
        sigma.3.sq[s+1] = 1/rgamma(n = 1, shape = temp.1, rate = temp.2)
      }
    }

    remainder = s%%1000
    if(remainder == 0){
      print(paste("no of iterations is", s))
      mc.index = seq(from = burn + 1, to = burn + nmc, by = thin)
      print(paste("no of curves with negative slope is", sum(rowMeans(theta.2.vec[,mc.index]) < 0)))
    }

  }

  # Print-out results
  {
    mc.index = seq(from = burn + 1, to = burn + nmc, by = thin)
    thinned.theta.1.vec = theta.1.vec[,mc.index]
    thinned.theta.2.vec = theta.2.vec[,mc.index]
    thinned.theta.3.vec = theta.3.vec[,mc.index]
    thinned.xi.vec = xi.vec[,mc.index]
    thinned.alpha.1 = alpha.1[mc.index]
    thinned.alpha.2 = alpha.2[mc.index]
    thinned.alpha.3 = alpha.3[mc.index]
    thinned.sigma.sq = sigma.sq[mc.index]
    thinned.sigma.1.sq = sigma.1.sq[mc.index]
    thinned.sigma.2.sq = sigma.2.sq[mc.index]
    thinned.sigma.3.sq = sigma.3.sq[mc.index]

    res = list(thinned.theta.1.vec = thinned.theta.1.vec,
               thinned.theta.2.vec = thinned.theta.2.vec,
               thinned.theta.3.vec = thinned.theta.3.vec,
               thinned.xi.vec = thinned.xi.vec,
               thinned.alpha.1 = thinned.alpha.1,
               thinned.alpha.2 = thinned.alpha.2,
               thinned.alpha.3 = thinned.alpha.3,
               thinned.sigma.sq = thinned.sigma.sq,
               thinned.sigma.1.sq = thinned.sigma.1.sq,
               thinned.sigma.2.sq = thinned.sigma.2.sq,
               thinned.sigma.3.sq = thinned.sigma.3.sq,
               mu = mu,
               rho.sq = rho.sq)

    return(res)
  }
}
