###############################################################################
##
## Functions used to simulate data
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: August 2024
##
###############################################################################

# This script contains internal functions used to generate the exported 'data'.
# Additional details on data generation can also be found in the script
# 'DATASET.R' in the 'data-raw' folder


#' Simulate from a simple ARDL(1,1)
#'
#' This function encloses the DGP for the exported data. If 'linear = TRUE',
#' it generates a panel of N units from the model: Yit = rho*Yt-1 + beta*Xt-1 + eps_it,
#' where eps_it ~ iid(0, sigma). If 'linear = FALSE' the function generates a panel dataset
#' of N units from the following non-linear model: Y_it = sin(rho*Yt-1 + beta*Xt-1) + eps_it.
#'
#'
#' @param seed              Numeric, random seed
#' @param X                 List of covariates
#' @param beta              Numeric vector of coefficients for 'X'
#' @param ar1_coef          Numeric, autoregressive coefficient
#' @param N                 Number of units in the panel
#' @param sigma             Standard deviation of the error term
#' @param impact            Numeric vector of fictional additive causal effects
#' @param impact_constant   Logical, if \code{TRUE} the values in \code{impact} are additive constant effects
#'                          otherwise they are first scaled by the standard deviation of the outcome
#' @param linear            Logical, defaults to \code{TRUE}. Set to \code{FALSE} for nonlinear specification
#' @param post_per          Number of post-intervention periods after the intervention date which, by default,
#'                          is included in the post-int period. So, for a single post-intervention, set \code{post_per = 0}
#'
#' @return                  A list with two components: \code{dat} (the generated dataset) and \code{true_ate}
#'                          (the true ATE, as results from the simulated effects added in the data)
#' @noRd
#'
.sim_ardl <- function(seed, X, beta, ar1_coef, N, sigma, impact, impact_constant, post_per, linear = TRUE){

  # Param checks
  if(!is.numeric(seed) | seed < 0 ) stop (" 'seed' must be a number greater than 0") # include test in testthat
  if(!is.list(X)) stop(" 'X' must be a list") # testthat it works with both
  if(!is.numeric(beta)) stop (" 'beta' must be numeric") # test with testthat
  if(!is.numeric(ar1_coef)) stop (" 'ar1_coef' must be numeric ")
  if(N %% 1 != 0 | N < 0) stop (" 'N' must be an integer") # test with testthat
  if(!is.numeric(impact) | length(impact) != post_per +1) stop ("'impact' is non-numeric or it doesn't match 'post_per' ")
  if(!all(sapply(list(impact_constant, linear), is.logical))) stop ("supplied non-logical values to logical parameters")

  # Settings
  set.seed(seed)
  t <- NROW(X[[1]])

  # Generating eps_it.
  eps_it <- lapply(1:N, function(x)(rnorm(n = t, mean = 0, sd = sigma)))

  # Generating Yit recursively (using 'y_recursive') for each i in 1:N
  Yt <- mapply(FUN = .y_recursive, X = X, eps = eps_it, MoreArgs = list(beta, ar1_coef), linear = linear, SIMPLIFY = F)

  # Adding fictional intervention
  Y0 <- sapply(Yt, function(x) x[(t - post_per):t])
  # Y0 <- sapply(Yt, function(x) x[t]) # old


  if(impact_constant){

    Yt <- lapply(Yt, function(x){x[(t - post_per):t] <- x[(t - post_per):t] + impact; x})
    # Yt <- lapply(Yt, function(x){x[t] <- x[t] + impact; x}) # old

  }  else {

    Yt <- lapply(Yt, function(x){x[(t - post_per):t] <- x[(t - post_per):t] + impact*sd(x); x})
    # Yt <- lapply(Yt, function(x){x[t] <- x[t] + impact*sd(x); x}) # old

  }

  # Dataset in long format
  dat <- data.frame(ID = rep(1:N, each = t),
                    Y = unlist(Yt),
                    # Xlag1 = do.call(rbind, lapply(X, function(x)(apply(x,2, .true_lag)))),
                    X = do.call(rbind, X),
                    year = rep(seq(2020-t+1,2020), times = N))
  colnames(dat) <- gsub(colnames(dat), pattern = ".x", replacement = ".")

  # Removing initial NAs
  # ind <- which(is.na(data$Xlag1.x1))
  # dat <- dat[-ind, ]

  # Returning results (dataset + true ATE)
  return(list(dat = dat,
              true_ate = rowMeans(matrix(sapply(Yt, function(x) x[(t-post_per):t]) - Y0, nrow = post_per + 1))))
  # true_ate = mean(sapply(Yt, function(x) x[t]) - Y0))) # old

}


#' Auxiliary function, used in '.sim_ardl()'
#'
#' This function generates Yt for a single unit, starting from its covariates 'X'
#' the given autoregressive coefficient and the error terms 'eps' (eps can include the
#' individual heterogeneity)
#'
#' @param beta       Numeric vector of coefficients for 'X'
#' @param X          Matrix of covairates
#' @param eps        Numeric, simulated error terms
#' @param ar1_coef   Numeric, autoregressive coefficient
#' @param linear     Logical, defaults to \code{TRUE}.
#'
#' @return           Numerical vector, generated time series from the ARDL (linear or non-linear) model
#' @noRd
#'
.y_recursive <- function(beta, X, eps, ar1_coef, linear){

  # Settings
  Yt <- c()
  t <- NROW(X)

  if(!linear){

    # Initializing
    Yt[1] <- eps[1]

    # Looping
    for(i in 2:t){

      Yt[i] <- sin(as.vector(X[i-1,]%*%beta) + ar1_coef*Yt[i-1]) + eps[i]

    }

  } else {

    # Initializing, Y0 = 0 and Y1 = beta*X1 + eps_it
    Yt[1] = eps[1]

    # Looping
    for(i in 2:t){

      Yt[i] <- as.vector(X[i-1,]%*%beta) + ar1_coef*Yt[i-1] + eps[i]

    }
  }

  return(Yt)

}
