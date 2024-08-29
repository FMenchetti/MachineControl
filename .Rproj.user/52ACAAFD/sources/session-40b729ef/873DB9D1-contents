###############################################################################
## 
## GENERATING PANEL DATA from different DGPs
##  
## Author:      Cerqua A., Letta M., Menchetti F.
## Last update: Feb 2023
## Functions:   sim_ardl , sim_np, sim_gen_np
##
##############################################################################

library(MASS)

##############################################################################
## MAIN SIMULATION FUNCTION
##############################################################################

#' Inferring causal effects from simulated panel datasets under different DGP
#'
#' @param gen_mod Character, takes values in c("ardl", "np", "nl")
#' @param niter Number of datasets to generate
#' @param X Matrix or data.frame of exogenous regressors
#' @param beta Numeric vector containing regressors' coefficients 
#' @param N Numeric, number of units in the panel
#' @param sigma Standard deviation of the error term 
#' @param ar1_coef The AR1 coefficient (needed when y.lag > 0)
#' @param inf_type Character, type of inference to be performed. Possible choices are 
#'                 'classic', 'block', 'bc classic', 'bc block' and 'bca'.
#' @param impact Numeric, assumed additive effect. Notice that internally, 'impact'
#'        multiplies the standard deviation of each unit time-series, so ideally should be
#'        either 1 or 2, meaning that we are assuming an intervention increasing the response
#'        of 1 or 2 standard deviations (we can test multiple impacts of course). 
#' @param linear Logic, whether the DGP should be linear or not.        
#' @param impact_constant Logic, whether to include a constant additive effect or an additive effect 
#'        proportional to the sample stardard deviation.
#' @param post_per Number of post intervention periods       
#' @param y.lag Optional, numeric indicating how many lags of Y to include
#' @param setseed Seed to replicate simulation results. Defaults to 1.
#' @param pcv_block Number of pre-intervention times to block for the panel cross validation. 
#' @param nboot Number of bootstrap replications, defaults to 1000.
#'
#' @return A data.frame with the following information: best method, estimated ATE and its variance
#'         Upper and lower confidence interval bound, true effect, interval length, coverage, bias and power


PanelSIM <- function(gen_mod, niter, X, beta, N, sigma, ar1_coef, inf_type, impact, linear, impact_constant, post_per, y.lag, setseed = 1, pcv_block = 1, nboot = 1000){
 
  # Generating data
  if(gen_mod == "ardl"){
    
    data_list <- lapply(1:niter, FUN = function(x)(sim_ardl(seed = x*setseed, beta = beta, X = X, N = N, sigma = sigma,
                                                            ar1_coef = ar1_coef, impact = impact, impact_constant = impact_constant,  
                                                            ylag = FALSE, linear = linear, post_per = post_per)))
  } 
  
  if(gen_mod == "np"){
    
    stop("Not yet optimized")
  }
  
  # Estimation and inference with the MLCM function 
  sim <- lapply(data_list, FUN = function(x)(MLCM(data = x$dat, y = "Y", id = "ID", timevar = "year", int_date = 2020 - post_per, inf_type = inf_type, pcv_block = pcv_block, nboot = nboot, y.lag = y.lag)))

  # Tabulating results 
  postimes <- names(sim[[1]]$ate)
  dims <- c(length(sim),1,length(postimes))
  tab <- aperm(array(mapply(x = rep(sim, times = length(unique(postimes))), y = rep(postimes, each = length(sim)), 
                     FUN = function(x,y){c(x$ate[y], x$var.ate[y], x$conf.ate[2,y], x$conf.ate[1,y])}), 
               dim = c(4, length(sim), length(postimes)), dimnames = list(c("ATE", "VAR_ATE", "Upper", "Lower"), 1:length(sim), postimes)), 
         perm = c(2, 1, 3))
  true_ate <- array(t(sapply(data_list, function(x) x$true_ate)), dims)
  tab <- abind(tab, True_effect = true_ate, along = 2)
  int_length <- array(tab[, "Upper", ] - tab[, "Lower",], dims)
  cover <- array(ifelse(tab[,"True_effect",] >= tab[,"Lower",] & tab[, "True_effect",] <=  tab[,"Upper",], 1, 0), dims)
  bias <- array(abs(tab[, "True_effect",] - tab[, "ATE",]), dims)
  power <- array(ifelse(0 >= tab[, "Lower",] & 0 <= tab[, "Upper",], 0, 1), dims)
  tab <- abind(tab, Int_length = int_length, Coverage = cover, Bias = bias, Power = power, along = 2)
  
  # Returning results
  return(list(tab = tab, Best_method = sapply(sim, FUN = function(x)(x$best_method$method))))
  
}

###########################################################
## Function to simulate from an ARDL(1,1)
###########################################################

# If 'linear = TRUE', the function generates a panel of N units from the model
# Yit = rho*Yt-1 + beta*Xt-1 + eps_it, eps_it ~ iid(0, sigma). If 'linear = FALSE'
# the function generates a panel dataset of N units from the following non-linear model

sim_ardl <- function(seed, X, beta, ar1_coef, N, sigma, impact, impact_constant, ylag, linear, post_per){
  
  # Settings
  set.seed(seed)
  t <- NROW(X[[1]])
  
  # Do we want different covariates'values for different units?
  # if(varying){
    
  #  ran <- lapply(1:N, FUN = mvrnorm, n = NROW(X) , mu = rep(1, ncol(X)), Sigma = 0.1*diag(1, nrow = ncol(X), ncol = ncol(X)))
  #  X <- lapply(ran, function(x)(x + X))
    
  # }
  
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
                    Xlag1 = do.call(rbind, lapply(X, function(x)(apply(x,2, .true_lag)))),
                    # X = do.call(rbind, X),
                    year = rep(seq(2020-t+1,2020), times = N))
  # if(ylag){
  #   
  #   # Lagged dependent variable
  #   dat$Ylag1 <- unlist(lapply(Yt, FUN = .true_lag))
  #   ind <- which(is.na(dat$Ylag1))
  #   
  # } else {
  #   
  #   ind <- which(is.na(dat$Xlag1.x1))
  #   
  # }
  
  # Removing initial NAs 
  # dat <- dat[-ind, ]
  
  # Returning results (dataset + true ATE)
  return(list(dat = dat, 
              true_ate = rowMeans(matrix(sapply(Yt, function(x) x[(t-post_per):t]) - Y0, nrow = post_per + 1))))
              # true_ate = mean(sapply(Yt, function(x) x[t]) - Y0))) # old
  
}


###########################################################
## Function to simulate from a non-parametric model with 2 covariates
###########################################################

# DA MODIFICARE PER IL POST INTERVENTO
# The function simulates a panel of N units from the model
# Yit = 10*I(x1<q1 & x2<q1) + 15*I(q1<x1<me & q1<x2<me) + 20*I(x1>me & x2>me) + eps_it, eps_it ~ iid(0, sigma)

sim_np <- function(seed, X, N, impact, varying = T){
  
  # Settings
  set.seed(seed)
  X_list <- list()
  Y_list <- list()
  
  for(i in 1:N){
    
    # Do we want different covariates'values for different units?
    if(varying){
      X <- apply(X, 2, function(x)(x + rnorm(NROW(x), 0, 1)))
    }
    
    # Quantiles
    Q <- apply(X, 2, quantile, probs = c(0.25,0.5))
    
    # Generating Y and adding the fictional intervention 
    Yt <- 10*(X[,1] < Q[1,1] & X[,2] < Q[1,2]) + 15*(X[,1] >= Q[1,1] & X[,1] <= Q[2,1] & X[,2] >= Q[1,2] & X[,2] <= Q[2,2]) +
      20*(X[,1]> Q[2,1] & X[,2] > Q[2,2])
    Yt[t] <- Yt[t]*impact
    
    # Saving in a list
    Y_list[[i]] <- Yt
    X_list[[i]] <- X
  }
  
  # Returning results in long format
  res <- data.frame(ID = rep(1:N, each = t),
                    Y = unlist(Y_list),
                    X = do.call(rbind, X_list))
  return(res)
} 

####### Auxiliary functions, used in 'sim_ardl' 

.y_recursive <- function(beta, X, eps, ar1_coef, linear){
 
  ## This function generates Yt for one unit, starting from its covariates 'X'
  ## the given autoregressive coef and the error terms 'eps' (eps can include the
  ## individual heterogeneity)
  
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

# -----------------------------------------------------------------------------
.true_lag <- function(x, lag = 1){
  
  # x is a vector (not necessarily a ts object). This function shift x 
  # backward of the given number of 'lag' and returns a vector of the same length
  # of x
  x_lag <- c(rep(NA, times = lag), head(x, n = length(x)-lag))
  return(x_lag)
}