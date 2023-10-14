###############################################################################
##
## Bootstrap inference
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: July 2023
##
###############################################################################

#' Bootstrap inference for ATE
#'
#' Internal function, used within the MLCM routine for the estimation of ATE
#' standard error and 95% confidence interval.
#'
#' @param data A 'PanelMLCM' object from a previous call to \code{as.PanelMLCM}.
#' @param ind  int_date    The date of the intervention, treatment, policy introduction or shock. It must be contained in 'timevar'.
#' @param bestt Object of class \code{train}, the best-performing ML method as selected
#'              by panel cross validation.
#' @param type Character, type of inference to be performed. Possible choices are 'classic', 'block', 'bc classic', 'bc block', 'bca'.
#' @param nboot Number of bootstrap replications.
#' @param alpha Confidence interval level to report for the ATE. Defaulting to 0.05 for a two sided
#'              95\% confidence interval
#' @param metric Character, the performance metric that should be used to select the optimal model.
#'               Possible choices are either \code{"RMSE"} (the default) or \code{"Rsquared"}.
#' @param ate Numeric, the estimated ATE in the sample.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{type}: the inference type that has been performed
#'   \item \code{conf.ate}: bootstrap confidence interval at the 95% level
#'   \item \code{var.ate}: estimated variance for ATE
#'   \item \code{ate.lower}: lower confidence interval bound
#'   \item \code{ate.upper}: upper confidence interval bound
#' }
#' @noRd
#' @import caret
#' @importFrom bcaboot bcajack2
#' @importFrom stats quantile
#' @importFrom stats qnorm
#' @importFrom stats var
#' @importFrom stats pnorm
#' @importFrom stats rnorm
#' @importFrom stats sd


boot_ate <- function(data, int_date, bestt, type, nboot, alpha, metric, y.lag, ate = NULL){

  ### Param checks
  if(!any(class(data) %in% "PanelMLCM")) stop("Invalid class in the PanelCrossValidation function, something is wrong with as.PanelMLCM")
  if(!any(type %in% c("classic", "block", "bc classic", "bc block", "bca"))) stop("Inference type not allowed, check the documentation")
  if(nboot < 1 | all(!class(nboot) %in% c("numeric", "integer")) | nboot%%1 != 0) stop("nboot must be an integer greater than 1")
  if(alpha < 0 | alpha > 1) stop("Invalid confidence interval level, alpha must be positive and less than 1")
  if(!is.numeric(ate)) stop ("ATE must be numeric")


  ### Setting pre/post intervention periods
  pre <- which(data[, "Time"] < int_date)
  post <- which(data[, "Time"] >= int_date)

  ### Classic bootstrap
  if(type %in% c("classic", "bc classic", "bca")){

    # Step 1. Sampling the indices (time-id pairs)
    ii<- matrix(sample(pre, size = nboot*length(pre), replace = T), nrow = nboot, ncol = length(pre))

    # Step 2. Estimating individual effects and ATE for each bootstrap iteration
    ate_boot <- apply(ii, 1, function(i){ate_est(data = data[c(i, post),], int_date = int_date, best = bestt, metric = metric, y.lag = y.lag, ran.err = TRUE)$ate})


  }


  ### Block bootstrap
  if(type %in% c("block", "bc block")){

    # Step 1. Sampling the units
    ids <- unique(data$ID)
    ind1 <- sapply(1:nboot, function(boot){set.seed(boot); sample(ids, size = length(ids), replace = TRUE)}, simplify = FALSE)

    # Step 2. Indices pre/post interventions corresponding to the sampled units
    ii0 <- sapply(ind1, FUN = function(x){

      unlist(sapply(x, function(y)(intersect(which(data[, "ID"] %in% y), pre)), simplify = FALSE))

    }, simplify = FALSE)

    ii1 <- mapply(ind1, ii0, FUN = function(ind1, ii0){

      unlist(sapply(ind1, function(y)(setdiff(which(data[, "ID"] %in% y), ii0)), simplify = FALSE))

    }, SIMPLIFY = FALSE)

    # Step 3. Estimating individual effects and ATE for each bootstrap iteration
    ate_boot <- mapply(x = ii0, y = ii1, function(x, y){ate_est(data = data[c(x, y),], int_date = int_date, best = bestt, metric = metric, y.lag = y.lag, ran.err = TRUE)$ate})


  }

  ### Estimating ATE from the matrix of individual effects
  ate_boot <- matrix(ate_boot, nrow = length(unique(data[post, "Time"])))
  rownames(ate_boot) <- unique(data[post, "Time"])
  conf.ate <- apply(ate_boot, 1, quantile, probs = c(alpha/2, 1 - alpha/2))
  var.ate <- apply(ate_boot, 1, var)

  ### Adjusting for bias and/or skewness (if 'type' is "bc classic", "bc block" or "bca")

  if(type %in% c("bc classic", "bc block")){

    # Bias correction
    z0 <- mapply(x = apply(ate_boot, 1, as.list), y = ate, FUN = function(x,y)(qnorm(sum(x < y)/nboot)), SIMPLIFY = TRUE)
    lower <- pnorm(2*z0 + qnorm(alpha/2))
    upper <- pnorm(2*z0 + qnorm(1 - alpha/2))
    conf.ate <- mapply(x = as.list(lower), y = as.list(upper), z = apply(ate_boot, 1, as.list),
                       FUN = function(x,y,z){quantile(unlist(z), probs = c(x,y))}, SIMPLIFY = TRUE)

  }

  if(type == "bca"){

    counts <- t(apply(ii, 1, FUN = function(x)(table(c(x, pre))-1)))
    Blist <- mapply(x = c(1,2,3), y = ate, FUN = function(x,y){
      list(Y = counts, tt = ate_boot[x, ], t0 = y)}, SIMPLIFY = FALSE)
    out2 <- mapply(B = Blist, FUN = bcajack2, MoreArgs = list(alpha = alpha), SIMPLIFY = FALSE)
    conf.ate <- sapply(out2, FUN = function(x)(x$lims[c(1,3), "bca"]))

  }

  # Returning results
  return(list(type = type, conf.ate = conf.ate, var.ate = var.ate, ate.lower = conf.ate[1, ], ate.upper = conf.ate[2, ]))
}

#' Bootstrap inference for CATE
#'
#' Internal function, used within the MLCM routine for the estimation of CATE
#' standard errors and 95% confidence intervals in each terminal node of the tree.
#' It works by resampling the observations at each final node of the tree.
#' Note that the observations are the estimated individual
#' causal effects (computed by comparing the observed data with the ML predictions).
#' Our estimand of interest is the average of the individual effects, so at each bootstrap
#' iteration we average the individual effects, obtaining a bootstrap distribution for the ATE
#' (which is in fact a CATE as we do that in each terminal node, i.e., conditionally on covariates).
#'
#' @param effect Numeric vector of estimated individual causal effects.
#' @param cate Object of class \code{rpart}, the estimated regression-tree-based CATE.
#' @param nboot Number of bootstrap replications.
#' @param alpha Confidence interval level to report for the ATE. Defaulting to 0.05 for a two sided
#'              95\% confidence interval
#'
#' @return A matrix containing the following information: estimated CATE within
#' each node, estimated variance and confidence interval (upper and lower bound)
#' estimated by bootstrap. Each column corresponds to a terminal node of the tree.
#' @noRd

boot_cate <- function(effect, cate, nboot, alpha){

  ### Param checks
  if(!is.numeric(effect)) stop("effect must be a numeric vector")
  if(class(cate) != "rpart") stop ("cate must be an 'rpart' object")
  if(nboot < 1 | all(!class(nboot) %in% c("numeric", "integer")) | nboot%%1 != 0) stop("nboot must be an integer greater than 1")
  if(alpha < 0 | alpha > 1) stop("Invalid confidence interval level, alpha must be positive and less than 1")

  ### Bootstrapping
  terminal.nodes <- cate$where
  x <- unique(terminal.nodes)
  node.inf <- mapply(x, FUN = function(x){y <- effect[which(terminal.nodes == x)];
                                          boot.dist <- matrix(sample(y, size = nboot*length(y), replace = TRUE),
                                                              nrow = nboot, ncol = length(y));
                                          mean.cate <- rowMeans(boot.dist);
                                          var.cate <- var(mean.cate);
                                          conf.cate <- quantile(mean.cate, probs = c(alpha/2, 1 - alpha/2));
                                          c(cate = mean(y), var.cate = var.cate, cate.lower = conf.cate[1], cate.upper = conf.cate[2])},
                     SIMPLIFY = TRUE)
  colnames(node.inf) <- paste0("Node_", x)
  return(node.inf)

}
