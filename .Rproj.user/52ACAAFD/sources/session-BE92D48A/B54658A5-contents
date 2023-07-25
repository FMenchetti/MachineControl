###############################################################################
##
## Bootstrap inference
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: July 2023
##
###############################################################################

# Bootstrap inference for ATE
#
# Internal function, used within the MLCM routine for the estimation of ATE
# standard error and 95% confidence interval.
#
# @param data A 'PanelMLCM' object from a previous call to \code{as.PanelMLCM}.
# @param ind  Vector indicating pre-intervention time points.
# @param bestt Object of class \code{train}, the best-performing ML method as selected
#              by panel cross validation.
# @param type Character, type of inference to be performed. Possible choices are 'classic', 'block', 'bc classic', 'bc block', 'bca'.
# @param nboot Number of bootstrap replications.
# @param ate Numeric, the estimated ATE in the sample.
#
# @return A list with the following components:
# \itemize{
#   \item \code{type}: the inference type that has been performed
#   \item \code{ate.boot}: the bootstrap distribution for ATE
#   \item \code{conf.ate}: bootstrap confidence interval at the 95% level
#   \item \code{var.ate}: estimated variance for ATE
#   \item \code{ate.lower}: lower confidence interval bound
#   \item \code{ate.upper}: upper confidence interval bound
# }


boot_ate <- function(data, ind, bestt, type, nboot, alpha, ate = NULL){

  ### Param checks
  if(!any(class(data) %in% "PanelMLCM")) stop("Invalid class in the PanelCrossValidation function, something is wrong with as.PanelMLCM")
  if(length(ind)<1) stop("A zero-length pre-intervention period was selected, please check your data and your definition of 'post_period'")
  if(!any(type %in% c("classic", "block", "bc classic", "bc block", "bca"))) stop("Inference type not allowed, check the documentation")
  if(nboot < 1 | all(!class(nboot) %in% c("numeric", "integer")) | nboot%%1 != 0) stop("nboot must be an integer greater than 1")
  if(alpha < 0 | alpha > 1) stop("Invalid confidence interval level, alpha must be positive and less than 1")
  if(!is.numeric(ate)) stop ("ATE must be numeric")

  ### Step 1. Sampling indices
  if(type %in% c("classic", "bc classic", "bca")){

    ii<- matrix(sample(ind, size = nboot*length(ind), replace = T), nrow = nboot, ncol = length(ind))

  }

  if(type %in% c("block", "bc block")){

    ids <- unique(data$ID)
    ii <- t(sapply(1:nboot, function(boot){set.seed(boot); ind1 <- sample(ids, size = length(ids), replace = T);
                                           unlist(sapply(ind1, function(x)(intersect(which(data[, "ID"] %in% x), ind)), simplify = F))}))
  }

  ### Step 2. Bootstrapping (estimating the causal effect on the units resampled in each bootstrap iteration)
  ate_boot <- apply(ii, 1, function(i){
    fitb <- train(Y ~ .,
                  data = data[i, !(names(data) %in% c("ID", "Time"))],
                  method = bestt$method,
                  metric = "RMSE",
                  trControl = trainControl(method="none"),
                  tuneGrid = bestt$bestTune) ;
    obs <- data[-c(ind, i), "Y"] ;
    eps <- data[i, "Y"] - predict(fitb)
    # error <- sample(eps, size = nrow(data[-ind,]), replace = T) # verifica che sia così anche per block boot
    error <- rnorm(n = nrow(data[-ind,]), mean = mean(eps), sd = sd(eps))
    pred <- predict(fitb, newdata = data[-c(ind, i), ]) + error ;
    #pred <- predict(fitb, newdata = data[-c(ind, i), ]) ;
    obs - pred
  })

  mean_ate_boot <- colMeans(ate_boot)
  conf.ate <- quantile(mean_ate_boot, probs = c(alpha/2, 1 - alpha/2))

  ### Step 3. Adjusting for bias and/or skewness (if 'type' is "bc classic", "bc block" or "bca")

  if(type %in% c("bc classic", "bc block")){

    # Bias correction
    z0 <- qnorm(sum(mean_ate_boot < ate)/nboot)
    lower <- pnorm(2*z0 + qnorm(alpha/2))
    upper <- pnorm(2*z0 + qnorm(1 - alpha/2))
    conf.ate <- quantile(mean_ate_boot, probs = c(lower, upper))

  }

  if(type == "bca"){

    counts <- t(apply(ii, 1, FUN = function(x)(table(c(x, ind))-1)))
    Blist <- list(Y = counts, tt = colMeans(ate_boot), t0 = ate)
    out2 <- bcajack2(B = Blist)
    conf.ate <- out2$lims[c(1,9),"bca"]

  }


  # Returning results
  return(list(type = type, ate.boot = ate_boot, conf.ate = conf.ate, var.ate = var(mean_ate_boot), ate.lower = conf.ate[1], ate.upper = conf.ate[2]))
}

# Bootstrap inference for CATE
#
# Internal function, used within the MLCM routine for the estimation of CATE
# standard errors and 95% confidence intervals in each terminal node of the tree.
# It works by resampling the observations at each final node of the tree.
# Note that the observations are the estimated individual
# causal effects (computed by comparing the observed data with the ML predictions).
# Our estimand of interest is the average of the individual effects, so at each bootstrap
# iteration we average the individual effects, obtaining a bootstrap distribution for the ATE
# (which is in fact a CATE as we do that in each terminal node, i.e., conditionally on covariates).
#
# @param effect Numeric vector of estimated individual causal effects.
# @param cate Object of class \code{rpart}, the estimated regression-tree-based CATE.
# @param nboot Number of bootstrap replications.
#
# @return A matrix containing the following information: estimated CATE within
# each node, estimated variance and confidence interval (upper and lower bound)
# estimated by bootstrap. Each column corresponds to a terminal node of the tree.

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
