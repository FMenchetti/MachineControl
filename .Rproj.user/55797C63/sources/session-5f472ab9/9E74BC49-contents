###############################################################################
##
## MLCM function
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: July 2023
##
###############################################################################

#' Machine Learning Control Method
#'
#' This function is the main workhorse of the package 'MachineControl'. It takes as
#' input a panel dataset, i.e., multiple units observed at several points in time and
#' exposed simultaneously to some policy (indicated by post_period). It then performs
#' the Panel Cross Validation (PCV) by comparing the predictive performance of several Machine
#' Learning (ML) methods in the pre-intervention periods and then outputs the estimated
#' Average Treatment Effect (ATE) of the policy and its confidence interval estimated by
#' bootstrap. Details on the causal assumptions, the estimation process and inference
#' can be found in Cerqua A., Letta M., and Menchetti F. (2023). The method is especially
#' suited when there are no control units available, but if there are control units these
#' can be easily added as control series among the covariates.
#'
#'
#' @param data        A panel dataset in long form, having one column for the time variable, one column for the units'
#'                    unique IDs, one column for the outcome variable and one or more columns for the covariates.
#' @param y           Character, name of the column containing the outcome variable. It can be omitted for \code{PanelMLCM} objects.
#' @param timevar     Character, name of the column containing the time variable. It can be omitted for \code{PanelMLCM} objects.
#' @param id          Character, name of the column containing the ID's. It can be omitted for \code{PanelMLCM} objects.
#' @param post_period The post-intervention period where the causal effect should be computed. It must be contained in 'timevar'.
#' @param inf_type    Character, type of inference to be performed. Possible choices are 'classic', 'block', 'bc classic', 'bc block', 'bca'
#' @param nboot       Number of bootstrap replications, defaults to 1000.
#' @param pcv_block   Number of pre-intervention times to block for panel cross validation. Defaults to 1, see Details.
#' @param metric      Character, the performance metric that should be used to select the optimal model.
#'                    Possible choices are either \code{"RMSE"} (the default) or \code{"Rsquared"}.
#' @param PCV         Optional, best performing ML method as selected from a previous call to \code{PanelCrossValidation}.
#' @param CATE        Whether the function should estimate also CATE (defaults to \code{FALSE}). See Details.
#' @param alpha       Confidence interval level to report for the ATE. Defaulting to 0.05 for a two sided
#'                    95\% confidence interval.
#' @details
#' The panel \code{data} must at least include the response variable, a column of the time variable,
#' and a column with the unit identifiers. Lagged values of the outcome, lagged or contemporaneous
#' exogenous covariates, control series can also be added in the panel dataset. The ML algorithms
#' will automatically treat every column that is left as additional information to improve the
#' prediction of the counterfactual post-intervention outcome.
#'
#' By default, the function internally uses \code{as.PanelMLCM} to organize a messy panel dataset and
#' then compares the performance of two linear models (Partial Least Squares and LASSO)
#' and two non-linear models (Random Forests and Stochastic Gradient Boosting) in one-step ahead
#' predictions by running internally the panel cross validation routine. The comparison is performed based on
#' the \code{metric} provided by the users and lasts until the end of the pre-intervention period.
#' Different ML algorithms can be selected by the user among all the options available in the \code{caret}
#' package. In this case, the users must start from a tidy panel dataset (from a previous call to \code{as.PanelMLCM})
#' and must execute their own panel cross validation. See the examples below.
#'
#' The user can choose among many different bootstrap algorithms. For details see the documentation
#' of the function \code{boot_ate}. For details on the PCV see the documentation of the function
#' \code{PanelCrossValidation}.
#'
#' To speed up the PCV process, the user can 'block' some pre-intervention periods by increasing the
#' \code{pcv_block} parameter, e.g., \code{pcv_block = 1} (the default) indicates to use the observations
#' in the first time period as the first training sample and to test on the next period. Then, the second
#' training sample will be formed by the observations on the first two time periods. Validation will be
#' performed on the third period and so on. For longer series, specifying \code{pcv_block > 1} reduces computational time.
#' For example, by setting \code{pcv_block = 4} when the length of the pre-intervention time series is 7 reduces the number
#' of validation sets to 3 instead of 6.
#'
#' By default, the function estimates an Average Treatment Effect (ATE), but when \code{CATE = TRUE} it also estimates
#' a Conditional Average Treatment Effect (CATE) as described in Cerqua A., Letta M. & Menchetti F.
#' <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4315389>. In short, CATE is estimated with a regression tree
#' and it is of interest when researchers suspect that the policy ('treatment' or 'intervention') has produced
#' heterogeneous effects on the units in the panel. For additional details on the causal estimands, estimation process and
#' underlying assumptions see the paper cited above.
#'
#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{best_method}: the best-performing ML algorithm as selected by the PCV routine
#'   \item \code{fit}: the final result of the training step of the best-performing ML algorithm
#'   on all pre-intervention data
#'   \item \code{ate}: the estimated ATE
#'   \item \code{var.ate}: the variance of ATE as estimated by the bootstrap algorithm selected
#'   by the user in 'inf_type'
#'   \item \code{ate.lower}: lower bound of a 95% bootstrap confidence interval
#'   \item \code{ate.upper}: upper bound of a 95% bootstrap confidence interval
#'   \item \code{cate}: a list of class \code{rpart} with the estimated regression-tree-based CATE
#'   \item \code{cate.inf}: a matrix containing the variance of CATE and the confidence interval
#'   estimated by bootstrap
#' }
#' @export
#'
#' @examples
#'
#' ### Example 1. Estimating ATE and CATE (with default ML methods)
#'
#' # Estimation
#' fit <- MLCM(data = data, y = "Y", timevar = "year", id = "ID", post_period = 2020, inf_type = "classic", nboot = 10, CATE = TRUE)
#'
#' # ATE & CATE
#' fit$ate
#' plot(fit$cate) ; text(fit$cate)
#'
#' # Bootstrap confidence interval
#' c(fit$ate.lower, fit$ate.upper)
#'

MLCM <- function(data, y = NULL, timevar = NULL, id = NULL, post_period, inf_type, nboot = 1000, pcv_block = 1, metric = "RMSE", PCV = NULL, CATE = FALSE, alpha = 0.05){

  ### Parameter checks
  if(!any(class(data) %in% c("matrix", "data.frame", "PanelMLCM"))) stop("data must be a matrix, a data.frame or a PanelMLCM object")
  if(!"PanelMLCM" %in% class(data) & any(c(is.null(y), is.null(timevar), is.null(id)))) stop("Unspecified columns in 'matrix' or 'data.frame' data")
  if(!(is.null(y) | class(y) == "character")) stop("y must be a character")
  if(!(is.null(timevar) | class(timevar) == "character")) stop("timevar must be a character")
  if(!(is.null(id) | class(id) == "character")) stop("id must be a character")
  if(!is.null(y)){if(!y %in% colnames(data)) stop (paste("there is no column called", y, "in 'data'"))}
  if(!is.null(timevar)){if(!timevar %in% colnames(data)) stop (paste("there is no column called", timevar, "in 'data'"))}
  if(!is.null(id)){if(!id %in% colnames(data)) stop (paste("there is no column called", id, "in 'data'"))}
  if(!any(class(post_period) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt", "numeric", "integer"))) stop("post_period must be integer, numeric or Date")
  if(is.null(timevar)){if(!post_period %in% data[, "Time"]) stop ("post_period must be contained in the 'Time' column")}
  if(!is.null(timevar)){if(!post_period %in% data[, timevar]) stop ("post_period must be contained in timevar")}
  if(!any(inf_type %in% c("classic", "block", "bc classic", "bc block", "bca"))) stop("Inference type not allowed, check the documentation")
  if(nboot < 1 | all(!class(nboot) %in% c("numeric", "integer")) | nboot%%1 != 0) stop("nboot must be an integer greater than 1")
  if(!metric %in% c("RMSE", "Rsquared")) stop("Metric not allowed, check documentation")
  if(!is.null(PCV)){if(!"train" %in% class(PCV)) stop ("Invalid PCV method, it should be an object of class 'train'")}
  if(alpha < 0 | alpha > 1) stop("Invalid confidence interval level, alpha must be positive and less than 1")

  ### Structuring the panel dataset in the required format
  if("PanelMLCM" %in% class(data)){

    data_panel <- data

  } else {

    data_panel <- as.PanelMLCM(y = data[, y], timevar = data[, timevar], id = data[, id],
                               x = data[, !(names(data) %in% c(y, id, timevar))])

  }


  ### Panel cross-validation
  if(is.null(PCV)){

    best <- PanelCrossValidation(data = data_panel, post_period = post_period, pcv_block = pcv_block, metric = metric)

  } else {

    best <- PCV

  }


  ### Fit the best (optimized) ML algorithm on all pre-intervention data and make predictions in the post-intervention period
  ind <- which(data_panel[, "Time"] < post_period)

  set.seed(1)

  invisible(capture.output(
    fit <- train(Y ~ .,
                 data = data_panel[ind, !(names(data_panel) %in% c("ID", "Time"))],
                 method = best$method,
                 metric = metric,
                 trControl = trainControl(method="none"),
                 tuneGrid = best$bestTune)
  ))

  obs <- data_panel[-ind, "Y"]
  pred <- predict(fit, newdata = data_panel[-ind, ])

  ### ATE
  ate <- mean(obs - pred)

  ### Inference
  invisible(capture.output(

    boot_inf <- boot_ate(data = data_panel, ind = ind, bestt = fit, type = inf_type, nboot = nboot, ate = ate, alpha = alpha)

  ))

  ### CATE
  if(CATE){

    data_cate <- data.frame(effect = obs - pred, data_panel[-ind, !names(data_panel) %in% c("Y","Time","ID")])
    cate <- rpart(effect ~ ., method="anova", data = data_cate, cp = 0, minbucket = 0.05*length(obs))
    cate.inf <- boot_cate(effect = obs - pred, cate = cate, nboot = nboot, alpha = alpha)

  } else {

    cate <- NULL
    cate.inf <- NULL
  }


  #options(scipen=999)
  #prp(model.rpart, fallen.leaves = FALSE, box.col="lightgray", type = 3, branch=1, branch.lty= 1,   main="Data-driven CATEs")

  ### Saving results
  return(list(best_method = best, fit = fit, ate = ate, var.ate = boot_inf$var.ate, ate.lower = boot_inf$ate.lower, ate.upper = boot_inf$ate.upper, cate = cate, cate.inf = cate.inf))
  # return(list(best_method = best, fit = fit, ate = ate, ate.lower = conf.ate[1], ate.upper = conf.ate[2]), conf.individual = conf.individual)
}

#' Structuring the panel dataset
#'
#' This function takes as input the panel dataset given by the user and changes the
#' ordering and the names of the columns to obtain an object of class 'PanelMLCM'
#' to be used by the function 'MLCM'.
#'
#' @param y Numeric, the outcome variable.
#' @param x Matrix or data.frame of covariates to include in the model.
#' @param timevar  The column containing the time variable. It can be numeric, integer or
#'                 a date object
#' @param id Numeric, the column containing the ID's.
#'
#' @details This function is mainly for internal use of \code{MLCM}. It is exported to give
#' users full flexibility in the choice of the ML algorithms and during the panel cross validation.
#' See the documentation of the function \code{PanelCrossValidation}.
#'
#' @return An object of class 'data.frame' and 'PanelMLCM'.
#' @export
#' @examples
#'
#' # Start from a disorganized panel dataset
#' head(data)
#'
#' # Run as.PanelMLCM
#' newdata <- as.PanelMLCM(y = data[, "Y"], timevar = data[, "year"], id = data[, "ID"],
#'                         x = data[, !(names(data) %in% c("Y", "ID", "year"))])
#'
#' # Results
#' head(newdata)

as.PanelMLCM <- function(y, x, timevar, id){

  # Parameter checks
  if(!(is.numeric(y) & length(y)>1)) stop("y must be a numeric vector of length greater than 1")
  if(!any(class(x) %in% c("numeric", "matrix", "data.frame"))) stop("x must be a vector, matrix or data.frame")
  if(NROW(x) != length(y)) stop("NROW(x) != length(y)")
  if(!any(class(timevar) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt", "integer", "numeric")) | length(timevar) != length(y)) stop("timevar must be a numeric vector or a 'Date' object of the same length as y")
  if(!(is.numeric(id) & length(id) == length(y))) stop("id must be a numeric vector of the same length as y")
  if(length(unique(id))*length(unique(timevar)) != length(y)) warning("The panel is unbalanced")

  # Structuring the panel dataset
  panel <- data.frame(Time = timevar, ID = id, Y = y, x)

  # Returning results
  class(panel) <- c("data.frame", "PanelMLCM")
  return(panel)
}
