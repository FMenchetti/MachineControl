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
#' @param y           Character, name of the column containing the outcome variable.
#' @param timevar     Character, name of the column containing the time variable.
#' @param id          Character, name of the column containing the ID's.
#' @param post_period The post-intervention period where the causal effect should be computed. It must be contained in 'timevar'.
#' @param inf_type    Character, type of inference to be performed. Possible choices are 'classic', 'block', 'bc classic', 'bc block', 'bca'
#' @param nboot       Number of bootstrap replications, defaults to 1000.
#' @param pcv_block   Number of pre-intervention times to block for panel cross validation. Defaults to 1, see Details.
#' @param metric      Character, the performance metric that should be used to select the optimal model.
#'                    Possible choices are either \code{"RMSE"} (the default) or \code{"Rsquared"}.
#'
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
#' of the function \code{boot_fun}. For details on the PCV see the documentation of the function
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
#' }
#' @export
#'
#' @examples
#'
#' ### Example 1. Running MLCM with the default options
#'
#' # Causal effect estimation of a policy occurred in 2020
#' fit <- MLCM(data = data, y = "Y", timevar = "year", id = "ID", post_period = 2020, inf_type = "classic", nboot = 10)
#' fit$ate
#'
#' # Bootstrap confidence interval
#' c(fit$ate.lower, fit$ate.upper)
#'
MLCM <- function(data, y, timevar, id, post_period, inf_type, nboot = 1000, pcv_block = 1, metric = "RMSE"){

  ### Parameter checks
  if(!any(class(data) %in% c("matrix", "data.frame"))) stop("data must be a matrix or a data.frame")
  if(class(y) != "character") stop("y must be a character")
  if(!(y %in% colnames(data))) stop (paste("there is no column called", y, "in 'data'"))
  if(class(timevar) != "character") stop("timevar must be a character")
  if(!(timevar %in% colnames(data))) stop (paste("there is no column called", timevar, "in 'data'"))
  if(class(id) != "character") stop("id must be a character")
  if(!(id %in% colnames(data))) stop (paste("there is no column called", id, "in 'data'"))
  if(!any(class(post_period) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt", "numeric", "integer"))) stop("post_period must be integer, numeric or Date")
  if(!any(inf_type %in% c("classic", "block", "bc classic", "bc block", "bca"))) stop("Inference type not allowed, check the documentation")

  ### Structuring the panel dataset in the required format
  data_panel <- as.PanelMLCM(y = data[, y], timevar = data[, timevar], id = data[, id],
                             x = data[, !(names(data) %in% c(y, id, timevar))])

  ### Panel cross-validation
  best <- PanelCrossValidation(data = data_panel, post_period = post_period, pcv_block = pcv_block, metric = metric)

  ### Fit the best (optimized) ML algorithm on all pre-intervention data and make predictions in the post-intervention period
  ind <- which(data_panel[, "Time"] < post_period)

  set.seed(1)
  fit <- train(Y ~ .,
               data = data_panel[ind, !(names(data_panel) %in% c("ID", "Time"))],
               method = best$method,
               metric = metric,
               trControl = trainControl(method="none"),
               tuneGrid = best$bestTune)
  obs <- data_panel[-ind, "Y"]
  pred <- predict(fit, newdata = data_panel[-ind, ])

  ### ATE
  ate <- mean(obs - pred)

  ### Inference
  boot_inf <- boot_fun(data = data_panel, ind = ind, bestt = fit, type = inf_type, nboot = nboot, ate = ate)

  ### Saving results
  return(list(best_method = best, fit = fit, ate = ate, var.ate = boot_inf$var.ate, ate.lower = boot_inf$ate.lower, ate.upper = boot_inf$ate.upper))
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
