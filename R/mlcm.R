###############################################################################
##
## MLCM function
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: July 2023
##
###############################################################################

#' MLCM
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
#'
#' pcv_block = 1 (the default) indicates to use the observations in the first time period as
#' the first training sample and to test on the next period. Then, the second training sample
#' will be formed by the observations on the first two time periods. Validation will be
#' performed on the the third period and so on. For longer series, specifying pcv_block > 1 reduces computational time.
#' For example, by setting pcv_block = 4 when the length of the pre-intervention time series is 7 reduces the number
#' of vlidation sets to 3 instead of 6.
#'
#' @return
#' @export
#'
#' @examples
#'
MLCM <- function(data, y, timevar, id, post_period, inf_type, nboot = 1000, pcv_block = 1){

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
  best <- PanelCrossValidation(data = data_panel, post_period = post_period, blocked = pcv_block)

  ### Fit the best (optimized) ML algorithm on all pre-intervention data and make predictions in the post-intervention period
  ind <- which(data_panel[, "Time"] < post_period)

  set.seed(1)
  fit <- train(Y ~ .,
               data = data_panel[ind, !(names(data_panel) %in% c("ID", "Time"))],
               method = best$method,
               metric = "RMSE",
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
