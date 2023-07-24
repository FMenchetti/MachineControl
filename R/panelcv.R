###############################################################################
##
## Panel Cross Validation
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: July 2023
##
###############################################################################

#' Panel Cross Validation
#'
#' @param data A 'PanelMLCM' object from a previous call to \code{as.PanelMLCM}.
#' @param post_period The post-intervention period where the causal effect should be computed.
#' @param blocked Number of pre-intervention times to block for panel cross validation. Defaults to 1, see Details.
#' @param metric Character, the performance metric that should be used to select the optimal model.
#'               Possible choices are either \code{"RMSE"} (the default) or \code{"Rsquared"}.
#' @param trControl Optional, used to customize the training step. It must be the output from a call to \code{trainControl} from the \code{caret} package.
#' @param ML_methods Optional list of ML methods to be used as alternatives to the default methods. Each method must be supplied
#'                   as a named list of two elements: a character defining the method name from all the ones available in \code{caret}
#'                   and the grid of parameter values to tune via the panel cross validation. See Details and the examples for additional explanations.#'
#' @return A list of class \code{train} with the best-performing ML method.
#' @export
#'
#' @examples
#'
#' ### Example 1. Changing the default training method
#'
#' # Organizing the dataset
#' newdata <- as.PanelMLCM(y = data[, "Y"], timevar = data[, "year"], id = data[, "ID"],
#'                         x = data[, !(names(data) %in% c("Y", "ID", "year"))])
#'
#' # Using the first two years for training and the last two years for testing
#' indices <- CreateSpacetimeFolds(newdata, timevar = "Time", k = 5)
#' trainx <- indices$indexOut[1:2]
#' testx <- indices$indexOut[3:4]
#' ctrl <- trainControl(index = trainx, indexOut = testx)
#'
#' # Customized panel cross validation
#' pcv <- PanelCrossValidation(data = newdata, post_period = 2020, trControl = ctrl)
#'
#' ### Example 2. Changing ML methods
#' enet <- list(method = "enet",
#'             tuneGrid = expand.grid(
#'               fraction = seq(0.1, 0.9, by = 0.1),
#'               lambda = seq(0.1, 0.9, by = 0.1)))
#'
#' linreg <- list(method = "lm",
#'               tuneGrid = expand.grid(
#'                 intercept = seq(0, 10, by = 0.5)))
#'
#' pcv <- PanelCrossValidation(data = newdata, post_period = 2020, ML_methods = list(enet, linreg))
#'
#' ### Example 3. Changing ML methods and trControl
#' pcv <- PanelCrossValidation(data = newdata, post_period = 2020, trControl = ctrl, ML_methods = list(enet, linreg))

PanelCrossValidation <- function(data, post_period, blocked = 1, metric = "RMSE", trControl = NULL, ML_methods = NULL){

  ### Parameter checks
  if(!any(class(data) %in% "PanelMLCM")) stop("Invalid class in the PanelCrossValidation function, something is wrong with as.PanelMLCM")
  if(!(post_period %in% data[, "Time"])) stop("post_period must be contained in timevar")
  if(length(unique(data[, "Time"])) - 2 - blocked < 1) stop("Panel cross validation must be performed in at least one time period")
  if(blocked <= 0) stop("The number of 'blocked' time periods for panel cross validation must be at least 1")
  if(!metric %in% c("RMSE", "Rsquared")) stop("Metric not allowed, check documentation")
  if(is.list(ML_methods)){if(any(sapply(ML_methods, FUN = length) != 2)) stop("'ML_methods' must be a list of methods, each of length 2")}
  if(is.list(ML_methods)){
    if(any(sapply(ML_methods, FUN = function(x)(any(!names(x) %in% c("method", "tuneGrid")))))) stop("Each method in 'ML_methods' must be a named list, check documentation")}

  ### STEP 1. The CAST package is used to generate separate testing sets for each year
  Tt <- length(unique(data[, "Time"]))
  post_period <- post_period
  indices <- CreateSpacetimeFolds(data, timevar = "Time", k = Tt)
  trainx <- lapply(blocked:(Tt-2), FUN = function(x) unlist(indices$indexOut[1:x]))
  testx <- lapply(blocked:(Tt-2), FUN = function(x) unlist(indices$indexOut[[x+1]]))

  ### STEP 2. Set control function by specifying the training and testing folds that caret will use
  ###         for cross-validation and tuning of the hyperparameters (i.e., the combination of folds defined above)
  if(is.null(trControl)){

    ctrl <- trainControl(index = trainx, indexOut = testx)

  } else {

    ctrl <- trControl # da capire come fare i param checks

  }


  ### STEP 3.  Tune the hyperparameters of each of the ML algorithms via temporal cross-validation

  if(is.null(ML_methods)){

    # STOCHASTIC GRADIENT BOOSTING
    gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3),
                            n.trees=c(500, 1000, 1500, 2000),
                            shrinkage = seq(0.01, 0.1, by = 0.01),
                            n.minobsinnode = c(10,20))

    set.seed(1)
    bo <- train(Y ~ .,
                data = data[, !(names(data) %in% c("ID", "Time"))],
                method = "gbm",
                metric = metric,
                trControl = ctrl,
                tuneGrid = gbmGrid,
                verbose = FALSE)

    # RANDOM FOREST
    set.seed(1)
    rf <- train(Y ~ .,
                data = data[, !(names(data) %in% c("ID", "Time"))],
                method = "rf",
                metric = metric,
                search = "grid",
                trControl = ctrl,
                tuneGrid = expand.grid(mtry = (2:(ncol(data)-3))),
                ntree=500)
    # LASSO
    lasso <- train(Y ~ .,
                   data = data[, !(names(data) %in% c("ID", "Time"))],
                   method = "lasso",
                   metric = metric,
                   trControl = ctrl,
                   tuneGrid = expand.grid(fraction = seq(0.1, 0.9, by = 0.1)),
                   preProc=c("center", "scale"))

    # PLS
    pls <- train(Y ~ .,
                 data = data[, !(names(data) %in% c("ID", "Time"))],
                 method = "pls",
                 metric = metric,
                 trControl = ctrl,
                 tuneGrid = expand.grid(ncomp = c(1:10)),
                 preProc=c("center", "scale"))

    # Storing results in a list
    m_list <- list(bo = bo, rf = rf, lasso = lasso, pls = pls)

  } else {

   m_list <- lapply(ML_methods, FUN = function(x){train(Y ~.,
                                                  data = data[, !(names(data) %in% c("ID", "Time"))],
                                                  method = x$method,
                                                  metric = metric,
                                                  trControl = ctrl,
                                                  tuneGrid = x$tuneGrid)})

  }


  ### STEP 4. Selecting the "best" ML algorithm based on the provided performance metric
  rmse_min <- sapply(m_list, FUN = function(x) min(x$results[, metric]), simplify = T)
  ind <- which(rmse_min == min(rmse_min))

  ### Returning result
  return(best = m_list[[ind]])

}
