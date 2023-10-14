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
#' This function implements the panel cross validation technique as described in
#' Cerqua A., Letta M. & Menchetti F., (2023) <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4315389>.
#' It works as a rolling window training and testing procedure where the ML methods are trained recursively on the
#' the observations in the first $t$ time periods and tested in the prediction of the next periods.
#'
#' @param data A 'PanelMLCM' object from a previous call to \code{as.PanelMLCM}.
#' @param int_date The date of the intervention, treatment, policy introduction or shock.
#' @param pcv_block Number of pre-intervention times to block for panel cross validation. Defaults to 1, see Details.
#' @param metric Character, the performance metric that should be used to select the optimal model.
#'               Possible choices are either \code{"RMSE"} (the default) or \code{"Rsquared"}.
#' @param trControl Optional, used to customize the training step. It must be the output from a call to \code{trainControl} from the \code{caret} package.
#' @param ML_methods Optional list of ML methods to be used as alternatives to the default methods. Each method must be supplied
#'                   as a named list of two elements: a character defining the method name from all the ones available in \code{caret}
#'                   and the grid of parameter values to tune via the panel cross validation. See Details and the examples for additional explanations.#'
#' @return A list of class \code{train} with the best-performing ML method.
#'
#' @details
#' To speed up computational time in case of longer time series, users can increase the \code{pcv_block} parameter: \code{pcv_block = 1}
#' (the default) indicates to use the observations in the first time period as the first training sample and to test on the next period.
#' Then, the second training sample will be formed by the observations on the first two time periods and validation will be performed
#' on the third period and so on. For longer series, specifying \code{pcv_block > 1} reduces computational time. For example,
#' by setting \code{pcv_block = 4} when the length of the pre-intervention time series is 7 reduces the number of validation sets to 3 instead of 6.
#'
#' By default, the panel cross validation routine compares two linear models (Partial Least Squares and LASSO) and two non-linear models
#' (Random Forests and Stochastic Gradient Boosting) in the prediction of the outcome in the testing sets. The default specification
#' is used internally by the \code{MLCM} function that records the best-performing ML algorithm and use it to estimate the ATE of the treatment.
#' Users can also change the default ML methods by the argument \code{ML_method}, a list of as many components as the ML methods that the users want
#' to try. Each element of \code{ML_method} should be a list of two elements: \code{method} must be a character with the method name from all the
#' ones available in \code{caret}; \code{tuneGrid} is a grid of parameter values to tune via the panel cross validation. See Example 2 for an
#' illustration.
#'
#' @export
#' @import caret
#' @importFrom rpart rpart
#' @importFrom bcaboot bcajack2
#' @importFrom CAST CreateSpacetimeFolds
#' @importFrom gbm gbm
#' @importFrom elasticnet enet
#' @importFrom pls plsr
#' @importFrom randomForest randomForest
#' @importFrom utils capture.output
#' @importFrom stats quantile
#' @importFrom stats qnorm
#' @importFrom stats var
#' @importFrom stats pnorm
#' @importFrom stats rnorm
#' @importFrom stats sd
#'
#' @examples
#'
#' ### Example 1. Changing the default training method
#'
#' # Organizing the dataset
#' newdata <- as.PanelMLCM(y = data[, "Y"], timevar = data[, "year"], id = data[, "ID"],
#'                         x = data[, !(names(data) %in% c("Y", "ID", "year"))], y.lag = 1)
#'
#' # Using the first two years for training and the last two years for testing
#' indices <- CAST::CreateSpacetimeFolds(newdata, timevar = "Time", k = length(unique(newdata$Time)))
#' trainx <- indices$indexOut[1:2]
#' testx <- indices$indexOut[3:4]
#' ctrl <- trainControl(index = trainx, indexOut = testx)
#'
#' # Customized panel cross validation
#' pcv <- PanelCrossValidation(data = newdata, int_date = 2019, trControl = ctrl)
#'
#' ### Example 2. Changing ML methods and estimating ATE
#'
#' enet <- list(method = "enet",
#'             tuneGrid = expand.grid(
#'               fraction = seq(0.1, 0.9, by = 0.1),
#'               lambda = seq(0.1, 0.9, by = 0.1)))
#'
#' linreg <- list(method = "lm",
#'               tuneGrid = expand.grid(
#'                 intercept = seq(0, 10, by = 0.5)))
#'
#' pcv <- PanelCrossValidation(data = newdata, int_date = 2019, trControl = ctrl,
#'                             ML_methods = list(enet, linreg))
#'
#' causal <- MLCM(data = newdata, int_date = 2019, inf_type = "classic", PCV = pcv,
#'                nboot = 10, CATE = FALSE)
#'
#' causal$ate
#' causal$conf.ate

PanelCrossValidation <- function(data, int_date, pcv_block = 1, metric = "RMSE", trControl = NULL, ML_methods = NULL){

  ### Parameter checks
  if(!any(class(data) %in% "PanelMLCM")) stop("Invalid class in the PanelCrossValidation function, something is wrong with as.PanelMLCM")
  if(!(int_date %in% data[, "Time"])) stop("int_date must be contained in timevar")
  if(sum(unique(data[, "Time"]) < int_date) - pcv_block < 1) stop("Panel cross validation must be performed in at least one time period")
  #if(length(unique(data[, "Time"])) - 2 - pcv_block < 1) stop("Panel cross validation must be performed in at least one time period")
  if(pcv_block <= 0) stop("The number of 'pcv_block' time periods for panel cross validation must be at least 1")
  if(!metric %in% c("RMSE", "Rsquared")) stop("Metric not allowed, check documentation")
  if(is.list(ML_methods)){if(any(sapply(ML_methods, FUN = length) != 2)) stop("'ML_methods' must be a list of methods, each of length 2")}
  if(is.list(ML_methods)){
    if(any(sapply(ML_methods, FUN = function(x)(any(!names(x) %in% c("method", "tuneGrid")))))) stop("Each method in 'ML_methods' must be a named list, check documentation")}

  ### STEP 1. The CAST package is used to generate separate testing sets for each year
  Tt <- length(unique(data[, "Time"]))
  indices <- CreateSpacetimeFolds(data, timevar = "Time", k = Tt)
  end <- sum(unique(data[, "Time"]) < (int_date - 1) )
  trainx <- lapply(pcv_block:end, FUN = function(x) unlist(indices$indexOut[1:x]))
  testx <- lapply(pcv_block:end, FUN = function(x) unlist(indices$indexOut[[x+1]]))

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

