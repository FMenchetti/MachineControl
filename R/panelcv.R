###############################################################################
##
## Panel Cross Validation
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: August 2024
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
#' @param default_par List of parameters for the default ML algorithms. See Details.
#' @param trControl Optional, used to customize the training step. It must be the output from a call to \code{trainControl} from the \code{caret} package.
#' @param ML_methods Optional list of ML methods to be used as alternatives to the default methods. Each method must be supplied
#'                   as a named list of two elements: a character defining the method name from all the ones available in \code{caret}
#'                   and the grid of parameter values to tune via the panel cross validation. See Details and the examples for additional explanations.
#' @return A list containing the following objects:
#' \itemize{
#'   \item \code{best}: the best-performing ML algorithm
#'   \item \code{best.metric}: the value of \code{metric} for the best-performing ML algorithm
#'   \item \code{all_methods}: detailed results for all the methods tried during the PCV routine
#'  }
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
#' set.seed(1)
#' enet <- list(method = "enet",
#'             tuneGrid = expand.grid(
#'               fraction = seq(0.1, 0.9, by = 0.1),
#'               lambda = seq(0.1, 0.9, by = 0.1)))
#'
#' linreg <- list(method = "lm",
#'               tuneGrid = expand.grid(
#'                 intercept = seq(0, 10, by = 0.5)))
#'
#' pls <- list(method = "pls", tuneGrid = expand.grid(ncomp = c(1:5)))
#'
#' pcv <- PanelCrossValidation(data = newdata, int_date = 2019, trControl = ctrl,
#'                             ML_methods = list(enet, linreg, pls))
#'
#' causal <- MLCM(data = newdata, int_date = 2019, inf_type = "classic", PCV = pcv,
#'                nboot = 10, CATE = FALSE)
#'
#' causal$estimate$ate
#' causal$estimate$conf.ate

PanelCrossValidation <- function(data, int_date, pcv_block = 1, metric = "RMSE",
                                 default_par = list(gbm = list(depth = c(1,2,3),
                                                               n.trees = c(500, 1000),
                                                               shrinkage = seq(0.01, 0.1, by = 0.02),
                                                               n.minobsinnode = c(10,20)),
                                                    rf = list(ntree = 500),
                                                    pls = list(ncomp = c(1:5)),
                                                    lasso = list(fraction = seq(0.1, 0.9, by = 0.1))),
                                 trControl = NULL, ML_methods = NULL){

  ### Parameter checks
  if(!any(class(data) %in% "PanelMLCM")) stop("Invalid class in the PanelCrossValidation function, something is wrong with as.PanelMLCM")
  if(!(int_date %in% data[, "Time"])) stop("int_date must be contained in timevar")
  if(sum(unique(data[, "Time"]) < int_date) - pcv_block < 1) stop("Panel cross validation must be performed in at least one time period")
  if(pcv_block <= 0) stop("The number of 'pcv_block' time periods for panel cross validation must be at least 1")
  if(!metric %in% c("RMSE", "Rsquared")) stop("Metric not allowed, check documentation")
  if(!is.null(ML_methods)){if(!is.list(ML_methods)) stop("ML_methods must be a list")}
  if(is.list(ML_methods)){
    if(any(sapply(ML_methods, FUN = function(x)(!all(c("method", "tuneGrid") %in% names(x)))))) stop("Each method in 'ML_methods' must be a named list, check documentation")}
  if(is.null(ML_methods) & !is.list(default_par)) stop(" 'default_par' must be a list, check the documentation")
  if(is.null(ML_methods) & any(!names(default_par) %in% c("gbm", "rf", "lasso", "pls"))) stop("'default_par' must be a list of parameters for the default gbm, rf, lasso and pls algorithm, check the documentation")

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

    ctrl <- trControl

  }


  ### STEP 3.  Tune the hyperparameters of each of the ML algorithms via temporal cross-validation

  if(is.null(ML_methods)){

    # STOCHASTIC GRADIENT BOOSTING
    gbmGrid <-  expand.grid(interaction.depth = default_par$gbm$depth,
                            n.trees = default_par$gbm$n.trees,
                            shrinkage = default_par$gbm$shrinkage,
                            n.minobsinnode = default_par$gbm$n.minobsinnode)

    #set.seed(1)
    bo <- train(Y ~ .,
                data = data[, !(names(data) %in% c("ID", "Time"))],
                method = "gbm",
                metric = metric,
                trControl = ctrl,
                tuneGrid = gbmGrid,
                verbose = FALSE)

    # RANDOM FOREST
    #set.seed(1)
    rf <- train(Y ~ .,
                data = data[, !(names(data) %in% c("ID", "Time"))],
                method = "rf",
                metric = metric,
                search = "grid",
                trControl = ctrl,
                tuneGrid = expand.grid(mtry = (2:(ncol(data)-3))),
                ntree = default_par$rf$ntree )
    # LASSO
    lasso <- train(Y ~ .,
                   data = data[, !(names(data) %in% c("ID", "Time"))],
                   method = "lasso",
                   metric = metric,
                   trControl = ctrl,
                   tuneGrid = expand.grid(fraction = default_par$lasso$fraction),
                   preProc=c("center", "scale"))

    # PLS
    pls <- train(Y ~ .,
                 data = data[, !(names(data) %in% c("ID", "Time"))],
                 method = "pls",
                 metric = metric,
                 trControl = ctrl,
                 tuneGrid = expand.grid(ncomp = default_par$pls$ncomp),
                 preProc=c("center", "scale"))

    # Storing results in a list
    m_list <- list(bo = bo, rf = rf, lasso = lasso, pls = pls)

  } else {

    invisible(capture.output(

      m_list <- lapply(ML_methods, FUN = function(x)(do.call(train, c(list(Y ~ .,
                                                                           data = data[, !(names(data) %in% c("ID", "Time"))],
                                                                           metric = metric,
                                                                           trControl = ctrl), x))))
    ))

  }


  ### STEP 4. Selecting the "best" ML algorithm based on the provided performance metric
  rmse_min <- sapply(m_list, FUN = function(x) min(x$results[, metric]), simplify = T)
  ind <- which(rmse_min == min(rmse_min))

  ### Returning result
  return(list(best = m_list[[ind]], best.metric = min(rmse_min), all_methods = m_list))

}

#' Panel Cross Validation Staggered
#'
#' Panel Cross Validation in a staggered adoption setting, i.e., when there are groups of units
#' treated at different times.
#'
#' @param data A 'PanelMLCM' object from a previous call to \code{as.PanelMLCM}.
#' @param pcv_block Number of pre-intervention times to block for panel cross validation. Defaults to 1, see Details.
#' @param metric Character, the performance metric that should be used to select the optimal model.
#'               Possible choices are either \code{"RMSE"} (the default) or \code{"Rsquared"}.
#' @param default_par List of parameters for the default ML algorithms. See the Details of the function
#'                    \code{PanelCrossValidation}
#' @param trControl Optional, used to customize the training step. It must be the output from a call to \code{trainControl} from the \code{caret} package.
#' @param ML_methods Optional list of ML methods to be used as alternatives to the default methods. Each method must be supplied
#'                   as a named list of two elements: a character defining the method name from all the ones available in \code{caret}
#'                   and the grid of parameter values to tune via the panel cross validation. See Details and the examples for additional explanations.
#' @return For each intervention group, a list containing the best-performing ML method, the value of the given \code{metric} and the detailed results
#'         of all methods tried during the PCV routine (see also \code{PanelCrossValidation}).
#' @details
#' The function works in the same way as \code{PanelCrossValidation} as it repeats independently the PCV procedure for each group of unit defined by
#' their treatment date. Note that the \code{int_date} argument must be a column of \code{data}.

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
#' ### Example 1. Using the default ML methods in a staggered setting
#' set.seed(1)
#'
#' # Assume the following intervention dates
#' int_year_i <- c(rep(2019, times = 60), rep(2020, times = 40))
#' int_year <- rep(int_year_i, each = length(unique(data$year)))
#'
#' # Define data_stag
#' data_stag <- data.frame(int_year, data)
#'
#' # Organizing the dataset with as.PanelMLCM
#' newdata <- as.PanelMLCM(y = data_stag[, "Y"], timevar = data_stag[, "year"], id = data_stag[, "ID"],
#'                         int_date = data_stag[, "int_year"],
#'                         x = data_stag[, !(names(data_stag) %in% c("Y", "ID", "year", "int_year"))], y.lag = 2)
#'
#' # Panel Cross Validation in a staggered setting with different ML methods
#' pcv_multi <- PanelCrossValidationStag(data = newdata, ML_methods = list(enet, linreg, pls))
#'
#' # ATE estimation
#' causal <- StagMLCM(data = newdata, int_date = "int_date", inf_type = "block", nboot = 10, PCV = pcv_multi)
#'
PanelCrossValidationStag <- function(data, pcv_block = 1, metric = "RMSE",
                                     default_par = list(gbm = list(depth = c(1,2,3),
                                                                    n.trees = c(500, 1000),
                                                                    shrinkage = seq(0.01, 0.1, by = 0.02),
                                                                    n.minobsinnode = c(10,20)),
                                                         rf = list(ntree = 500),
                                                         pls = list(ncomp = c(1:5)),
                                                         lasso = list(fraction = seq(0.1, 0.9, by = 0.1))),
                                      trControl = NULL, ML_methods = NULL){

  ### Parameter checks
  if(!any(class(data) %in% "PanelMLCM")) stop("Invalid class in the PanelCrossValidation function, something is wrong with as.PanelMLCM")
  if(any(sapply(unique(data[, "int_date"]), FUN = function(x)(sum(unique(data[, "Time"]) < x) - pcv_block < 1)))) stop("Panel cross validation must be performed in at least one time period")
  if(pcv_block <= 0) stop("The number of 'pcv_block' time periods for panel cross validation must be at least 1")
  if(!metric %in% c("RMSE", "Rsquared")) stop("Metric not allowed, check documentation")
  if(!is.null(ML_methods)){if(!is.list(ML_methods)) stop("ML_methods must be a list")}
  if(is.list(ML_methods)){
    if(any(sapply(ML_methods, FUN = function(x)(!all(c("method", "tuneGrid") %in% names(x)))))) stop("Each method in 'ML_methods' must be a named list, check documentation")}
  if(is.null(ML_methods) & !is.list(default_par)) stop(" 'default_par' must be a list, check the documentation")
  if(is.null(ML_methods) & any(!names(default_par) %in% c("gbm", "rf", "lasso", "pls"))) stop("'default_par' must be a list of parameters for the default gbm, rf, lasso and pls algorithm, check the documentation")

  ### Applying PCV for each intervention group
  int_date <- data[, "int_date"]
  data[, "int_date"] <- NULL
  res <- mapply(unique(int_date), FUN = function(x){datax <- data[int_date == x,];
                                                    PanelCrossValidation(data = datax, int_date = x, pcv_block = pcv_block, metric = metric, default_par = default_par,
                                                                         trControl = trControl, ML_methods = ML_methods)}, SIMPLIFY = FALSE)
  ### Returning results
  names(res) <- paste0("int_", unique(int_date))
  return(res)

}
