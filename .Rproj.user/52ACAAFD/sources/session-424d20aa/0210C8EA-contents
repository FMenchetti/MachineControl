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
#' exposed simultaneously to some policy (indicated by int_date). It then performs
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
#' @param int_date    The date of the intervention, treatment, policy introduction or shock. It must be contained in 'timevar'.
#'                    By default, this is the first period that the causal effect should be computed for. See Details.
#' @param inf_type    Character, type of inference to be performed. Possible choices are 'classic', 'block', 'bc classic', 'bc block', 'bca'
#' @param y           Character, name of the column containing the outcome variable. It can be omitted for \code{PanelMLCM} objects.
#' @param timevar     Character, name of the column containing the time variable. It can be omitted for \code{PanelMLCM} objects.
#' @param id          Character, name of the column containing the ID's. It can be omitted for \code{PanelMLCM} objects.
#' @param y.lag       Optional, number of lags of the dependent variable to include in the model. Defaults to zero.
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
#' By default, the function estimates the ATE by averaging the individual causal effects across units. This is done
#' for each time period after \code{int_date} (the first ATE will be computed exactly at \code{int_date}).
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
#'   \item \code{conf.ate}: a (1-\code{alpha})% bootstrap confidence interval
#'   \item \code{cate}: a list of class \code{rpart} with the estimated regression-tree-based CATE
#'   \item \code{cate.inf}: a matrix containing the variance of CATE and the confidence interval
#'   estimated by bootstrap
#' }
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
#' ### Example 1. Estimating ATE and CATE (with default ML methods)
#'
#' # Estimation
#' fit <- MLCM(data = data, y = "Y", timevar = "year", id = "ID", int_date = 2019,
#'             inf_type = "classic", nboot = 10, CATE = TRUE, y.lag = 2)#'
#'
#' # ATE & CATE
#' fit$ate
#' # plot(fit$cate) ; text(fit$cate)
#'
#' # Bootstrap confidence interval
#' fit$conf.ate
#'
#' # Individual effects
#' head(fit$ind.effects)
#' head(fit$conf.individual)
#'

MLCM <- function(data, int_date, inf_type, y = NULL, timevar = NULL, id = NULL, y.lag = 0, nboot = 1000, pcv_block = 1, metric = "RMSE", PCV = NULL, CATE = FALSE, alpha = 0.05){

  ### Parameter checks
  if(!any(class(data) %in% c("matrix", "data.frame", "PanelMLCM"))) stop("data must be a matrix, a data.frame or a PanelMLCM object")
  if(!"PanelMLCM" %in% class(data) & any(c(is.null(y), is.null(timevar), is.null(id)))) stop("Unspecified columns in 'matrix' or 'data.frame' data")
  if(!(is.null(y) | class(y) == "character")) stop("y must be a character")
  if(!(is.null(timevar) | class(timevar) == "character")) stop("timevar must be a character")
  if(!(is.null(id) | class(id) == "character")) stop("id must be a character")
  if(!is.numeric(y.lag) | y.lag < 0 ) stop("y.lag must be numeric and strictly positive")  # should be integer
  if(!is.null(y)){if(!y %in% colnames(data)) stop (paste("there is no column called", y, "in 'data'"))}
  if(!is.null(timevar)){if(!timevar %in% colnames(data)) stop (paste("there is no column called", timevar, "in 'data'"))}
  if(!is.null(id)){if(!id %in% colnames(data)) stop (paste("there is no column called", id, "in 'data'"))}
  if(!any(class(int_date) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt", "numeric", "integer"))) stop("int_date must be integer, numeric or Date")
  if(is.null(timevar)){if(!int_date %in% data[, "Time"]) stop ("int_date must be contained in the 'Time' column")}
  if(!is.null(timevar)){if(!int_date %in% data[, timevar]) stop ("int_date must be contained in timevar")}
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
                               x = data[, !(names(data) %in% c(y, id, timevar))], y.lag = y.lag)

  }


  ### Panel cross-validation
  if(is.null(PCV)){

    best <- PanelCrossValidation(data = data_panel, int_date = int_date, pcv_block = pcv_block, metric = metric)$best

  } else {

    best <- PCV

  }

  ### Fit the best (optimized) ML algorithm on all pre-intervention data and make predictions in the post-intervention period
  # ind <- which(data_panel[, "Time"] < int_date)
  # set.seed(1)
  # invisible(capture.output(
  #  fit <- train(Y ~ .,
  #               data = data_panel[ind, !(names(data_panel) %in% c("ID", "Time"))],
  #               method = best$method,
  #               metric = metric,
  #               trControl = trainControl(method="none"),
  #               tuneGrid = best$bestTune)
  # ))

  ### ATE & individual effects estimation
  effects <- ate_est(data = data_panel, int_date = int_date, best = best, metric = metric, y.lag = y.lag, ran.err = FALSE)
  ate <- effects$ate
  ind_effects <- effects$ind_effects

  ### ATE & individual effects inference
  invisible(capture.output(

    boot_inf <- boot_ate(data = data_panel, int_date = int_date, bestt = best, type = inf_type, nboot = nboot, ate = ate,
                         alpha = alpha, metric = metric, y.lag = y.lag, ind.eff = ind_effects)

  ))

  ### CATE estimation & inference
  if(CATE){

    # Selecting post-intervention times
    post <- which(data_panel$Time >= int_date)
    postimes <- data_panel$Time[post]

    # Matrix containing the estimated individual effects and post-intervention covariates
    data_cate <- data.frame(Time = postimes, effect = c(t(ind_effects)), data_panel[post, !names(data_panel) %in% c("Y","Time","ID", "Ylag1")])

    # CATE estimation & inference
    cate <- lapply(unique(postimes), FUN = function(x){
      rpart(effect ~ ., method="anova", data = data_cate[data_cate$Time == x, -1], cp = 0, minbucket = 0.05*length(unique(data_panel$ID)))})
    mat <- data.frame(postimes, c(t(ind_effects)))
    cate.inf <- mapply(x = cate, y = unique(postimes), FUN = function(x,y)(
      boot_cate(effect = mat[mat$postimes == y, -1], cate = x, nboot = nboot, alpha = alpha)), SIMPLIFY = FALSE)
    names(cate.inf) <- unique(postimes)

  } else {

    cate <- NULL
    cate.inf <- NULL
  }


  ### Saving results
  return(list(best_method = best, fit = best, ate = ate, var.ate = boot_inf$var.ate, conf.ate = boot_inf$conf.ate,
              ind.effects = ind_effects, conf.individual = boot_inf$conf.individual, cate = cate, cate.inf = cate.inf))

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
#' @param y.lag Optional, number of lags of the dependent variable to include in the model.
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
#'                         x = data[, !(names(data) %in% c("Y", "ID", "year"))], y.lag = 2)
#'
#' # Results
#' head(newdata)

as.PanelMLCM <- function(y, x, timevar, id, y.lag = 0){

  # Parameter checks
  if(!(is.numeric(y) & length(y)>1)) stop("y must be a numeric vector of length greater than 1")
  if(!any(class(x) %in% c("numeric", "matrix", "data.frame"))) stop("x must be a vector, matrix or data.frame")
  if(NROW(x) != length(y)) stop("NROW(x) != length(y)")
  if(!any(class(timevar) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt", "integer", "numeric")) | length(timevar) != length(y)) stop("timevar must be a numeric vector or a 'Date' object of the same length as y")
  if(!(is.numeric(id) & length(id) == length(y))) stop("id must be a numeric vector of the same length as y")
  if(length(unique(id))*length(unique(timevar)) != length(y)) warning("The panel is unbalanced")
  if(!is.numeric(y.lag) | y.lag < 0 ) stop("y.lag must be numeric and strictly positive")
  if(length(unique(timevar)) <= y.lag) stop("The number of selected lags is greater or equal to the number of times, resulting in an empty dataset")

  ### STEP 1. Structuring the panel dataset
  panel <- data.frame(Time = timevar, ID = id, Y = y, x)

  ### STEP 2. Are there any past lags of 'y' to include?
  if(y.lag > 0){

    # Applying the internal '.true_lag' function to the Y variable of each unit in the panel
    ids <- unique(id)
    ylags <- sapply(1:y.lag, function(l){unlist(lapply(ids, function(x)(.true_lag(y[id == x], l))))})
    colnames(ylags) <- paste0("Ylag", 1:y.lag)
    panel <- data.frame(panel, ylags)

    # Removing initial NAs from the panel
    ind <- which(is.na(panel[, paste0("Ylag", y.lag)]))
    panel <- panel[-ind, ]

  }

  # Returning results
  class(panel) <- c("data.frame", "PanelMLCM")
  return(panel)
}

#' ATE estimation
#'
#' Internal function, used within the MLCM and bootstrap routines for the estimation of ATE
#'
#'
#' @param data A 'PanelMLCM' object from a previous call to \code{as.PanelMLCM}.
#' @param int_date The date of the intervention, treatment, policy introduction or shock.
#' @param best Object of class \code{train}, the best-performing ML method as selected
#'             by panel cross validation.
#' @param metric Character, the performance metric that should be used to select the optimal model.
#' @param ran.err Logical, whether to include a random error to the predicted/counterfactual
#'                post-intervention observations. It is set to \code{FALSE} when the objective is ATE
#'                estimation. It is set to \code{TRUE} when the objective is estimating bootstrap standard errors.
#'
#' @return A list with the following components:
#' \itemize{
#' \item ind_effects: matrix of estimated unit-level causal effects
#' \item ate: vector of estimated ATEs for each post-intervention period
#' }
#'
#' @noRd
#' @import caret
#' @importFrom stats rnorm
#' @importFrom stats sd

ate_est <- function(data, int_date, best, metric, ran.err, y.lag){

  ### Step 1. Settings (empty matrix)
  postimes <- data[which(data[, "Time"] >= int_date), "Time"]
  ind_effects <- matrix(NA, nrow = nrow(data[data$Time == int_date, ]) , ncol = length(unique(postimes)))
  colnames(ind_effects) <- unique(postimes)

  ### Step 2. Fit the best (optimized) ML algorithm on all pre-intervention data and make predictions
  ###         The following loops over the post-intervention periods and implements the recursive procedure
  ###         described in the paper
  for(i in 1:length(unique(postimes))){

    pre <- which(data[, "Time"] < postimes[i])
    post <- which(data[, "Time"] == postimes[i])
    invisible(capture.output(
      fit <- train(Y ~ .,
                   data = data[pre, !(names(data) %in% c("ID", "Time"))],
                   method = best$method,
                   metric = metric,
                   trControl = trainControl(method="none"),
                   tuneGrid = best$bestTune)
    ))

    ### Step 3. Counterfactual prediction, if the option 'ran.err' is active, a random error is added
    ###         to the prediction (recommended only during bootstrap to get reliable estimates of ATEs variance)
    if(ran.err){

      eps <- data[pre, "Y"] - predict.train(fit)
      error <- rnorm(n = nrow(data[post,]), mean = mean(eps), sd = sd(eps))
      pred <- predict.train(fit, newdata = data[post, ]) + error

    } else {

      pred <- predict.train(fit, newdata = data[post, ])
      error <- 0

    }

    ### STEP 4. ATE estimation (observed - predicted). Note that when there is more than 1 post-intervention
    ###         period and y.lag > 1, the MLCM routine will use the observed impacted series, disrupting all estimates.
    ###         e.g., int_date = 2020, 2 post-int periods, 2 lags: to predict Y_2020 MLCM will use Y_2018 (pre-int) and Y_2019 (pre-int),
    ###         but to predict Y_2021 MLCM will use Y_2019 (pre-int) and Y_2020 (post-int), which is not ok. With this last step,
    ###         we impute post-intervention Y's with their predicted counterfactual
    obs <- data[post, "Y"]
    ind_effects[,i] <- obs - pred

    if(length(unique(postimes)) > 1 & y.lag > 0 & i < length(unique(postimes))){

      # Substituting counterfactual Y (contemporaneous)
      data[post, "Y"] <- pred - error
      # Substituting counterfactual Y in future lags
      maxl <- max(1, y.lag-i+1)

      for(l  in 1:maxl){

        data[(post+l), paste0("Ylag",l)] <- pred - error
        data

      }
    }
  }

  ### Step 3. Returning the matrix of individual effects and the ATE
  return(list(ind_effects = ind_effects, ate = colMeans(ind_effects)))
}

#' Generation of lagged variables
#'
#' Internal function, used to generate lags of the dependent variable.
#'  This function shifts \code{x} backward of the given number of \code{lag}.
#' If \code{lag = 1} (the default), the function shifts the variable of 1 step.
#' For example, if \code{x[t]} is a time series, \code{lag = 1} generates \code{x[t-1]};
#' \code{lag = 2} generates \code{x[t-2]} and so on.
#'
#' @param x Numeric, variable to lag
#' @param lag Numeric, lag of \code{x} to generate, defaulting to 1.See Details.
#'
#' @return
#' A vector of the lagged variable having the same length as \code{x] (note that
#' there will be as many initial NAs as the number of \code{lag}.)
#'
#' @noRd
#' @importFrom utils head

.true_lag <- function(x, lag = 1){

  x_lag <- c(rep(NA, times = lag), head(x, n = length(x)-lag))
  return(x_lag)

}
