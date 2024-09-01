###############################################################################
##
## MLCM functions
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: August 2024
##
## Functions: MLCM, MLCMStag, as.PanelMLCM
## Internal : .ate_est, .cate_est, .true_lag, .groupavg_eff, .check_MLCM, .check_xcate
##
###############################################################################

#' Machine Learning Control Method
#'
#' This function is the main workhorse of the package 'MachineControl'. It takes as
#' input a panel dataset, i.e., multiple units observed at several points in time and
#' exposed simultaneously to some policy (indicated by \code{int_date}). It then performs
#' the Panel Cross Validation (PCV) by comparing the predictive performance of several Machine
#' Learning (ML) methods in the pre-intervention periods and then outputs the estimated
#' Average Treatment Effect (ATE) of the policy and its confidence interval estimated by
#' bootstrap. Details on the causal assumptions, the estimation process and inference
#' can be found in Cerqua A., Letta M., and Menchetti F. (2023). The method is especially
#' suited when there are no control units available and it can also accomodate
#' staggered adoption settings where groups of units are treated or shocked at
#' different times. See Details.
#'
#'
#' @param data        A panel dataset in long form, having one column for the time variable, one column for the units'
#'                    unique IDs, one column for the outcome variable and one or more columns for the covariates.
#' @param int_date    The date of the intervention, treatment, policy introduction or shock. It must be contained in 'timevar'.
#'                    By default, this is the first period that the causal effect should be computed for.
#' @param inf_type    Character, type of inference to be performed. Possible choices are 'classic', 'block', 'bc classic', 'bc block', 'bca'. Defaults to 'block'.
#' @param y           Character, name of the column containing the outcome variable. It can be omitted for \code{PanelMLCM} objects.
#' @param timevar     Character, name of the column containing the time variable. It can be omitted for \code{PanelMLCM} objects.
#' @param id          Character, name of the column containing the ID's. It can be omitted for \code{PanelMLCM} objects.
#' @param y.lag       Optional, number of lags of the dependent variable to include in the model. Defaults to zero.
#' @param nboot       Number of bootstrap replications, defaults to 1000.
#' @param pcv_block   Number of pre-intervention times to block for panel cross validation. Defaults to 1, see Details.
#' @param metric      Character, the performance metric that should be used to select the optimal model.
#'                    Possible choices are either \code{"RMSE"} (the default) or \code{"Rsquared"}.
#' @param default_par List of parameters for the default ML algorithms. See the Details of \code{PanelCrossValidation}.
#' @param PCV         Optional, list returned from a previous call to \code{PanelCrossValidation} or \code{PanelCrossValidationStag}.
#' @param CATE        Whether the function should estimate also CATE (defaults to \code{FALSE}). Currently not allowed for staggered adoption. See Details.
#' @param x.cate      Optional matrix or data.frame of external regressors to use as predictors of CATE. If missing, the
#'                    same covariates used to estimate ATE will be used. See Details.
#' @param alpha       Confidence interval level to report for the ATE. Defaulting to 0.05 for a two sided
#'                    95\% confidence interval.
#' @param fe          Logical, whether to include fixed effects dummy variables. Defaults to false.
#'
#' @details
#' The panel \code{data} must at least include the response variable, a column of the time variable,
#' and a column with the unit identifiers. Lagged values of the outcome, lagged or contemporaneous
#' exogenous covariates, control series can also be added in the panel dataset. The ML algorithms
#' will automatically treat every column that is left as additional information to improve the
#' prediction of the counterfactual post-intervention outcome. If \code{fe = TRUE}, N-1 dummy variables
#' will be added to the dataset, where N is the total number of units in the panel.
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
#'   \item \code{int_date}: the intervention date
#'   \item \code{ate}: a list with the estimated ATE and its confidence interval
#'   \item \code{individual}: a list with the estimated individual effects and their confidence intervals
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
#' ### Example 1. Estimating ATE (with default ML methods)
#' \dontrun{
#' # Estimation
#' fit <- MLCM(data = data, y = "Y", timevar = "year", id = "ID", int_date = 2019,
#'             inf_type = "classic", nboot = 10, y.lag = 2)
#'
#' # ATE
#' fit$ate$estimate
#' fit$ate$conf.interval
#'
#' # Individual effects
#' head(fit$individual$estimate)
#' head(fit$individual$conf.interval)
#' }
#' ### Example 2. Estimating ATE and CATE (with external regressors in CATE)
#' \dontrun{
#' # Simulating time-varying external regressors
#' x.cate <- cbind(ID = rep(1:100, each = 2), year = rep(2019:2020, times = 100), x1 = rnorm(200),
#'                 x2 = 2*rnorm(200), x3 = sample(1:5, size = 200, replace = TRUE))
#'
#' # Estimation
#' fit <- MLCM(data = data, y = "Y", timevar = "year", id = "ID", int_date = 2019,
#'             inf_type = "classic", nboot = 10, CATE = TRUE, y.lag = 2, x.cate = x.cate)
#'
#' # CATE
#' plot(fit, type = "cate")
#' }
#'
#' # Example 3. Estimating ATE with custom ML methods
#'
#' # Organizing the dataset
#' newdata <- as.PanelMLCM(y = data[, "Y"], timevar = data[, "year"], id = data[, "ID"],
#'                         x = data[, !(names(data) %in% c("Y", "ID", "year"))],
#'                         y.lag = 2)
#'
#' # PCV with enet, linreg and pls
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
#' pls <- list(method = "pls", tuneGrid = expand.grid(ncomp = c(1:5)))
#'
#' pcv <- PanelCrossValidation(data = newdata, int_date = 2019,
#'                             ML_methods = list(enet, linreg, pls))
#'
#' # ATE estimation and inference
#' causal <- MLCM(data = newdata, int_date = 2019, inf_type = "classic", PCV = pcv,
#'                nboot = 10, CATE = FALSE, y.lag = 2)
#' plot(causal, type = "ate")
#'
MLCM <- function(data, int_date, inf_type = "block", y = NULL, timevar = NULL, id = NULL, y.lag = 0, nboot = 1000, pcv_block = 1, metric = "RMSE",
                 default_par = list(gbm = list(depth = c(1,2,3),
                                               n.trees = c(500, 1000),
                                               shrinkage = seq(0.01, 0.1, by = 0.02),
                                               n.minobsinnode = c(10,20)),
                                    rf = list(ntree = 500),
                                    pls = list(ncomp = c(1:5)),
                                    lasso = list(fraction = seq(0.1, 0.9, by = 0.1))),
                 PCV = NULL, CATE = FALSE, x.cate = NULL, alpha = 0.05, fe = FALSE){

  ### Parameter checks

  ## ATE checks
  .check_MLCM(data = data, int_date = int_date, inf_type = inf_type, y = y, timevar = timevar, id = id, y.lag = y.lag, nboot = nboot, pcv_block = pcv_block, metric = metric, default_par = default_par, PCV = PCV, alpha = alpha, fe = fe)

  ## CATE checks
  if(CATE & is.character(int_date)) stop("CATE = TRUE non allowed in staggered settings")
  if(CATE & !is.null(x.cate)){

    x.cate <- .check_xcate(x.cate = x.cate, data = data, id = id, timevar = timevar)

  } else if (!is.null(x.cate) & !CATE){ stop("Inserted external data for CATE estimation but 'CATE' is set to FALSE")}

  ## 1. Structuring the panel dataset in the required format
  if("PanelMLCM" %in% class(data)){

    data_panel <- data

  } else {

    data_panel <- as.PanelMLCM(y = data[, y], timevar = data[, timevar], id = data[, id],
                               x = data[, !(names(data) %in% c(y, id, timevar))], y.lag = y.lag, fe = fe)

  }

  ## 2. Panel cross-validation
  if(is.null(PCV)){

  best <- PanelCrossValidation(data = data_panel, int_date = int_date, pcv_block = pcv_block, metric = metric)$best

  } else {

    best <- PCV$best

  }

  ## 3. ATE & individual effects estimation
  effects <- .ate_est(data = data_panel, int_date = int_date, best = best, metric = metric, y.lag = y.lag, ran.err = FALSE)
  ate <- effects$ate
  ind_effects <- effects$ind_effects
  placebo_effects <- effects$placebo_effects

  ## 4. ATE & individual effects inference
  invisible(capture.output(

  boot_inf <- .boot_ate(data = data_panel, int_date = int_date, bestt = best, type = inf_type, nboot = nboot, ate = ate,
                           alpha = alpha, metric = metric, y.lag = y.lag, ind.eff = ind_effects, fe = fe, placebo.eff = placebo_effects)

  ))

  ## 5. CATE estimation & inference
  if(CATE){

    cate_effects <- .cate_est(data = data_panel, int_date = int_date, ind_effects = ind_effects, x.cate = x.cate, nboot = nboot, alpha = alpha)
    cate <- cate_effects$cate
    cate.inf <- cate_effects$cate.inf

  } else {

    cate <- NULL
    cate.inf <- NULL
  }

  ## 6. Saving results
  res <- list(best_method = best,
              int_date = int_date,
              ate = list(estimate = c(placebo_effects, ate), conf.interval = cbind(boot_inf$conf.placebo, boot_inf$conf.ate),
                         ate.boot = boot_inf$ate_boot),
              individual = list(estimate = ind_effects, conf.interval = boot_inf$conf.individual),
              cate = cate, cate.inf = cate.inf)
  class(res) <- "MLCM"
  return(res)

}

#' Machine Learning Control Method with staggered adoption
#'
#' This function adapts the MLCM approach to a staggered adoption setting. It takes as
#' input a panel dataset where groups of units are treated or shocked at
#' different times. It then performs the Panel Cross Validation (PCV) by comparing
#' the predictive performance of several Machine Learning (ML) methods in the pre-intervention
#' periods and then outputs the estimated effect of the policy. Three main causal estimands
#' are computed (see Details) and their confidence intervals are estimated by
#' bootstrap.
#'
#'
#' @param data        A panel dataset in long form, having one column for the time variable, one column for the units'
#'                    unique IDs, one column for the outcome variable and one or more columns for the covariates. For staggered
#'                    adoption settings, the dataset must also have a column for the intervention dates.
#' @param int_date    Name of the column containing the dates of the interventions.
#' @param inf_type    Character, type of inference to be performed. Possible choices are 'classic', 'block', 'bc classic', 'bc block', 'bca'. Defaults to 'block'.
#' @param y           Character, name of the column containing the outcome variable. It can be omitted for \code{PanelMLCM} objects.
#' @param timevar     Character, name of the column containing the time variable. It can be omitted for \code{PanelMLCM} objects.
#' @param id          Character, name of the column containing the ID's. It can be omitted for \code{PanelMLCM} objects.
#' @param y.lag       Optional, number of lags of the dependent variable to include in the model. Defaults to zero.
#' @param nboot       Number of bootstrap replications, defaults to 1000.
#' @param pcv_block   Number of pre-intervention times to block for panel cross validation. Defaults to 1, see Details.
#' @param metric      Character, the performance metric that should be used to select the optimal model.
#'                    Possible choices are either \code{"RMSE"} (the default) or \code{"Rsquared"}.
#' @param default_par List of parameters for the default ML algorithms. See the Details of \code{PanelCrossValidation}.
#' @param PCV         Optional, list returned from a previous call to \code{PanelCrossValidation} or \code{PanelCrossValidationStag}.
#' @param alpha       Confidence interval level to report for the ATE. Defaulting to 0.05 for a two sided
#'                    95\% confidence interval.
#' @param fe          Logical, whether to include fixed effects dummy variables. Defaults to false.
#'
#' @details
#'
#' This function extends the MLCM method in a staggered adoption settings where groups of units are treated or shocked at
#' different times, e.g., one group of units becomes treated starting from 2022, another group of unit becomes treated
#' in 2024 and so on. In this case, the \code{data} must include a column with the units' intervention times.
#'
#' The target causal estimand is the global ATE estimated as follows: a counterfactual is forecasted for
#' the units in each treatment group or cohort. If those groups are treated for more than one time period, a temporal average ATE
#' is calculated for each group and these are averaged together (a weighted average based on group size is used).
#' The result is the global effect of the policy or shock on all units that were treated.
#'
#' In addition, the function computes an ATE for each cohort as well as a group-average effect of the policy at each point
#' in time following the intervention and the ATEs are averaged across the cohorts by their size (note that if the first cohort
#' becomed treated in 2022 and the others start being treated from 2023, the group-average ATE in 2022 is just the ATE of the
#' first cohort).

#'
#' @return A list with the following components:
#' \itemize{
#'   \item \code{best_method}: the best-performing ML algorithm as selected by the PCV routine
#'   \item \code{int_date}: a vector with the intervention dates
#'   \item \code{global_ate}: a list with the estimated global ATE and its confidence interval
#'   \item \code{tempavg_individual}: a list with the estimated individual effects (averaged over
#'                                    the post-intervention periods) and their confidence intervals
#'   \item \code{group_ate}: a list with the estimated ATE for each treatment cohort and their
#'                           confidence intervals
#'   \item \code{groupavg}: a list with the estimated ATE for each period following the intervention and
#'                          their confidence intervals
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
#' ### Example 1. Estimating ATE in staggered settings (with custom ML methods)
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
#'                         x = data_stag[, !(names(data_stag) %in% c("Y", "ID", "year", "int_year"))],
#'                         y.lag = 2)
#'
#' # Panel Cross Validation in a staggered setting with different ML methods
#' enet <- list(method = "enet",
#'              tuneGrid = expand.grid(
#'                fraction = seq(0.1, 0.9, by = 0.1),
#'                lambda = seq(0.1, 0.9, by = 0.1)))
#'
#' linreg <- list(method = "lm",
#'                tuneGrid = expand.grid(
#'                  intercept = seq(0, 10, by = 0.5)))
#'
#' pls <- list(method = "pls", tuneGrid = expand.grid(ncomp = c(1:5)))
#' pcv_stag <- PanelCrossValidationStag(data = newdata, ML_methods = list(enet, linreg, pls))
#'
#' # Causal effect estimation and inference
#' causal_stag <- MLCMStag(data = newdata, int_date = "int_date", inf_type = "block", nboot = 10, PCV = pcv_stag, y.lag = 2)
#'
#' # Plotting
#' g <- plot(causal_stag, type = "cohort")
#' g[[1]]
#' g[[2]]
#'
MLCMStag <- function(data, int_date, inf_type = "block", y = NULL, timevar = NULL, id = NULL, y.lag = 0, nboot = 1000, pcv_block = 1, metric = "RMSE",
                     default_par = list(gbm = list(depth = c(1,2,3),
                                                   n.trees = c(500, 1000),
                                                   shrinkage = seq(0.01, 0.1, by = 0.02),
                                                   n.minobsinnode = c(10,20)),
                                        rf = list(ntree = 500),
                                        pls = list(ncomp = c(1:5)),
                                        lasso = list(fraction = seq(0.1, 0.9, by = 0.1))),
                 PCV = NULL, alpha = 0.05, fe = FALSE){

  ### Parameter checks

  ## ATE checks
  .check_MLCM(data = data, int_date = int_date, inf_type = inf_type, y = y, timevar = timevar, id = id, y.lag = y.lag, nboot = nboot, pcv_block = pcv_block, metric = metric, default_par = default_par, PCV = PCV, alpha = alpha, fe = fe)

  ## 1. Structuring the panel dataset in the required format
  if("PanelMLCM" %in% class(data)){

     data_panel <- data

  } else {

     data_panel <- as.PanelMLCM(y = data[, y], timevar = data[, timevar], id = data[, id],
                                 int_date = data[, int_date],
                                 x = data[, !(names(data) %in% c(y, id, timevar, int_date))],
                                 y.lag = y.lag, fe = fe)
  }

  ## 2. Panel cross-validation in staggered settings
  if(is.null(PCV)){

    best <- lapply(PanelCrossValidationStag(data = data_panel, pcv_block = pcv_block, metric = metric, default_par = default_par),
                   FUN = function(x)(x[["best"]]))

  } else {

    best <- lapply(PCV, FUN = function(x)(x[["best"]]))

  }

  ## 3. Global ATE & individual effects estimation
  #  3.0. m-applying 'ate_est' to each intervention group
  nint <- unique(data_panel$int_date)
  ate_i <- mapply(x = nint, y = best, FUN = function(x,y){datax <- data_panel[data_panel$int_date == x,]
                                                          datax$int_date <- NULL
                                                          .ate_est(data = datax, int_date = x, best = y, metric = metric, ran.err = FALSE, y.lag = y.lag)}, SIMPLIFY = FALSE)
  names(ate_i) <- paste0("int_", nint)

  #  3.1. Global ATE: estimation
  #  This is a temporal average ATE -- avg. of individual effects over both units and times and then averaged across cohorts
  weights <- sapply(unique(data_panel$int_date), FUN = function(x)(sum(data_panel$int_date == x)))/nrow(data_panel)
  global_ate <- weighted.mean(sapply(ate_i, FUN = function(x)(mean(x$ate))), w = weights)

  #  3.2. Global individual effects: estimation
  #  This is a temporal average effect for each unit -- avg. of individual effects over time
  global_ind <- lapply(ate_i, FUN = function(x){temp.avg <- rowMeans(as.matrix(x$ind_effects[, -1]))
  data.frame(ID = x$ind_effects[, 1], temp.avg = temp.avg)})
  global_ind <- do.call("rbind", global_ind)
  ta_ind_effects <- merge(global_ind, unique(data_panel[, c("ID", "int_date")]))

  #  3.3. Group-average effects: estimation
  #  This is an average of individual effects in the pre- and post-intervention across cohorts
  #  Note that the pre-intervention effects are 'placebo' effects
  group <- .groupavg_eff(data = data_panel, alpha = alpha, ate_i = ate_i)

  ## 4. Global ATE & individual effects inference
  #  4.0. m-applying 'boot_ate' to each intervention group
  ind <- which(colnames(ta_ind_effects) == "int_date")
  boot_inf <- mapply(x = nint, y = best, FUN = function(x,y){datax <- data_panel[data_panel$int_date == x,]
                                                             datax$int_date <- NULL
                                                             .boot_ate(data = datax, int_date = x, bestt = y, type = inf_type, nboot = nboot,
                                                             ate = global_ate, ind.eff = ta_ind_effects[, -ind], alpha = alpha, metric = metric, y.lag = y.lag, fe = fe)}, SIMPLIFY = FALSE)

  names(boot_inf) <- paste0("int_", nint)

  #  4.1. Global ATE: confidence interval
  global_ate_boot <- apply(sapply(boot_inf, FUN = function(x)(colMeans(x$ate_boot))), 1,
                           FUN = weighted.mean, w = weights)
  conf.global.ate <- quantile(global_ate_boot, probs = c(alpha/2, 1- alpha/2))

  # 4.2. Global individual effects: confidence interval
  tempavg_ind <- lapply(boot_inf, FUN = function(x){n <- NROW(x$ate_boot)
                                                    asvec <- unlist(asplit(x$ind_boot, 1))
                                                    asarr <- array(asvec, dim = c(nboot, nrow(x$ind_boot)/n ,n))
                                                    apply(asarr, c(1,2), mean)})
  conf.tempavg.ind <- lapply(tempavg_ind, FUN = function(x)(apply(x, 2, quantile, probs = c(alpha/2, 1- alpha/2))))
  conf.tempavg.ind <- t(do.call("cbind", conf.tempavg.ind))
  conf.tempavg.ind <- data.frame(do.call("rbind", sapply(nint, FUN = function(x)(unique(data_panel[data_panel$int_date == x, c("ID", "int_date")])), simplify = FALSE)),
                                  conf.tempavg.ind)

  # 4.3. Group-average effects: confidence interval
  lower <- sapply(boot_inf, FUN = function(x)(apply(rbind(x$placebo_boot, x$ate), 1, FUN = quantile, probs = alpha/2)))
  lower <- lower[, order(colnames(lower))]
  group.lower <- sapply(1:nrow(lower), FUN = function(i)(weighted.mean(lower[i,], w = group$global.weights[i,])))
  upper <- sapply(boot_inf, FUN = function(x)(apply(rbind(x$placebo_boot, x$ate), 1, FUN = quantile, probs = 1- alpha/2)))
  upper <- upper[, order(colnames(upper))]
  group.upper <- sapply(1:nrow(upper), FUN = function(i)(weighted.mean(upper[i,], w = group$global.weights[i,])))
  conf.groupavg <- rbind(group.lower, group.upper)
  colnames(conf.groupavg) <- names(group$groupavg_eff)
  rownames(conf.groupavg) <- c(alpha/2, 1- alpha/2)

  ## 5. Saving results
  cohort_eff <- lapply(1:length(nint), function(i)(list(estimate = c(ate_i[[i]]$placebo_effects, ate_i[[i]]$ate),
                                          conf.interval = cbind(boot_inf[[i]]$conf.placebo, boot_inf[[i]]$conf.ate))))
  names(cohort_eff) <- names(boot_inf)
  Npre <- length(unique(data_panel[data_panel$Time < min(nint), "Time"]))

  res <- list(best_method = best, int_date = unique(data_panel$int_date),
              global_ate = list(estimate = c(group$groupavg_eff[1:Npre], global_ate = global_ate),
                                conf.interval = cbind(conf.groupavg[,1:Npre], conf.global.ate), boot = global_ate_boot),
              tempavg_individual = list(estimate = ta_ind_effects, conf.interval = conf.tempavg.ind),
              group_ate = cohort_eff,
              groupavg = list(estimate = group$groupavg_eff, conf.interval = conf.groupavg))

  class(res) <- "MLCMStag"
  return(res)

}

#' Structuring the panel dataset
#'
#' This function takes as input the panel dataset given by the user and changes the
#' ordering and the names of the columns to obtain an object of class 'PanelMLCM'
#' to be used by the function 'MLCM'.
#'
#' @param y Numeric, the outcome variable.
#' @param x Matrix or data.frame of covariates to include in the model.
#' @param timevar  The column containing the time variable. It can be \code{numeric}, \code{integer} or
#'                 \code{Date}.
#' @param id Numeric, the column containing the ID's.
#' @param int_date Optional, only needed in staggered adoption setting to identify the column
#'                 reporting the dates of intervention.
#' @param y.lag Optional, number of lags of the dependent variable to include in the model.
#' @param fe Logical, whether to include fixed effects dummy variables. Defaults to false.
#'
#' @details This function is mainly for internal use of \code{MLCM}. It is exported to give
#' users full flexibility in the choice of the ML algorithms and during the panel cross validation.
#' See the documentation of the function \code{PanelCrossValidation} or, for a staggered adoption
#' setting, \code{PanelCrossValidationStag}.
#'
#' If \code{fe = TRUE}, N-1 dummy variables will be added to the dataset, where N is the total number
#' of units in the panel).
#'
#' @return An object of class 'data.frame' and 'PanelMLCM'.
#' @export
#' @examples
#'
#' # Example 1. Start from a disorganized panel dataset
#' head(data)
#'
#' # Run as.PanelMLCM
#' newdata <- as.PanelMLCM(y = data[, "Y"], timevar = data[, "year"], id = data[, "ID"],
#'                         x = data[, !(names(data) %in% c("Y", "ID", "year"))], y.lag = 2)
#'
#' # Results
#' head(newdata)
#'
#' # Example 2. Disorganized panel dataset in a staggered adoption setting
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
#'                        int_date = data_stag[, "int_year"],
#'                        x = data_stag[, !(names(data_stag) %in% c("Y", "ID", "year", "int_year"))], y.lag = 2)
#' # Results
#' head(newdata)
#'
as.PanelMLCM <- function(y, x, timevar, id, int_date = NULL, y.lag = 0, fe = FALSE){

  # Parameter checks
  if(!(is.numeric(y) & length(y)>1)) stop("y must be a numeric vector of length greater than 1")
  if(!any(class(x) %in% c("numeric", "matrix", "data.frame"))) stop("x must be a vector, matrix or data.frame")
  if(NROW(x) != length(y)) stop("NROW(x) != length(y)")
  if(!any(class(timevar) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt", "integer", "numeric")) | length(timevar) != length(y)) stop("timevar must be a numeric vector or a 'Date' object of the same length as y")
  if(!(is.numeric(id) & length(id) == length(y))) stop("id must be a numeric vector of the same length as y")
  if(!is.null(int_date) & length(unique(id))*length(unique(timevar)) != length(y)) warning("The panel is unbalanced")
  if(!is.null(int_date) & length(int_date) != length(y)) stop("int_date must be NULL or, for staggered adoption, a vector of the same length as y ")
  if(any(!(unique(int_date) %in% timevar))) stop("all dates in 'int_date' must be contained in timevar")
  if(!is.numeric(y.lag) | y.lag < 0 ) stop("y.lag must be numeric and strictly positive")
  if(length(unique(timevar)) <= y.lag) stop("The number of selected lags is greater or equal to the number of times, resulting in an empty dataset")
  if(!is.logical(fe)) stop("fe must be logical, defaults to FALSE")

  ### STEP 1. Structuring the panel dataset
  if(is.null(int_date)){

    panel <- data.frame(Time = timevar, ID = id, Y = y, x)

  } else {

    panel <- data.frame(Time = timevar, ID = id, Y = y, int_date = int_date, x)
    #panel <- panel[order(panel$int_date), ]

  }

  if(fe){

    # Generating matrix of FE dummies
    dm <- model.matrix(~factor(panel$ID) -1)
    dm <- dm[,-1]
    colnames(dm) <- paste0("FE_", unique(panel$ID)[-1])
    panel <- data.frame(panel, dm)
  }

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

.ate_est <- function(data, int_date, best, metric, ran.err, y.lag){

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

    ### Step 4. ATE estimation (observed - predicted). Note that when there is more than 1 post-intervention
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
      minl <- min(y.lag, length(unique(postimes))-i)

      for(l  in 1:minl){

        data[(post+l), paste0("Ylag",l)] <- pred - error
        data

      }
    }
  }

  ### Step 5. Estimation of placebo effects (effects in the pre-intervention period)
  ind <- data$Time < int_date
  placebo <- data[ind, "Y"] - predict.train(fit, newdata = data[ind, ])
  npre <- unique(data[ind, "Time"])
  placebo_mean <- sapply(npre, FUN = function(x)(mean(placebo[data[ind, "Time"] == x])))
  names(placebo_mean) <- npre
  placebo_mean <- placebo_mean[order(names(placebo_mean))]

  ### Step 6. Returning the matrix of individual effects, the ATE and the placebo effects
  ind_effects <- cbind(ID = unique(data$ID), ind_effects)
  ind_effects <- ind_effects[order(ind_effects[, "ID"], decreasing = F), ]
  ate <- colMeans(as.matrix(ind_effects[, -1]))
  names(ate) <- colnames(ind_effects)[-1]

  return(list(ind_effects = ind_effects, ate = ate,
              placebo_effects = placebo_mean))
}

#' CATE estimation
#'
#' @param data A 'PanelMLCM' object from a previous call to \code{as.PanelMLCM}.
#' @param int_date The date of the intervention, treatment, policy introduction or shock.
#' @param ind_effects A matrix of estimated individual causal effects, returned from a previous
#'                    call to \code{ate_est}.
#' @param x.cate Optional, a matrix or data.frame of external regressors to use for CATE estimation.
#' @param nboot Number of bootstrap iterations.
#' @param alpha Confidence interval level to report for the ATE. Defaulting to 0.05 for a two sided
#'              95\% confidence interval.
#'
#' @return A list with the following components:
#' \itemize{
#' \item cate: a list with as many components as the post-intervention period,
#'       each containing an 'rpart' object.
#' \item cate.inf: a list with as many components as the post-intervention period,
#'       each containing the estimated CATE, its variance and confidence interval
#'       for the terminal nodes.
#' }
#' @noRd

.cate_est <- function(data, int_date, ind_effects, x.cate, nboot, alpha){

  ### Step 1. Selecting post-intervention times
  post <- which(data$Time >= int_date)
  postimes <- data$Time[post]

  ### Step 2. Matrix containing the estimated individual effects and post-intervention covariates
  if(is.null(x.cate)){

    data_cate <- data.frame(Time = postimes, effect = c(t(ind_effects[,-1])), data[post, !names(data) %in% c("Y","Time","ID", "Ylag1")])

  } else {

    data_cate <- data.frame(Time = postimes, ID = data[post, "ID"], effect = c(t(ind_effects[,-1])))
    data_cate <- merge(data_cate, x.cate, by = c("ID", "Time"))
    data_cate <- data_cate[order(data_cate[, "ID"], decreasing = FALSE),]
    data_cate$ID <- NULL

  }

  ### Step 3. CATE estimation & inference
  cate <- lapply(unique(postimes), FUN = function(x){
    rpart(effect ~ ., method="anova", data = data_cate[data_cate$Time == x, -1], cp = 0,
          minbucket = 0.05*length(unique(data$ID)), model = TRUE)})
  mat <- data.frame(postimes, c(t(ind_effects[,-1])))
  cate.inf <- mapply(x = cate, y = unique(postimes), FUN = function(x,y)(
    .boot_cate(effect = mat[mat$postimes == y, -1], cate = x, nboot = nboot, alpha = alpha)), SIMPLIFY = FALSE)
  names(cate) <- names(cate.inf) <- unique(postimes)

  ### Step 4. Returning estimated CATE
  return(list(cate = cate, cate.inf = cate.inf))
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
#' A vector of the lagged variable having the same length as \code{x} (note that
#' there will be as many initial NAs as the number of \code{lag}.)
#'
#' @noRd
#' @importFrom utils head

.true_lag <- function(x, lag = 1){

  x_lag <- c(rep(NA, times = lag), head(x, n = length(x)-lag))
  return(x_lag)

}


#' Checking CATE
#'
#' Internal function, used when external regressors are included in the estimation of CATE.
#' See Details in MLCM function. The function is purely used to checks the concordance of
#' the information included in 'data' and 'x.cate' (e.g., both datasets should contain
#' the same unique identifiers).
#'
#' @param x.cate Matrix or data.frame of external regressors used for CATE estimation
#' @param data Matrix, data.frame or PanelMLCM object
#' @param id Character, variable in 'data' containing unique identifiers
#' @param timevar Character, variable in 'data' containing times
#'
#' @noRd

.check_xcate <- function(x.cate, data, id, timevar){

  if(!(is.matrix(x.cate)|is.data.frame(x.cate))) stop("x.var must be 'matrix' or 'data.frame'")
  x.cate <- as.data.frame(x.cate)

  if(is.null(id)){ # later, change with is.PanelMLCM

    ind <- which(colnames(data) == "ID")
    ti <- which(colnames(data) == "Time")
    if(!colnames(data)[ind] %in% colnames(x.cate)) stop("'x.cate' must have a column named 'ID' with unique identifiers, like the one in 'data'")
    if(!colnames(data)[ti] %in% colnames(x.cate)) stop("'x.cate' must have a column named 'Time' with time identifiers, like the one in 'data'")
    if(!all.equal(unique(data[, "ID"]), unique(x.cate[, "ID"]))) stop ("unique identifiers differ for 'x.cate' and 'data'")

  } else {

    if(!id %in% colnames(x.cate)) stop("'x.cate' does not have unique identifiers or the colnames of identifiers do not match those in 'data'")
    if(!timevar %in% colnames(x.cate)) stop("'x.cate' does not have time identifiers or the colnames do not match those in 'data'")
    if(!all.equal(unique(data[, id]),unique(x.cate[, id]))) stop ("unique identifiers differ for 'x.cate' and 'data'")
    colnames(x.cate)[which(colnames(x.cate) == paste(id))] <- "ID"
    colnames(x.cate)[which(colnames(x.cate) == paste(timevar))] <- "Time"
  }

  return(x.cate)
}

.check_MLCM <- function(data, int_date, inf_type, y , timevar, id, y.lag, nboot, pcv_block, metric, default_par, PCV, alpha, fe){

  ### Parameter checks

  # Checking class of 'data'
  if(!any(class(data) %in% c("matrix", "data.frame", "PanelMLCM"))) stop("data must be a matrix, a data.frame or a PanelMLCM object")
  if(!"PanelMLCM" %in% class(data) & any(c(is.null(y), is.null(timevar), is.null(id)))) stop("y or id or timevar are missing with no default")

  # Checking class of 'y', 'id' and 'timevar'
  ck <- mapply(list(y, timevar, id), FUN = function(par)(!(is.null(par) | class(par) == "character")))
  if(any(ck)) stop(paste(c("y", "timevar", "id")[which(ck)], "must be 'character' "))
  ck <- unlist(mapply(list(y, timevar, id), FUN = function(par)(!(is.null(par) | par %in% colnames(data)))))
  if(any(ck)) stop (paste("there is no column called", c(y, timevar, id)[which(ck)], "in 'data'"))

  # Checking 'y.lag', 'fe'
  if(!is.numeric(y.lag) | y.lag < 0 ) stop("y.lag must be numeric and strictly positive")  # should be integer
  if(!is.logical(fe)) stop("fe must be logical, defaults to FALSE")

  # Checking 'int_date'
  if(!any(class(int_date) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt", "numeric", "integer", "character"))) stop("int_date must be integer, numeric, date or character")

  if (is.character(int_date)) {
    time_col <- if (is.null(timevar)) "Time" else timevar
    if (any(!(unique(data[, int_date]) %in% data[, time_col]))) {
      stop(paste("all dates in 'int_date' must be contained in the", time_col, "column"))
    }
  } else {
    time_col <- if (is.null(timevar)) "Time" else timevar
    if (!int_date %in% data[, time_col]) {
      stop(paste("int_date must be contained in the", time_col, "column"))
    }
  }

  # Checking inf_type, nboot, metric, PCV, alpha
  if(!any(inf_type %in% c("classic", "block", "bc classic", "bc block", "bca"))) stop("Inference type not allowed, check the documentation")
  if(nboot < 1 | all(!class(nboot) %in% c("numeric", "integer")) | nboot%%1 != 0) stop("nboot must be an integer greater than 1")
  if(!metric %in% c("RMSE", "Rsquared")) stop("Metric not allowed, check documentation")
  if(!is.null(PCV) & !is.character(int_date) & (!"train" %in% class(PCV$best))) stop("Invalid PCV method, it should be a list returned from a previous call to 'PanelCrossValidation()'")
  if(!is.null(PCV) & is.character(int_date)){if(any(sapply(PCV, function(x)(!"train" %in% class(x$best))))) stop("Invalid PCV method, it should be a list returned from a previous call to 'PanelCrossValidationStag()' ")}
  if(alpha < 0 | alpha > 1) stop("Invalid confidence interval level, alpha must be positive and less than 1")

  # Checking default_par
  if(is.null(PCV) & !is.list(default_par)) stop(" 'default_par' must be a list, check the documentation")
  if(is.null(PCV) & any(!names(default_par) %in% c("gbm", "rf", "lasso", "pls"))) stop("'default_par' must be a list of parameters for the default gbm, rf, lasso and pls algorithm, check the documentation")

}

.groupavg_eff <- function(data, alpha, ate_i){

  # Internal function that computes the average causal effect for each time
  # pre- and post-intervention. The average is done across cohorts, i.e., for each cohort
  # a vector of effect is estimated with the function ate_est. Before the cohort gets treated,
  # the effect is a 'placebo' effect and after treatment it is the causal effect of the treatment received
  # The result of this function is a vector of effects and their confidence interval.

  ### Settings
  int.dates <- sort(unique(data$int_date))
  pre.int <- sort(unique(data$Time[data$Time < min(int.dates)]))
  post.int <- sort(unique(data$Time[data$Time>=min(int.dates)]))

  ### STEP 1. Weights pre-intervention
  weights <- table(data$int_date)/nrow(data)
  placebo.pre <- matrix(weights, length(pre.int), length(int.dates), byrow = T)
  colnames(placebo.pre) <- paste0("int_", names(weights))
  rownames(placebo.pre) <- pre.int

  ### STEP 2. Weights post-intervention
  placebo.post <- sapply(post.int, FUN = function(x){w <- table(data$int_date[data$int_date<=x])/sum(data$int_date<=x)
                                                     add <- int.dates[which(!int.dates %in% names(w))]
                                                     w <- c(w, setNames(rep(0, length(add)), add))
                                                     w <- w[order(names(w))]})
  colnames(placebo.post) <- post.int
  rownames(placebo.post) <- gsub(rownames(placebo.post), patter="^", replacement = "int_")


  ### STEP 3. Matrix of global weights (pre- and post-intervention)
  global.weights <- rbind(placebo.pre, t(placebo.post))

  ### STEP 4. Average effects at each time post-intervention -- weighted average based on cohort numerosity
  mat <- sapply(ate_i, FUN = function(x)(c(x$placebo_effects, x$ate)))
  rownames(mat) <- c(pre.int, post.int)
  mat <- mat[, order(colnames(mat))]
  groupavg_eff <- sapply(1:nrow(mat), FUN = function(i)(weighted.mean(mat[i,], w = global.weights[i,])))
  names(groupavg_eff) <- rownames(mat)

  ### STEP 5. Returning results
  return(list(groupavg_eff = groupavg_eff, global.weights = global.weights, time_eff = mat))

}
