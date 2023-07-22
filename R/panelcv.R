###############################################################################
##
## Panel Cross Validation
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: July 2023
##
###############################################################################

PanelCrossValidation <- function(data, post_period, blocked = pcv_block){

  ### Parameter checks
  if(!any(class(data) %in% "PanelMLCM")){stop("Invalid class in the PanelCrossValidation function, something is wrong with as.PanelMLCM")}
  if(!(post_period %in% data[, "Time"]))(stop("post_period must be contained in timevar"))
  if(length(unique(data[, "Time"])) - 2 - blocked < 1) stop("Panel cross validation must be performed in at least one time period")
  if(blocked <= 0) stop("The number of 'blocked' time periods for panel cross validation must be at least 1")

  ### STEP 1. The CAST package is used to generate separate testing sets for each year

  Tt <- length(unique(data[, "Time"]))
  post_period <- post_period
  indices <- CreateSpacetimeFolds(data, timevar = "Time",
                                  k=Tt)
  trainx <- lapply(blocked:(Tt-2), FUN = function(x) unlist(indices$indexOut[1:x]))
  testx <- lapply(blocked:(Tt-2), FUN = function(x) unlist(indices$indexOut[[x+1]]))

  ### STEP 2. Set control function by specifying the training and testing folds that caret will use
  ###         for cross-validation and tuning of the hyperparameters (i.e., the combination of folds defined above)

  ctrl <- trainControl(index = trainx, indexOut = testx)

  ### STEP 3.  Tune the hyperparameters of each of the ML algorithms via temporal cross-validation

  # DOMANDA: di quanti dei successivi parametri vogliamo lasciare la scelta all'utente?
  # Vogliamo lasciare la scelta anche degli algoritmi?? Forse potremmo far sì che l'utente possa inserire le funzioni 'train' che vuole

  # Stochastic gradient boosting
  gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3),
                          n.trees=c(500, 1000, 1500, 2000),
                          # shrinkage = 0.1,
                          shrinkage = seq(0.01, 0.1, by = 0.01),
                          n.minobsinnode = c(10,20))

  set.seed(1)
  bo <- train(Y ~ .,
              data = data[, !(names(data) %in% c("ID", "Time"))],
              # data = data[, !(names(data) %in% c("Time"))],
              method = "gbm",
              metric = "RMSE",
              trControl = ctrl,
              tuneGrid = gbmGrid,
              verbose = FALSE)
  # Random forest
  set.seed(1)
  rf <- train(Y ~ .,
              data = data[, !(names(data) %in% c("ID", "Time"))],
              # data = data[, !(names(data) %in% c("Time"))],
              method = "rf",
              metric = "RMSE",
              search = "grid",
              trControl = ctrl,
              tuneGrid = expand.grid(mtry = (2:(ncol(data)-3))),
              ntree=500)
  # LASSO
  lasso <- train(Y ~ .,
                 data = data[, !(names(data) %in% c("ID", "Time"))],
                 # data = data[, !(names(data) %in% c("Time"))],
                 method = "lasso",
                 metric = "RMSE",
                 trControl = ctrl,
                 tuneGrid = expand.grid(fraction = seq(0.1, 0.9, by = 0.1)),
                 preProc=c("center", "scale"))

  # PLS
  pls <- train(Y ~ .,
               data = data[, !(names(data) %in% c("ID", "Time"))],
               method = "pls",
               metric = "RMSE",
               trControl = ctrl,
               tuneGrid = expand.grid(ncomp = c(1:10)),
               preProc=c("center", "scale"))

  ### STEP 4. Selecting the "best" ML algorithm based on the provided performance metric
  m_list <- list(bo = bo, rf = rf, lasso = lasso, pls = pls)
  rmse_min <- sapply(m_list, FUN = function(x) min(x$results[, "RMSE"]), simplify = T)
  ind <- which(rmse_min == min(rmse_min))

  ### Returning result
  return(best = m_list[[ind]])

}
