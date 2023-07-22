###############################################################################
##
## Panel Cross Validation
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: March 2023
##
###############################################################################

PanelCrossValidation <- function(data, post_period){

  ### Parameter checks
  if(!any(class(data) %in% "PanelMLCM")){stop("Invalid class in the PanelCrossValidation function, something is wrong with as.PanelMLCM")}
  if(!(post_period %in% data[, "Time"]))(stop("post_period must be contained in timevar"))

  ### STEP 1. The CAST package is used to generate separate testing sets for each year

  Tt <- length(unique(data[, "Time"]))
  post_period <- post_period
  indices <- CreateSpacetimeFolds(data, timevar = "Time",
                                  k=Tt)
  trainx <- lapply(1:(Tt-2), FUN = function(x) unlist(indices$indexOut[1:x]))
  testx <- lapply(1:(Tt-2), FUN = function(x) unlist(indices$indexOut[[x+1]]))

  ### STEP 2. Set control function by specifying the training and testing folds that caret will use
  ###         for cross-validation and tuning of the hyperparameters (i.e., the combination of folds defined above)

  ctrl <- trainControl(index = trainx, indexOut = testx)

  ### STEP 3.  Tune the hyperparameters of each of the ML algorithms via temporal cross-validation

  # DOMANDA: di quanti dei successivi parametri vogliamo lasciare la scelta all'utente?
  # Vogliamo lasciare la scelta anche degli algoritmi?? Forse potremmo far sì che l'utente possa inserire le funzioni 'train' che vuole

  # Stochastic gradient boosting
  gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3),
                          n.trees=c(500, 1000, 1500),
                          shrinkage = 0.1,
                          n.minobsinnode = 10)
  bo <- train(Y ~ .,
              data = data[, !(names(data) %in% c("ID", "Time"))],
              # data = data[, !(names(data) %in% c("Time"))],
              method = "gbm",
              metric = "RMSE",
              trControl = ctrl,
              tuneGrid = gbmGrid)
  # Random forest
  rf <- train(Y ~ .,
              data = data[, !(names(data) %in% c("ID", "Time"))],
              # data = data[, !(names(data) %in% c("Time"))],
              method = "rf",
              metric = "RMSE",
              trControl = ctrl,
              ntree=1000)
  # LASSO
  lasso <- train(Y ~ .*.,
                 data = data[, !(names(data) %in% c("ID", "Time"))],
                 # data = data[, !(names(data) %in% c("Time"))],
                 method = "lasso",
                 metric = "RMSE",
                 trControl = ctrl,
                 preProc=c("center", "scale"))

  ### STEP 4. Selecting the "best" ML algorithm based on the provided performance metric
  m_list <- list(bo = bo, rf = rf, lasso = lasso)
  rmse_min <- sapply(m_list, FUN = function(x) min(x$results[, "RMSE"]), simplify = T)
  ind <- which(rmse_min == min(rmse_min))

  ### Returning result
  return(best = m_list[[ind]])

}
