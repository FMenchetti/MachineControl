###############################################################################
##
## MLCM Package
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: March 2023
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
#' @param nboot       Number of bootstrap replications, defaults to 1000.
#'
#' @return A list
#' @export
#'
#'
#'
MLCM <- function(data, y, timevar, id, post_period, nboot = 1000){

  ### Parameter checks
  if(!any(class(data) %in% c("matrix", "data.frame"))) stop("data must be a matrix or a data.frame")
  if(class(y) != "character") stop("y must be a character")
  if(class(timevar) != "character") stop("timevar must be a character")
  if(class(id) != "character") stop("id must be a character")
  if(!any(class(post_period) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt", "numeric", "integer"))) stop("post_period must be integer, numeric or Date")

  ### Structuring the panel dataset in the required format
  data_panel <- as.PanelMLCM(y = data[, y], timevar = data[, timevar], id = data[, id],
                             x = data[, !(names(data) %in% c(y, id, timevar))])

  ### Panel cross-validation
  best <- PanelCrossValidation(data = data_panel, post_period = post_period)

  ### Fit the best (optimized) ML algorithm on all pre-intervention data and make predictions in the post-intervention period
  ind <- which(data_panel[, "Time"] < post_period)
  fit <- train(Y ~ .,
               data = data_panel[, !(names(data_panel) %in% c("ID", "Time"))],
               # data = data_panel[ind, !(names(data) %in% c("Time"))],
               method = best$method,
               metric = "RMSE",
               trControl = trainControl(method="none"),
               tuneGrid = best$bestTune)
  obs <- data_panel[-ind, "Y"]
  pred <- predict(fit, newdata = data_panel[-ind, ])

  ### ATE
  ate <- mean(obs - pred)

  ### Bootstrap

  # IL SEGUENTE CODICE PER IL BOOTSTRAP NON FUNZIONA, LA STATISTICA STIMATA (prova$t0) RISULTA NaN
  # MENTRE INVECE prova$t VIENE STIMATO MA PROSSIMO A ZERO PER OGNI REPLICATION
  # boot_fun <- function(data, ind=ind, bestt=fit){
  #
  #  fitb <- train(Y ~ .,
  #               # data = data[, !(names(data) %in% c("ID", "Time"))],
  #               data = data[ind, !(names(data) %in% c("Time"))],
  #               method = bestt$method,
  #               metric = "RMSE",
  #               trControl = trainControl(method="none"),
  #               tuneGrid = bestt$bestTune)
  #  obs <- data[-ind, "Y"]
  #  pred <- predict(fitb, newdata = data[-ind, ])
  #  ate <- mean(obs - pred)

  # }
  # prova <- boot(data = data_panel, boot_fun, R = 100)

  # SE PERO' HO BEN CAPITO COME FUNZIONA IL CODICE SOPRA, NON SI FA ALTRO CHE CAMPIONARE DEGLI INDICI
  # NEL PERIODO PRE-INTERVENTO E RISTIMARE CON IL METODO SELEZIONATO IN PRECEDENZA ( 'fit' nel nostro caso)
  # SE E' COSì, IL CODICE SEGUENTE DOVREBBE BASTARE:

  ate_boot <- unlist(lapply(1:nboot, function(x){set.seed(x); boot_fun( data = data_panel, ind = ind, bestt = fit)}))


  ### Saving results
  return(list(best_method = best, fit = fit, ate = ate, ate_boot = ate_boot))

}

as.PanelMLCM <- function(y, x, timevar, id){

  # Parameter checks
  if(!(is.numeric(y) & length(y)>1)) stop("y must be a numeric vector of length greater than 1")
  if(!any(class(x) %in% c("numeric", "matrix", "data.frame"))) stop("x must be a vector, matrix or data.frame")
  if(NROW(x) != length(y)) stop("NROW(x) != length(y)")
  if(!any(class(timevar) %in% c("Date", "POSIXct", "POSIXlt", "POSIXt", "integer", "numeric")) | length(timevar) != length(y)) stop("timevar must be a numeric vector or a 'Date' object of the same length as y")
  if(!(is.numeric(id) & length(id) == length(y))) stop("id must be a numeric vector of the same length as y")
  if(length(unique(id))*length(unique(timevar)) != length(y)) stop("The panel is unbalanced")

  # Structuring the panel dataset
  panel <- data.frame(Time = timevar, ID = id, Y = y, x)

  # Returning results
  class(panel) <- c("data.frame", "PanelMLCM")
  return(panel)
}

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

boot_fun <- function(data, ind, bestt){

  ind1<- sample(ind, length(ind)*0.9, replace = F)
  fitb <- train(Y ~ .,
                data = data[ind1, !(names(data) %in% c("ID", "Time"))],
                # data = data[ind1, !(names(data) %in% c("Time"))],
                method = bestt$method,
                metric = "RMSE",
                trControl = trainControl(method="none"),
                tuneGrid = bestt$bestTune)
  obs <- data[-c(ind, ind1), "Y"]
  pred <- predict(fitb, newdata = data[-c(ind, ind1), ])
  ate <- mean(obs - pred)

}
