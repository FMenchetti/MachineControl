###############################################################################
##
## Bootstrap inference
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: March 2023
##
###############################################################################

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
