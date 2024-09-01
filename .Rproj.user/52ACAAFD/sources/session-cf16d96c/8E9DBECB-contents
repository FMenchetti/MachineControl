load("./refmod.RData")

test_that("PanelCrossValidation works", {

  newdata <- as.PanelMLCM(y = data[, "Y"], timevar = data[, "year"], id = data[, "ID"],
                          x = data[, !(names(data) %in% c("Y", "ID", "year"))], y.lag = 1)
  set.seed(1)
  enet <- list(method = "enet",
               tuneGrid = expand.grid(
                 fraction = seq(0.1, 0.9, by = 0.1),
                 lambda = seq(0.1, 0.9, by = 0.1)))
  set.seed(1)
  linreg <- list(method = "lm",
                 tuneGrid = expand.grid(
                   intercept = seq(0, 10, by = 0.5)))
  set.seed(1)
  pls <- list(method = "pls", tuneGrid = expand.grid(ncomp = c(1:5)))

  set.seed(1)
  pcv.1 <- PanelCrossValidation(data = newdata, int_date = 2019,
                                ML_methods = list(enet, linreg, pls))

  expect_equal(pcv.1$best.metric, pcv$best.metric)
  expect_equal(pcv.1$best$results, pcv$best$results)
  expect_equal(pcv.1$best$bestTune, pcv$best$bestTune)
  expect_equal(pcv.1$best$trainingData, pcv$best$trainingData)
  expect_equal(pcv.1$best$resample, pcv$best$resample)

})


test_that("PanelCrossValidationStag works", {

  int_year_i <- c(rep(2019, times = 60), rep(2020, times = 40))
  int_year <- rep(int_year_i, each = length(unique(data$year)))
  data_stag <- data.frame(int_year, data)
  newdata <- as.PanelMLCM(y = data_stag[, "Y"], timevar = data_stag[, "year"], id = data_stag[, "ID"],
                          int_date = data_stag[, "int_year"],
                          x = data_stag[, !(names(data_stag) %in% c("Y", "ID", "year", "int_year"))], y.lag = 2)

  set.seed(1)
  enet <- list(method = "enet",
               tuneGrid = expand.grid(
                 fraction = seq(0.1, 0.9, by = 0.1),
                 lambda = seq(0.1, 0.9, by = 0.1)))
  set.seed(1)
  linreg <- list(method = "lm",
                 tuneGrid = expand.grid(
                   intercept = seq(0, 10, by = 0.5)))
  set.seed(1)
  pls <- list(method = "pls", tuneGrid = expand.grid(ncomp = c(1:5)))

  set.seed(1)
  pcv_stag.1 <- PanelCrossValidationStag(data = newdata, ML_methods = list(enet, linreg, pls))

  expect_equal(lapply(pcv_stag.1, FUN = function(x)(x$best.metric)), lapply(pcv_stag, FUN = function(x)(x$best.metric)))
  expect_equal(lapply(pcv_stag.1, FUN = function(x)(x$best$results)), lapply(pcv_stag, FUN = function(x)(x$best$results)))
  expect_equal(lapply(pcv_stag.1, FUN = function(x)(x$best$bestTune)), lapply(pcv_stag, FUN = function(x)(x$best$bestTune)))
  expect_equal(lapply(pcv_stag.1, FUN = function(x)(x$best$trainingData)), lapply(pcv_stag, FUN = function(x)(x$best$trainingData)))
  expect_equal(lapply(pcv_stag.1, FUN = function(x)(x$best$resample)), lapply(pcv_stag, FUN = function(x)(x$best$resample)))

})
