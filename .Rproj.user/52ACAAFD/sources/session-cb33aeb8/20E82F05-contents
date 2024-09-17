load("./refmod.RData")

test_that("MLCM works", {

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

  set.seed(1)
  causal.1 <- MLCM(data = newdata, int_date = 2019, inf_type = "classic", PCV = pcv.1,
                   nboot = 10, CATE = FALSE, y.lag = 2)

  expect_equal(causal.1$ate, causal$ate)
  expect_equal(causal.1$individual, causal$individual)

})

test_that("MLCMStag works", {

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

  set.seed(1)
  causal_stag.1 <- MLCMStag(data = newdata, int_date = "int_date", inf_type = "block", nboot = 10, PCV = pcv_stag.1, y.lag = 2)

  expect_equal(causal_stag.1$global_ate, causal_stag$global_ate)
  expect_equal(causal_stag.1$tempavg_individual, causal_stag$tempavg_individual)
  expect_equal(causal_stag.1$groupavg, causal_stag$groupavg)
  expect_equal(causal_stag.1$group_ate, causal_stag$group_ate)

})
