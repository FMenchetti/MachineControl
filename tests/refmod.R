### Code to prepare `refexamples.RData`

## Non-staggered setting: Example 2 of 'PanelCrossValidation()'

# Organizing the dataset
newdata <- as.PanelMLCM(y = data[, "Y"], timevar = data[, "year"], id = data[, "ID"],
                        x = data[, !(names(data) %in% c("Y", "ID", "year"))], y.lag = 1)

# PCV
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
pcv <- PanelCrossValidation(data = newdata, int_date = 2019,
                            ML_methods = list(enet, linreg, pls))

# Causal effect estimation and inference
set.seed(1)
causal <- MLCM(data = newdata, int_date = 2019, inf_type = "classic", PCV = pcv,
               nboot = 10, CATE = FALSE, y.lag = 2)

## Staggered setting : Example 1 of 'PanelCrossValidationStag()'

# Assume the following intervention dates
int_year_i <- c(rep(2019, times = 60), rep(2020, times = 40))
int_year <- rep(int_year_i, each = length(unique(data$year)))

# Define data_stag
data_stag <- data.frame(int_year, data)

# Organizing the dataset with as.PanelMLCM
newdata <- as.PanelMLCM(y = data_stag[, "Y"], timevar = data_stag[, "year"], id = data_stag[, "ID"],
                        int_date = data_stag[, "int_year"],
                        x = data_stag[, !(names(data_stag) %in% c("Y", "ID", "year", "int_year"))], y.lag = 2)

# Panel Cross Validation in a staggered setting with different ML methods
set.seed(1)
pcv_stag <- PanelCrossValidationStag(data = newdata, ML_methods = list(enet, linreg, pls))

# Causal effect estimation and inference
set.seed(1)
causal_stag <- MLCMStag(data = newdata, int_date = "int_date", inf_type = "block", nboot = 10, PCV = pcv_stag, y.lag = 2)

## Saving
save(pcv, pcv_stag, causal, causal_stag, file = "./tests/testthat/refmod.RData")
