## Code to prepare `DATA` dataset

#### STEP 1. Setting arguments for function .sim_ardl()

## Parameters
seed <- 1                       # Seed to ensure exact replication of the dataset
post_per <- 1                   # Number of post-intervention periods (excluding int_date)
N <- 100                        # Number of units in the panel
beta <- c(0, 1, 2.5, 0.1,
          0, 0, 2, 1.5)         # Coefficient of Xt-1
rho <- 0.8                      # Coefficient of Yt-1
sigma <- 2                      # St.dev of the error term

## Linear or non linear model?
linear <- TRUE

## Impact settings
impact <- c(2,1.5)              # Fictional additive effect (it can be fixed or vary)
impact_constant <- FALSE        # Set to FALSE for additive impact proportional to std.dev(Yt).


#### STEP 2. Simulating covariates

t <- 7

## Generating continuous covariates, varying across units
set.seed(seed)
x1 <- seq(0.2, by = 0.1, length.out = t) + rnorm(t, 0, 1)
xm <- MASS::mvrnorm(n = t, mu = c(1,2,3), Sigma = matrix(c(1,0.5,0.7,0.5, 1, 0.3, 0.7, 0.3, 1), nrow = 3, ncol = 3))
X <- cbind(x1 = x1, x2 = xm[,1], x3 = xm[,2], x4 = xm[,3])
ran <- lapply(1:N, FUN = MASS::mvrnorm, n = NROW(X) , mu = rep(1, ncol(X)), Sigma = 0.1*diag(1, nrow = ncol(X), ncol = ncol(X)))
X <- lapply(ran, function(x){y <- x + X; colnames(y) <- colnames(X); y})

## Generating categorical covariates, varying across units
x5 <- lapply(1:N, FUN = function(x)(sample(0:1, t, replace = TRUE)))
x6 <- lapply(1:N, FUN = function(x)(sample(1:3, t, replace = TRUE)))

## Binding and adding interactions
X <- mapply(X, x5, x6, FUN = function(X, x5, x6)(cbind(X, x5, x6, x7 = X[,3]*x6, x8 = X[,2]*x5)), SIMPLIFY = F)

#### STEP 3. Single dataset generation

data <- .sim_ardl(seed = seed, beta = beta, X = X, N = N, sigma = sigma,
                  ar1_coef = rho, impact = impact, impact_constant = impact_constant,
                  linear = linear, post_per = post_per)
data <- data$dat

usethis::use_data(data, overwrite = TRUE)
