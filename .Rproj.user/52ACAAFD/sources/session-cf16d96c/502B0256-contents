test_that(".sim_ardl stops when wrong param are set", {

  ## param
  post_per <- 1                   # Number of post-intervention periods (excluding int_date)
  N <- 10                         # Number of units in the panel
  beta <- c(1, 2.5, 0.1)          # Coefficient of Xt-1
  rho <- 0.8                      # Coefficient of Yt-1
  sigma <- 2                      # St.dev of the error term
  impact <- c(2,1.5)              # Fictional additive effect (it can be fixed or vary)
  t <- 7
  X <- lapply(1:N, FUN = function(x)(cbind(x1 = seq(0.2, by = 0.1, length.out = t) + rnorm(t, 0, 1),
                                           x2 = sample(0:1, t, replace = TRUE),
                                           x3=sample(1:3, t, replace = TRUE))))
  ## tests
  expect_error(.sim_ardl(seed = 1, beta = beta, X = do.call(cbind, X), N = N, sigma = sigma,
                         ar1_coef = rho, impact = impact, impact_constant = impact_constant,
                         linear = linear, post_per = post_per))
  expect_error(.sim_ardl(seed = 1, beta = as.factor(beta), X = X, N = N, sigma = sigma,
                         ar1_coef = rho, impact = impact, impact_constant = impact_constant,
                         linear = linear, post_per = post_per))
  expect_error(.sim_ardl(seed = 1, beta = beta, X = X, N = N, sigma = sigma,
                         ar1_coef = rho, impact = impact, impact_constant = impact_constant,
                         linear = linear, post_per = 2))
})
