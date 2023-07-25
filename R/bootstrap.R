###############################################################################
##
## Bootstrap inference
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: July 2023
##
###############################################################################

boot_ate <- function(data, ind, bestt, type, nboot, ate = NULL){

  ### Step 1. Sampling indices
  if(type %in% c("classic", "bc classic", "bca")){

    ii<- matrix(sample(ind, size = nboot*length(ind), replace = T), nrow = nboot, ncol = length(ind))
    # Nota: non sto ricampionando un sottoinsieme (es. 60% del totale) ma tutto
  }

  if(type %in% c("block", "bc block")){

    ids <- unique(data$ID)
    ii <- t(sapply(1:nboot, function(boot){set.seed(boot); ind1 <- sample(ids, size = length(ids), replace = T);
                                           unlist(sapply(ind1, function(x)(intersect(which(data[, "ID"] %in% x), ind)), simplify = F))}))
    # ind1 <- sample(unique(data$ID), round(length(unique(data$ID))*0.6), replace = T)
    # ind1 <- intersect(which(data[, "ID"] %in% ind1), ind)
    # Nota: anche qui non sto ricampionando un sottoinsieme ma tutte le unità con replacement
  }

  ### Step 2. Bootstrapping (estimating the causal effect on the units resampled in each bootstrap iteration)
  ate_boot <- apply(ii, 1, function(i){
    fitb <- train(Y ~ .,
                  data = data[i, !(names(data) %in% c("ID", "Time"))],
                  method = bestt$method,
                  metric = "RMSE",
                  trControl = trainControl(method="none"),
                  tuneGrid = bestt$bestTune) ;
    obs <- data[-c(ind, i), "Y"] ;
    eps <- data[i, "Y"] - predict(fitb)
    # error <- sample(eps, size = nrow(data[-ind,]), replace = T) # verifica che sia così anche per block boot
    error <- rnorm(n = nrow(data[-ind,]), mean = mean(eps), sd = sd(eps))
    pred <- predict(fitb, newdata = data[-c(ind, i), ]) + error ;
    #pred <- predict(fitb, newdata = data[-c(ind, i), ]) ;
    obs - pred
  })

  mean_ate_boot <- colMeans(ate_boot)
  conf.ate <- quantile(mean_ate_boot, probs = c(0.025, 0.975))

  ### Step 3. Adjusting for bias and/or skewness (if 'type' is "bc classic", "bc block" or "bca")

  if(type %in% c("bc classic", "bc block")){

    # Bias correction
    z0 <- qnorm(sum(mean_ate_boot < ate)/nboot)
    lower <- pnorm(2*z0 + qnorm(0.025))
    upper <- pnorm(2*z0 + qnorm(0.975))
    conf.ate <- quantile(mean_ate_boot, probs = c(lower, upper))

  }

  if(type == "bca"){

    counts <- t(apply(ii, 1, FUN = function(x)(table(c(x, ind))-1)))
    Blist <- list(Y = counts, tt = colMeans(ate_boot), t0 = ate)
    out2 <- bcajack2(B = Blist) ## doesn't work
    conf.ate <- out2$lims[c(1,9),"bca"]
    # or
    # func <- function(x, bestt, ind){

    # fitb <- train(Y ~ .,
    #              data = x[ind, !(names(data) %in% c("ID", "Time"))],
    #             method = bestt$method,
    #              metric = "RMSE",
    #             trControl = trainControl(method="none"),
    #             tuneGrid = bestt$bestTune) ;
    # obs <- x[-ind, "Y"] ;
    # pred <- predict(fitb, newdata = x[-ind, ]) ;
    # mean(obs - pred)
    # }
    # out <- bcapar(t0 = ate, tt = rowMeans(ate_boot), bb = counts)
  }


  #fitb <- train(Y ~ .,
  #              data = data[ind1, !(names(data) %in% c("ID", "Time"))],
  #              method = bestt$method,
  #              metric = "RMSE",
  #              trControl = trainControl(method="none"),
  #              tuneGrid = bestt$bestTune)
  # obs <- data[-c(ind, ind1), "Y"]
  # pred <- predict(fitb, newdata = data[-c(ind, ind1), ])

  # Returning results
  return(list(type = type, ate_boot = ate_boot, conf.ate = conf.ate, var.ate = var(mean_ate_boot), ate.lower = conf.ate[1], ate.upper = conf.ate[2]))
}

boot_cate <- function(effect, cate, nboot){

  terminal.nodes <- cate$where
  x <- unique(terminal.nodes)
  node.inf <- mapply(x, FUN = function(x){y <- effect[which(terminal.nodes == x)];
                                          boot.dist <- matrix(sample(y, size = nboot*length(y), replace = TRUE),
                                                              nrow = nboot, ncol = length(y));
                                          mean.cate <- rowMeans(boot.dist);
                                          var.cate <- var(mean.cate);
                                          conf.cate <- quantile(mean.cate, probs = c(0.025, 0.975));
                                          c(cate = mean(y), var.cate = var.cate, cate.lower = conf.cate[1], cate.upper = conf.cate[2])},
              SIMPLIFY = TRUE)
  colnames(node.inf) <- paste0("Node_", x)
  return(node.inf)

}
