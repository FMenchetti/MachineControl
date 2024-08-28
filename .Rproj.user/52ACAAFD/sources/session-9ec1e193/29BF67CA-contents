###############################################################################
##
## Plotting functions
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: August 2024
##
###############################################################################

#' Plot of estimated causal effects
#'
#'
#' @param x     Object returned from a previous call to \code{MLCM()}
#' @param type  Character string indicating the plot to be produced. Possible values
#'              in \code{c("ate", "cate")}
#' @param ...   Optional arguments from the \code{rpart.plot()} function to be used while
#'              plotting the regression tree. See \code{help(rpart.plot)} for details.
#' @details
#' The option \code{"ate"} plots the ATE at each point in time after the intervention.
#' Note that the effects plotted before the introduction of the policy are placebo effects (they
#' should be non-significant effects centered around zero otherwise either the intervention has
#' some anticipation effects or the model is misspecified). The option \code{"cate"} draws
#' the estimated regression tree using the \code{rpart.plot()} function from the \code{rpart.plot}
#' package. Users can use the \code{...} argument for full flexibility in drawing CATE regression
#' by exploiting all arguments of \code{rpart.plot()}.
#'
#'
#' @return NULL, invisibly
#' @export
#' @importFrom rpart.plot rpart.plot
#'
plot.MLCM <- function(x, type = c("ate", "cate"), ...){

  # Param checks
  if(class(x) != "MLCM") stop("this function is a plotting method for 'MLCM' objects")
  if(!all(type %in% c("ate", "cate")))
    stop("allowed 'type' values are 'ate' or 'cate'")

  # ATE
  if("ate" %in% type){

    ate <- .plot_effects(x$ate, int_date = x$int_date) +
      ggtitle("ATE")

  } else {ate <- NULL}


  # CATE

  if("cate" %in% type){

    cate <- lapply(1:length(x$cate), function(i)(rpart.plot(x$cate[[i]], main = paste("CATE,", names(x$cate)[i]), ...)))

  } else {cate <- NULL}

  # Returning results
  res <- list(ate = ate, cate = cate)
  return(res)
}

#' Plot of estimated causal effects in staggered adoption setting
#'
#'
#' @param x     Object returned from a previous call to \code{StagMLCM()}
#' @param type  Character string indicating the plot to be produced. Possible values
#'              in \code{c("global", "cohort", "time")}
#' @details
#' The option \code{"global"} plots the global ATE, i.e., a weighted average (by cohort size)
#' across all treatment groups of the temporal average causal effect of the policy. The option
#' \code{"cohort} plots the ATE over the entire post-intervention period separately for the
#' different cohorts. The option \code{"time} draws a weighted average (by group size) of the
#' across all treatment groups of the estimated ATEs in each time point after the intervention.
#' Note that the effects plotted before the introduction of the policy are placebo effects (they
#' should be non-significant effects centered around zero otherwise either the intervention has
#' some anticipation effects or the model is misspecified).
#'
#'
#' @return NULL, invisibly
#' @export
#'
plot.StagMLCM <- function(x, type = c("global", "cohort", "time")){

  # Param checks
  if(class(x) != "StagMLCM") stop("this function is a plotting method for 'StagMLCM' objects")
  if(!all(type %in% c("global", "cohort", "time")))
    stop("allowed 'type' values are 'global', 'cohort' and 'time'")

  # Global ATE
  if("global" %in% type){

    global <- .plot_effects(x$global_ate, int_date = min(x$int_date)) +
              ggtitle("Global ATE")

  } else {global <- NULL}


  # ATE for each cohort

  if("cohort" %in% type){

    cohort <- mapply(X = x$group_ate, Y = as.list(x$int_date),
                     FUN = function(X,Y)(.plot_effects(x = X, int_date = Y) +
                                         ggtitle(paste("ATE for treatment cohort", Y))),
                     SIMPLIFY = FALSE)

  } else {cohort <- NULL}

  # ATE at each time post-intervention
  if("time" %in% type){

    time <- .plot_effects(x$groupavg, int_date = min(x$int_date)) +
            ggtitle("ATE, weighted average across treatment cohorts")

  } else {time <- NULL}

  # Returning results
  res <- list(global = global, cohort = cohort, time = time)
  return(res)
}

#' Internal function for plotting
#'
#' @param x         Object returned from a previous call to \code{MLCM()} or \code{StagMLCM()}
#' @param int_date  Intervention date (or dates in case of staggered adoption of the policy)
#'
#' @return NULL, invisibly
#' @noRd
#'
.plot_effects <- function(x, int_date){

  # Param checks
  # if(class(x) != "StagMLCM") stop("this function is a plotting method for 'StagMLCM' objects")

  # Step 1. Creating a dataframe
  df <- data.frame(times = names(x$estimate), y = x$estimate, lower = x$conf.interval[1,],
                  upper = x$conf.interval[2,])
  df$col <- ifelse(df$times < int_date, "Placebo","Post-int")

  # Step 2. Gg-plotting
  g <- ggplot(df, aes(x = times, y = y, ymin = min(lower), ymax = max(upper))) +
              geom_point(size = 2, aes(colour = col)) +
              geom_errorbar(aes(ymin = lower, ymax = upper, colour = col), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
    labs(colour = " ") + xlab("") + ylab("") +
    scale_color_manual(values = c("Placebo" = "brown1", "Post-int" = "deepskyblue"))

  # Step 3. Returning results
  return(g)

}
