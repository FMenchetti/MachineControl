###############################################################################
##
## Plotting functions
##
## Authors: Cerqua A., Letta M., Menchetti F.
##
## Date last modified: August 2024
##
###############################################################################

.plot_cohorts <- function(x){

  # Param checks
  if(class(x) != "StagMLCM") stop("this function is a plotting method for 'StagMLCM' objects")

  # Step 1.

}

.plot_groupavg <- function(x){

  # Param checks
  if(class(x) != "StagMLCM") stop("this function is a plotting method for 'StagMLCM' objects")

  # Step 1. Creating a dataframe
  df <- data.frame(times = names(x$groupavg), y = x$groupavg, lower = x$conf.groupavg[1,],
                  upper = x$conf.groupavg[2,])
  df$col <- ifelse(df$times < 2004, "Pre-int","Post-int")

  # Step 2. Ggplotting
  g <- ggplot(df, aes(x = times, y = y, ymin = min(lower), ymax = max(upper))) +
              geom_point(size = 2, aes(colour = col)) +
              geom_errorbar(aes(ymin = lower, ymax = upper, colour = col), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
    labs(title = "Group-average effects over time", colour = " ") + xlab("") + ylab("") +
    scale_color_manual(values = c("Pre-int" = "brown1", "Post-int" = "deepskyblue"))

  # Step 3. Returning results
  return(g)

}
