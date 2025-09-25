library(dplyr); library(ggplot2)

plot_tracker <- function(
    dataList,
    tracker,
    food = NULL,
    plot_trials_n = 30,
    y_lims = c(0,1),
    summarize = FALSE,
    alpha_individual = 0.8,
    title = NULL,
    ylab  = NULL
){
  obj <- if (is.character(tracker)) {
    dataList[[tracker]]
  } else tracker
  
  if (is.array(obj) && length(dim(obj)) == 3) {
    stopifnot(!is.null(food))
    mat <- obj[ , , food, drop = FALSE][ , , 1]
  } else mat <- as.matrix(obj)
  
  Tn <- ncol(mat)
  xmax <- min(plot_trials_n, Tn)
  
  df <- as.data.frame(mat)
  colnames(df) <- paste0("trial_", seq_len(ncol(df)))
  df$id <- seq_len(nrow(df))
  df$treatment <- dataList$treatment
  
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2)
  })
  df_long <- df %>%
    tidyr::pivot_longer(cols = dplyr::starts_with("trial_"),
                        names_to = "trial", values_to = "value") %>%
    dplyr::mutate(trial = as.integer(sub("trial_", "", trial)))
  
  if (is.null(ylab)) ylab <- "Value"
  if (is.null(title)) title <- "Tracker across test trials"
  
  p <- ggplot(df_long, aes(x = trial, y = value, group = id, color = treatment)) +
    geom_line(alpha = alpha_individual) +
    coord_cartesian(xlim = c(1, xmax), ylim = y_lims) +
    labs(x = "Test trial", y = ylab, title = title) +
    theme_minimal()
  
  if (summarize) {
    p <- p + stat_summary(aes(group = treatment),
                          fun = mean, geom = "line", linewidth = 1.2, alpha = 0.95)
  }
  p
}
