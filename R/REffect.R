#' Bootstrap confidence bands for two-group nonparametric relative effects
#'
#' Fits a nonparametric multivariate test with [npmv::nonpartest()] comparing
#' two groups on several dependent variables, extracts the two-group relative
#' effects, and estimates 95% percentile bootstrap confidence intervals for
#' each relative effect by resampling the rows of `data` with replacement.
#' Optionally draws a profile plot of the relative effects with the bootstrap
#' bands for each group.
#'
#' A relative effect of 0.5 means no tendency for that group to have larger or
#' smaller values than the pooled sample; values above (below) 0.5 indicate a
#' tendency towards larger (smaller) values.
#'
#' @param data A data frame containing the dependent variables and the
#'   grouping variable.
#' @param dep_vars Character vector with the names of the (numeric) dependent
#'   variables.
#' @param indep_var Character string with the name of the independent
#'   (grouping) variable. It must contain only the values `0` and `1`.
#' @param n_boot Number of bootstrap resamples. Default `1000`.
#' @param plot Logical. If `TRUE` (default), the profile plot is built and
#'   printed.
#' @param plot_colors Character vector of length 2 with the colors for group
#'   `0` and group `1`, respectively.
#' @param plot_labels Named list with elements `y` and `x` giving the axis
#'   titles.
#' @param plot_theme A ggplot2 theme object applied to the plot. Default
#'   [ggplot2::theme_classic()].
#' @param plot_legend_position Legend position passed to
#'   [ggplot2::theme()] (e.g. `"top"`, `"bottom"`, `"none"`).
#' @param plot_legend_labels Character vector of length 2 with the legend
#'   labels for group `0` and group `1`, respectively.
#' @param verbose Logical. If `TRUE` (default), prints the progress of the
#'   bootstrap loop.
#'
#' @return Invisibly, a list with three elements:
#' \describe{
#'   \item{obs_manova}{Data frame with the observed two-group relative effects
#'     (one row per group, one column per dependent variable).}
#'   \item{df_perc_boot}{List with one data frame per group containing, for
#'     each dependent variable (`Var2`), the 2.5% (`X2.5.`) and 97.5%
#'     (`X97.5.`) bootstrap percentiles.}
#'   \item{plot}{The `ggplot` object, or `NULL` when `plot = FALSE`.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' dat <- data.frame(
#'   y1    = rnorm(60),
#'   y2    = rnorm(60, mean = rep(c(0, 1), each = 30)),
#'   y3    = rnorm(60),
#'   group = rep(c(0, 1), each = 30)
#' )
#' res <- REffect(dat, dep_vars = c("y1", "y2", "y3"),
#'                indep_var = "group", n_boot = 200)
#' res$obs_manova
#' res$df_perc_boot
#' res$plot
#' }
#'
#' @export
REffect <- function(data, dep_vars, indep_var, n_boot = 1000,
                    plot = TRUE, plot_colors = c("grey20", "darkblue"),
                    plot_labels = list(y = "Relative effect", x = "Variable"),
                    plot_theme = ggplot2::theme_classic(), plot_legend_position = "top",
                    plot_legend_labels = c("Group 0", "Group 1"), verbose = TRUE) {

  missing_vars <- setdiff(c(dep_vars, indep_var), names(data))
  if (length(missing_vars) > 0) {
    stop("Variable(s) not found in `data`: ", paste(missing_vars, collapse = ", "))
  }

  if (!all(data[[indep_var]] %in% c(0, 1))) {
    stop("The independent variable must contain only 0 and 1 (two groups).")
  }

  n_obs <- nrow(data)
  list_boot <- list()

  formula_str <- paste(dep_vars, collapse = " | ")
  formula <- stats::as.formula(paste(formula_str, "~", indep_var))

  # Bootstrap loop (i = 0 is the observed, non-resampled data)
  for (i in 0:n_boot) {
    if (i == 0) {
      df_sample <- 1:n_obs
    } else {
      df_sample <- sample(1:n_obs, n_obs, replace = TRUE)
    }

    manova_fit <- npmv::nonpartest(formula, data = data[df_sample, ], permtest = FALSE,
                                   plots = FALSE, permreps = 1000, tests = c(1, 0, 0, 0))

    list_boot[[i + 1]] <- manova_fit

    if (verbose) cat("Bootstrap iteration", i, "of", n_boot, "\n")
  }

  obs_manova <- list_boot[[1]]$twogroupreleffects
  row_manova <- rownames(obs_manova)

  releffects <- lapply(list_boot[-1], "[[", 2)

  # 2.5% and 97.5% bootstrap percentiles for each group (i) and variable (j)
  df_perc_boot <- lapply(seq_len(nrow(obs_manova)), function(i) {
    q <- sapply(seq_len(ncol(obs_manova)), function(j) {
      v <- sapply(releffects, function(r) {
        row <- rownames(r)
        if (row[1] != row_manova[1]) r <- r[2:1, ]
        r[i, j]
      })
      stats::quantile(v, probs = c(0.025, 0.975))
    })
    colnames(q) <- colnames(obs_manova)
    df <- data.frame(Var2 = colnames(q), t(q),
                     group = ifelse(i == 1, row_manova[1], row_manova[2]))
    df$Var2 <- factor(df$Var2, levels = sort(unique(df$Var2)), ordered = TRUE)
    return(df)
  })

  names(df_perc_boot) <- row_manova

  # Long-format point estimates, ordered by group and decreasing relative effect
  df <- reshape2::melt(
    cbind(obs_manova, group = row_manova),
    id.vars = "group"
  )
  df$Var1 <- "point_est"
  df <- df[, c(4, 2, 3, 1)]
  colnames(df)[2] <- "Var2"

  df <- dplyr::arrange(df, group, dplyr::desc(value), Var2)
  df$Var2 <- as.character(df$Var2)
  df$Var2 <- factor(df$Var2, levels = unique(df$Var2), ordered = TRUE)

  for (k in seq_len(2)) {
    df_perc_boot[[k]]$Var2 <- factor(df_perc_boot[[k]]$Var2, levels = levels(df$Var2),
                                     ordered = TRUE)
  }

  p <- NULL
  if (plot) {
    p <- ggplot2::ggplot(data = df, ggplot2::aes(x = Var2, y = value, group = group,
                                                 color = group)) +
      ggplot2::geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[1]],
                           ggplot2::aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group),
                           fill = plot_colors[1], alpha = 0.3) +
      ggplot2::geom_ribbon(inherit.aes = FALSE, data = df_perc_boot[[2]],
                           ggplot2::aes(x = Var2, ymin = X2.5., ymax = X97.5., group = group),
                           fill = plot_colors[2], alpha = 0.3) +
      ggplot2::geom_line() +
      ggplot2::geom_point(ggplot2::aes(size = value)) +
      ggplot2::scale_colour_manual("", values = c("1" = plot_colors[2], "0" = plot_colors[1]),
                                   labels = plot_legend_labels) +
      plot_theme +
      ggplot2::theme(
        legend.text = ggplot2::element_text(size = 12),
        axis.text.x = ggplot2::element_text(size = 12, angle = 90),
        axis.text.y = ggplot2::element_text(size = 12),
        axis.line.y = ggplot2::element_blank(),
        axis.title = ggplot2::element_text(size = 15),
        panel.grid.major = ggplot2::element_line(),
        legend.position = plot_legend_position,
        legend.box.background = ggplot2::element_rect()
      ) +
      ggplot2::labs(y = plot_labels$y, x = plot_labels$x)

    print(p)
  }

  invisible(list(obs_manova = obs_manova, df_perc_boot = df_perc_boot, plot = p))
}
