# REffect: An R package for bootstrap-based estimation and visualization of nonparametric relative effects in biomedical research

**REffect** is an R package for estimating **nonparametric relative effects** between
two groups across several dependent variables, with **percentile bootstrap confidence
intervals** and a profile plot of the effects.

It works with any set of numeric (or ordinal) outcomes — biomarkers, questionnaire
items, gene expression, test scores, etc. — measured in two groups.

## What it does

`REffect()`:

1. Fits a nonparametric multivariate test with `npmv::nonpartest()` comparing two
   groups (independent variable coded `0`/`1`) on several dependent variables, and
   extracts the two-group relative effects (`twogroupreleffects`).
2. Resamples the rows of the data with replacement `n_boot` times and recomputes the
   relative effects in each resample.
3. Computes the 2.5% and 97.5% percentiles (95% CI) of each relative effect, per group.
4. Optionally draws (with `ggplot2`) the profile of the observed relative effects with
   the bootstrap confidence bands of each group.

A relative effect of 0.5 means no tendency for that group to have larger or smaller
values than the pooled sample. Values above (below) 0.5 indicate that the group tends
to have larger (smaller) values on that variable.

## Installation

```r
# install.packages("devtools")
devtools::install_github("oh-denny/relative_effect")
```

## Usage

```r
library(REffect)

set.seed(123)
dat <- data.frame(
  y1    = rnorm(60),
  y2    = rnorm(60, mean = rep(c(0, 1), each = 30)),
  y3    = rnorm(60),
  group = rep(c(0, 1), each = 30)
)

res <- REffect(
  data      = dat,
  dep_vars  = c("y1", "y2", "y3"),
  indep_var = "group",
  n_boot    = 500,
  verbose   = FALSE
)

res$obs_manova    # observed relative effects (rows = groups)
res$df_perc_boot  # 95% bootstrap CI per group
res$plot          # ggplot object
```

### Customizing the plot

```r
res <- REffect(
  dat, c("y1", "y2", "y3"), "group",
  n_boot               = 500,
  plot_colors          = c("grey40", "firebrick"),
  plot_labels          = list(y = "Relative effect", x = "Outcome"),
  plot_theme           = ggplot2::theme_minimal(),
  plot_legend_position = "bottom",
  plot_legend_labels   = c("Control", "Treatment")
)

# Save the plot
ggplot2::ggsave("releffects.png", res$plot, width = 8, height = 5)
```

Since the plot is returned as a `ggplot` object, you can keep adding layers to it
(e.g. `res$plot + ggplot2::ggtitle("My study")`).

## Return value

An (invisible) list with:

| Element        | Content                                                                |
|----------------|------------------------------------------------------------------------|
| `obs_manova`   | `data.frame` with the observed relative effects (one group per row).   |
| `df_perc_boot` | list with one `data.frame` per group: `Var2`, `X2.5.`, `X97.5.`, `group`. |
| `plot`         | `ggplot` object (or `NULL` if `plot = FALSE`).                         |

## Notes

- The independent variable must contain only `0` and `1`; otherwise the function
  stops with an error.
- With small samples, `npmv::nonpartest()` may emit the warning
  *"Rank covariance matrix is singular, only ANOVA test returned"* in some bootstrap
  resamples. This does not affect the relative effects.
- Use `set.seed()` before calling `REffect()` for reproducible confidence intervals.

## License

MIT © Denny Leandro
