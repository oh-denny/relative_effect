sim_data <- function(n = 40, group = rep(c(0, 1), length.out = n)) {
  set.seed(42)
  data.frame(
    y1 = stats::rnorm(n),
    y2 = stats::rnorm(n, mean = group),
    y3 = stats::rnorm(n),
    grp = group
  )
}

test_that("REffect runs and returns the expected list", {
  dat <- sim_data()
  # npmv may warn about a singular rank covariance matrix in some small
  # bootstrap resamples; this does not affect the relative effects.
  res <- suppressWarnings(
    REffect(dat, dep_vars = c("y1", "y2", "y3"), indep_var = "grp",
            n_boot = 20, plot = FALSE, verbose = FALSE)
  )

  expect_type(res, "list")
  expect_named(res, c("obs_manova", "df_perc_boot", "plot"))
  expect_null(res$plot)

  expect_s3_class(res$obs_manova, "data.frame")
  expect_equal(dim(res$obs_manova), c(2L, 3L))
  expect_equal(colnames(res$obs_manova), c("y1", "y2", "y3"))

  expect_length(res$df_perc_boot, 2)
  expect_named(res$df_perc_boot, rownames(res$obs_manova))
  for (d in res$df_perc_boot) {
    expect_s3_class(d, "data.frame")
    expect_true(all(c("Var2", "X2.5.", "X97.5.", "group") %in% names(d)))
    expect_equal(nrow(d), 3L)
    expect_true(all(d$X2.5. <= d$X97.5.))
  }
})

test_that("REffect returns a ggplot object when plot = TRUE", {
  dat <- sim_data()
  pdf(NULL)
  on.exit(grDevices::dev.off())
  res <- suppressWarnings(
    REffect(dat, dep_vars = c("y1", "y2", "y3"), indep_var = "grp",
            n_boot = 20, plot = TRUE, verbose = FALSE)
  )
  expect_s3_class(res$plot, "ggplot")
})

test_that("REffect errors when the independent variable is not 0/1", {
  dat <- sim_data(group = rep(c(1, 2), length.out = 40))
  expect_error(
    REffect(dat, dep_vars = c("y1", "y2", "y3"), indep_var = "grp",
            n_boot = 20, plot = FALSE, verbose = FALSE),
    "only 0 and 1"
  )
})

test_that("REffect errors when a variable is missing from data", {
  dat <- sim_data()
  expect_error(
    REffect(dat, dep_vars = c("y1", "y2", "nope"), indep_var = "grp",
            n_boot = 20, plot = FALSE, verbose = FALSE),
    "not found in `data`: nope"
  )
})
