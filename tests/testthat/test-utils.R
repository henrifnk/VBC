test_that("calc_wasserstein returns zero distance for identical distributions", {
  # Two identical datasets (observed and model)
  obs <- data.frame(x = 1:5)
  mod <- data.frame(x = 1:5)  # same values as obs
  dist <- calc_wasserstein(obs, mod)
  # Should return a numeric vector of length 2 (Wasserstein_1 and Wasserstein_2)
  expect_type(dist, "double")
  expect_equal(length(dist), 2)
  # Both distances should be 0 for identical distributions
  expect_equal(unname(dist["Wasserstein_1"]), 0)
  expect_equal(unname(dist["Wasserstein_2"]), 0)
})

test_that("calc_wasserstein captures distribution differences (simple shifted distribution)", {
  # Observed vs model data where model is a constant shift of +1
  obs <- data.frame(x = 1:5)
  mod <- data.frame(x = 2:6)  # model values are observed + 1
  dist <- calc_wasserstein(obs, mod)
  # Expect Wasserstein-1 and Wasserstein-2 distances to be approximately 1 (the shift)
  expect_equal(unname(dist["Wasserstein_1"]), 1, tolerance = 1e-8)
  expect_equal(unname(dist["Wasserstein_2"]), 1, tolerance = 1e-8)
})

test_that("plot_tails returns a ggplot (or patchwork) object for given data", {
  # Use the sample climate data (if available) to test plot_tails
  data("climate", package = "VBC")
  # Plot tails for precipitation ("pr") in model projection data
  plt <- plot_tails(climate$mp, var = "pr", scale_d = 0.5, mult = 2, xmin = 0)
  # The returned object should inherit from ggplot
  expect_true(inherits(plt, "ggplot") || inherits(plt, "gg"))
})
