# Tests for evaluation_measures (MCI and Wasserstein)

test_that("MCI returns expected values in simple cases", {
  set.seed(1)
  rc <- data.table(x = 1:5, y = 1:5)
  mp <- copy(rc) # perfectly concordant
  mci_val <- calc_mci(rc, mp)
  
  expect_equal(mci_val$mci, rep(0, 5)) # perfect concordance
})

test_that("calc_wasserstein behaves sensibly", {
  rc <- data.frame(x = 1:5)
  mp <- data.frame(x = 1:5)
  dist_same <- calc_wasserstein(rc, mp)
  expect_equal(unname(dist_same["Wasserstein_1"]), 0)
  expect_equal(unname(dist_same["Wasserstein_2"]), 0)
  
  mp_shift <- data.frame(x = 2:6)
  dist_shift <- calc_wasserstein(rc, mp_shift)
  expect_equal(unname(dist_shift["Wasserstein_1"]), 1, tolerance = 1e-8)
  expect_equal(unname(dist_shift["Wasserstein_2"]), 1, tolerance = 1e-8)
  
  rc <- data.frame(x = 1:5, y = 1:5)
  mp <- data.frame(x = 1:5, y = 1:5)
  dist_same <- calc_wasserstein(rc, mp)
  expect_equal(unname(dist_same["Wasserstein_1"]), 0)
  expect_equal(unname(dist_same["Wasserstein_2"]), 0)
  dist_same <- calc_wasserstein(rc, mp, scale_dta = "both")
  expect_equal(unname(dist_same["Wasserstein_1"]), 0)
  expect_equal(unname(dist_same["Wasserstein_2"]), 0)
  dist_same <- calc_wasserstein(rc, mp, scale_dta = "sd")
  expect_equal(unname(dist_same["Wasserstein_1"]), 0)
  expect_equal(unname(dist_same["Wasserstein_2"]), 0)

})
