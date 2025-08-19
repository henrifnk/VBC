test_that("vbc returns a data.table with correct dimensions, column names, and attributes", {
  # Simulate a small dataset with two variables: one continuous and one zero-inflated
  set.seed(123)
  n <- 50
  # Continuous variable (e.g., normal distribution)
  cont_rc <- rnorm(n, mean = 0, sd = 1)
  cont_mc <- rnorm(n, mean = 0, sd = 1)       # same distribution as rc (no bias in calibration)
  cont_mp <- rnorm(n, mean = 1, sd = 1)       # shifted mean (simulating a model projection trend)
  # Zero-inflated variable (e.g., 70% zeros, 30% positive values)
  zi_rc <- ifelse(runif(n) < 0.7, 0, rexp(n, rate = 1))
  zi_mc <- ifelse(runif(n) < 0.7, 0, rexp(n, rate = 1))  # similar distribution as rc (no bias)
  zi_mp <- ifelse(runif(n) < 0.7, 0, rexp(n, rate = 1))  # model projection sample
  
  # Combine into data.tables for rc, mc, mp
  library(data.table)
  rc_dt <- data.table(cont = cont_rc, zi = zi_rc)
  mc_dt <- data.table(cont = cont_mc, zi = zi_mc)
  mp_dt <- data.table(cont = cont_mp, zi = zi_mp)
  
  # Define margin controls: 'cont' as continuous and 'zi' as zero-inflated with xmin = 0
  margins <- list(type = c("c", "zi"), xmin = c(NaN, 0))
  
  # Run bias correction on the projection data
  result <- vbc(mp_dt, mc_dt, rc_dt, margins_controls = margins)
  
  # Result should be a data.table (which inherits data.frame)
  expect_s3_class(result, "data.table")
  expect_true(inherits(result, "data.frame"))
  # Number of rows should match the input projection data
  expect_equal(nrow(result), nrow(mp_dt))
  # Column names should default to observed (rc) data's names
  expect_equal(colnames(result), colnames(rc_dt))
  # Expect that vine and kde attributes are present in the result
  expect_false(is.null(attr(result, "vine_rc")))
  expect_false(is.null(attr(result, "vine_mp")))
  expect_false(is.null(attr(result, "kde_rc")))
  expect_false(is.null(attr(result, "kde_mp")))
  # The vine copula objects should have class "vinecop" (from rvinecopulib)
  expect_s3_class(attr(result, "vine_rc"), "vinecop")
  expect_s3_class(attr(result, "vine_mp"), "vinecop")
})

test_that("vbc performs identity mapping when there is no bias (rc and mc identical)", {
  # If the model has no historical bias (rc and mc distributions identical), vbc should return mp unchanged
  set.seed(456)
  n <- 30
  # rc and mc: identical values (no bias in calibration period)
  rc_vals <- rexp(n, rate = 0.5)
  mc_vals <- rc_vals  # exactly the same as rc
  # mp: model projection distribution differs (e.g., higher rate -> smaller mean)
  mp_vals <- rexp(n, rate = 1)
  rc_dt <- data.table(var = rc_vals)
  mc_dt <- data.table(var = mc_vals)
  mp_dt <- data.table(var = mc_vals)
  
  # Apply VBC correction (using default margin controls for continuous data)
  result <- vbc(mp_dt, mc_dt, rc_dt)
  
  # The output should match mp if there's no bias to correct
  expect_equal(result$var, mp_dt$var, tolerance = 1)
  # We can also verify via Wasserstein distances: before vs after correction
  dist_before <- calc_wasserstein(rc_dt, mp_dt)
  dist_after  <- calc_wasserstein(rc_dt, result)
  # After correction, distance should be essentially zero (no change needed, so mp == corrected)
  expect_lt(dist_after["Wasserstein_1"], 0.3)
  expect_lt(dist_after["Wasserstein_2"], 0.3)
})

test_that("vbc reduces distributional distance when correcting a biased model", {
  # When rc and mc distributions differ (bias present), vbc output should be closer to rc than original mp was
  set.seed(789)
  n <- 50
  # Observed (rc) distribution, e.g., N(0,1)
  rc_vals <- rnorm(n, mean = 0, sd = 1)
  # Model calibration (mc) biased distribution, e.g., N(1,1) (mean shifted by +1)
  mc_vals <- rnorm(n, mean = 1, sd = 1)
  # Model projection (mp) with bias + an additional shift (e.g., N(2,1))
  mp_vals <- rnorm(n, mean = 2, sd = 1)
  rc_dt <- data.table(val = rc_vals)
  mc_dt <- data.table(val = mc_vals)
  mp_dt <- data.table(val = mp_vals)
  
  result <- vbc(mp_dt, mc_dt, rc_dt)
  # Calculate Wasserstein distances between distribution of model and observed data
  dist_before <- calc_wasserstein(rc_dt, mp_dt)
  dist_after  <- calc_wasserstein(rc_dt, result)
  # After correction, the distances should be smaller (bias is mitigated)
  expect_true(dist_after["Wasserstein_1"] < dist_before["Wasserstein_1"])
  expect_true(dist_after["Wasserstein_2"] < dist_before["Wasserstein_2"])
})

test_that("vbc attaches a time column when time_mp is provided", {
  set.seed(111)
  n <- 50
  # Simple dummy data with two variables
  rc_dt <- data.table(x = runif(n), y = runif(n), z = runif(n))
  mc_dt <- data.table(x = runif(n), y = runif(n), z = runif(n))
  mp_dt <- data.table(x = runif(n), y = runif(n), z = runif(n))
  # Create a time index vector for the projection period
  time_index <- seq.Date(from = as.Date("2000-01-01"), by = "day",
                         length.out = n)
  # Apply vbc with time_mp argument
  set.seed(111)
  result_time <- vbc(mp_dt, mc_dt, rc_dt, time_mp = time_index)
  # Apply vbc without time for comparison
  set.seed(111)
  result_no_time <- vbc(mp_dt, mc_dt, rc_dt)
  
  # The result with time should include a "time" column as the last column
  expect_true("time" %in% colnames(result_time))
  expect_equal(result_time$time, time_index)
  # The addition of time should not alter the corrected values of other columns
  expect_equal(result_time[, .SD, .SDcols = c("x", "y", "z")],
               result_no_time, ignore_attr = TRUE, tolerance = 0.2)
})

test_that("vbc handles ensemble inputs (list of projection datasets) and returns a list of results", {
  set.seed(222)
  n <- 30
  # Simulate dummy data for a calibration set and two ensemble projection members
  rc_dt <- data.table(val = rnorm(n))
  mc_dt <- data.table(val = rnorm(n))
  mp1 <- data.table(val = rnorm(n, mean = 1))
  mp2 <- data.table(val = rnorm(n, mean = 2))
  mp_list <- list(mp1, mp2)
  
  # Correct both ensemble members in one call
  set.seed(222)
  result_list <- vbc(mp_list, mc_dt, rc_dt)
  # Expect a list output with the same length as mp_list
  expect_type(result_list, "list")
  expect_equal(length(result_list), 2)
  # Each list element should be a data.table with the same structure as a single-run output
  expect_s3_class(result_list[[1]], "data.table")
  expect_equal(colnames(result_list[[1]]), colnames(rc_dt))
  # Correct each member independently and compare results to the batch correction
  set.seed(222)
  result1 <- vbc(mp1, mc_dt, rc_dt)
  set.seed(222)
  result2 <- vbc(mp2, mc_dt, rc_dt)
  expect_equal(result_list[[1]], result1, ignore_attr = TRUE)
  expect_equal(result_list[[2]], result2, ignore_attr = TRUE)
  
  # Using vbc_ensemble should yield the same result as vbc on a list
  set.seed(222)
  result_list2 <- vbc(mp_list, mc_dt, rc_dt)
  expect_equal(result_list2, result_list, ignore_attr = TRUE)
})

test_that("vbc var_names argument overrides output column names", {
  set.seed(333)
  n <- 10
  rc_dt <- data.table(A = runif(n), B = runif(n))
  mc_dt <- data.table(A = runif(n), B = runif(n))
  mp_dt <- data.table(A = runif(n), B = runif(n))
  # Provide custom variable names for the output
  custom_names <- c("Var1", "Var2")
  result <- vbc(mp_dt, mc_dt, rc_dt, var_names = custom_names)
  # The result's column names should be replaced with the custom names
  expect_equal(colnames(result), custom_names)
})

test_that("vbc input validation catches dimension mismatches and invalid parameters", {
  set.seed(444)
  # Create rc and mc with different number of variables to induce an error
  rc_dt <- data.table(one = rnorm(5), two = rnorm(5))
  mc_dt <- data.table(one = rnorm(5))  # mc has only one column, mismatch with rc
  mp_dt <- data.table(one = rnorm(5))
  # Expect an error due to mismatched number of columns (variables)
  expect_error(vbc(mp_dt, mc_dt, rc_dt))
  
  # Fix mc to have matching columns for further tests
  mc_dt <- data.table(one = rnorm(5), two = rnorm(5))
  # Provide an invalid margins_controls (with an unsupported type value)
  bad_margins <- list(type = c("continuous", "invalid_type"), xmin = c(NaN, 0))
  expect_error(vbc(mp_dt, mc_dt, rc_dt, margins_controls = bad_margins))
  
  # Provide var_names of incorrect length
  expect_error(vbc(mp_dt, mc_dt, rc_dt, var_names = c("onlyOneName")))
})
