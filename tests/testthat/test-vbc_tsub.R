# vbc_tsub(): time-based subsetting via month/hour filters

test_that("vbc_tsub subsets by month and hour on single table", {
  library(data.table)
  
  set.seed(10)
  
  n   <- 24 * 15  # ~15 days 6 hourly
  tm  <- as.POSIXct("2010-01-01 00:00:00", tz = "UTC") + 3600 * seq_len(n)
  rc  <- data.table(x = rnorm(n), y = rnorm(n, 1))
  mc  <- data.table(x = rnorm(n, 2), y = rnorm(n, 3))
  mp  <- data.table(x = rnorm(n, 3), y = rnorm(n, 4))
  # Ensure each table has a POSIXct 'time' column
  if (!"time" %in% names(rc)) rc[, time := tm]
  if (!"time" %in% names(mc)) mc[, time := tm]
  if (!"time" %in% names(mp)) mp[, time := tm]
  # choose January (1) and early hours 0..6
  out <- vbc_tsub(mp, mc, rc,
                  t_subs = list(list(month = 1, hours = 0:6)))
  
  expect_s3_class(out, "data.table")
  expect_true(all(c("x", "time") %in% names(out)))
  # only hours 0..6 in January
  expect_true(all(as.integer(format(out$time, "%m")) == 1L))
  expect_true(all(as.integer(format(out$time, "%H")) %in% 0:6))
})

test_that("vbc_tsub errors on empty selection (no times match)", {
  library(data.table)
  
  set.seed(11)
  n   <- 24 * 5
  tm  <- as.POSIXct("2010-04-01 00:00:00", tz = "UTC") + 3600 * seq_len(n) # February
  rc  <- data.table(x = rnorm(n))
  mc  <- data.table(x = rnorm(n, 1))
  mp  <- data.table(x = rnorm(n, 2))
  if (!"time" %in% names(rc)) rc[, time := tm]
  if (!"time" %in% names(mc)) mc[, time := tm]
  if (!"time" %in% names(mp)) mp[, time := tm]
  
  # Ask for January hours in a February series -> empty
  expect_error(
    vbc_tsub(mp, mc, rc, t_subs = list(list(month = 1, hours = 0:6)))
  )
})
