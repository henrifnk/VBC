test_that("data", {
  expect_no_error(data("climate"))
  expect_list(climate, types = "data.frame", len = 4)
})
