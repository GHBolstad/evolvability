test_that("round_and_format rounds to a fixed number of decimals", {
  expect_equal(round_and_format(c(1.234, 5.6789)), c("1.23", "5.68"))
  expect_equal(round_and_format(1.5, digits = 3), "1.500")
})

test_that("round_and_format keeps significant digits", {
  expect_equal(round_and_format(0.123, sign_digits = 2), "0.12")
  expect_equal(round_and_format(0.0456, sign_digits = 2), "0.046")
  expect_equal(round_and_format(0.5, sign_digits = 2), "0.50")
  expect_equal(round_and_format(1.5, sign_digits = 2), "1.5")
  expect_equal(round_and_format(12.34, sign_digits = 3), "12.3")
  expect_equal(round_and_format(1234.5, sign_digits = 2), "1200")
  expect_equal(round_and_format(-0.123, sign_digits = 2), "-0.12")
})

test_that("round_and_format with sign_digits is vectorized", {
  expect_equal(round_and_format(c(0.123, 1.5, 100, 0), sign_digits = 2),
               c("0.12", "1.5", "100", "0"))
  expect_equal(round_and_format(c(NA, 0.5), sign_digits = 2), c("NA", "0.50"))
})
