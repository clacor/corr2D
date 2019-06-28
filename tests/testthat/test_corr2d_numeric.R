context("corr2d numeric correctness")
library("corr2D")

test_that("corr2d gives correct numeric result", {
  expect_equal(corr2d(matrix(c(1:3, 3:1), 3, 2), Wave1 = 1:2, Ref1 = rep(0, 2), corenumber = 1)$FT[1, 2], 30/(2*pi) + 0i)
  expect_equal(corr2d(matrix(c(1:3, 1:3), 3, 2), Wave1 = 1:2, corenumber = 1)$FT[1, 2], 6/(2*pi) + 0i)
  expect_equal(corr2d(matrix(c(1:3, 1:3), 3, 2) * 10, Wave1 = 1:2, corenumber = 1)$FT[1, 2], corr2d(matrix(c(1:3, 1:3), 3, 2), Wave1 = 1:2, corenumber = 1)$FT[1, 2] * 100)
})