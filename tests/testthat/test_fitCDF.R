library(testthat)
library(usefr)

context("fitCDF tests")

test_that("fitCDF dummy test", {
    set.seed(1230)
    x1 <- rnorm(1000, mean = 0.5, sd = 1)
    cdfp <- fitCDF(x1, distNames = "Normal", only.info = TRUE, plot = FALSE)
    expect_true(cdfp$mean > 0.45)
})
