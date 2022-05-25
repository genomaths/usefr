library(testthat)
library(usefr)

context("mcgoftest tests")

test_that("mcgoftest dummy test", {
    set.seed(1)
    x <- rnorm(1000, mean = 1.5, sd = 2)
    ks <- mcgoftest(x,
        cdf = pnorm, pars = c(1.5, 2), num.sampl = 100,
        sample.size = 100, num.cores = 1
    )
    expect_true(ks[1] > 0.1)
})
