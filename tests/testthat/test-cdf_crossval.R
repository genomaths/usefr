test_that("cdf_crossval works", {
    set.seed(1)
    x1 = rnorm(100, mean = 1.5, sd = 2)
    cdfp <- fitCDF(dt$x1, distNames = "Normal", plot = F)
    expect_true(cdf_crossval(model = cdfp$bestfit, q = x1) > 0.99)
})
