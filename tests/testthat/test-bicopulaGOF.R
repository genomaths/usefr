test_that("bicopulaGOF works", {
    set.seed(12)
    margins <- c("norm", "norm")
    X <- rnorm(2 * 1e3, mean = 0, sd = 10)
    Y <- rnorm(2 * 1e3, mean = 0, sd = 10)
    parMargins <- list(
        list(mean = 0, sd = 10),
        list(mean = 0, sd = 10)
    )
    gof <- bicopulaGOF(
        x = X, y = Y, copula = "normalCopula", sample.size = 1e2,
        margins = margins, paramMargins = parMargins, nboots = 99,
        Rosenblatt = TRUE, approach = "adgamma", num.cores = 1L,
        seed = 12
    )
    expect_true(gof$gof[2] > 0.9)
})
