test_that("BICmodel works", {
    set.seed(77)
    x <- runif(100, 1, 5)
    y <- 2 * exp(-0.5 * x)
    nlm <- nls(Y ~ a * exp(b * X),
               data = data.frame(X = x, Y = y),
               start = list(a = 1.5, b = -0.7),
               control = nls.control(maxiter = 10^4, tol = 1e-05),
               algorithm = "port"
    )
    expect_true(round(BICmodel(nlm), 3) == round(BIC(nlm), 3))
})
