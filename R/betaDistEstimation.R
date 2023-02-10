## Copyright (C) 2019 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez
##
## This program is part 'usefr' R package.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## 'usefr' is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

#' @rdname betaDistEstimation
#'
#' @title Nonlinear Parameter Estimation for Beta Distribution
#' @description This function perform a rough estimation of the shape
#'     parameters of beta distribution
#' @details In order to obtain the estimates for shape parameters of beta
#' distribution, the squared of the difference between the empirical cumulative
#' distribution function (ecdf) & the theoretical cdf is minimized using the
#' Non-Linear Minimization function \code{\link[stats]{nlm}} 'stats' package.
#'
#' If \code{\link[stats]{nlm}} function fails, then an estimation using
#' \code{\link[stats]{optim}} function is tried.
#'
#' @param q prior probabilities
#' @param init.pars  initial parameter values. Defaults to alpha = 1 &
#' beta = 1, which imply the parsimony pseudo-counts greater than zero.
#' @param force.optim Whether to force the use of \code{\link[stats]{optim}}
#' function for the parameter estimation. Default is FALSE.
#' @param hessian if TRUE, the hessian of f at the minimum is returned.
#' @param method,control,lower,upper,gr (Optional). In the case that
#' \code{\link[stats]{nlm}} function fails, the methods, list of control
#' parameters, and bounders to be used with \code{\link[stats]{optim}} to
#' accomplish the parameter estimation.
#' @param control (Optional). In the case that \code{\link[stats]{nlm}}
#' function fails, a list of control parameters to be used with
#' \code{\link[stats]{optim}} function (see function help: ?optim) accomplish
#' the parameter estimation.
#' @param ... Further parameter for \code{\link[stats]{nlm}} function.
#' @return A list with components, which would vary depending on whether the
#' estimation was performed with \code{\link[stats]{nlm}} or
#' \code{\link[stats]{optim}}. In all the cases the list element carrying
#' the estimated parameters values is named \strong{\emph{parameters}}.
#'
#' @importFrom stats ecdf optim nlm pbeta rbeta var pt
#' @author Robersy Sanchez <https://genomaths.com>
#' @export
#' @seealso \code{\link{betaBinPost}} and \code{\link{estimateDirichDist}}
#' @examples
#' ### A random generation numerical values with Beta distribution
#' x1 <- rbeta(n = 1000, shape1 = 2, shape2 = 3)
#'
#' ### Parameter estimation with "nlm" function
#' betaDistEstimation(q = x1, gradtol = 1e-12, hessian = TRUE)
#'
#' ### Parameter estimation with "optim" function
#' betaDistEstimation(q = x1, force.optim = TRUE, hessian = TRUE)
#'
betaDistEstimation <- function(q,
    init.pars = c(1, 1),
    force.optim = FALSE,
    hessian = FALSE,
    method = "BFGS",
    gr = NULL,
    control = list(
        maxit = 500,
        abstol = (10^-8)
    ),
    lower = -Inf,
    upper = Inf,
    seed = 123,
    ...) {
    ## q: prior probabilities
    ## init.pars: initial parameter values. Defaults to alpha = 1 & beta = 1,
    ## which the parsimony pseudo-counts greater than zero.

    Q <- ecdf(q) ## Empirical Cumulative Distribution Function
    dat <- data.frame(Q = Q(q), q = q)

    #### Beta fitting. In order to obtain the estimates for
    #### shape paramaters the squared of the difference
    #### between the ecdf & the theoretical cdf is
    #### minimized
    min.RSS <- function(data, inits) {
        alpha <- inits[1] ## This corresponds to the initial
        ## starting parameter for alpha
        beta <- inits[2] ## This corresponds to the initial
        ## starting parameter for beta Because optim
        ## minimizes a function, the
        with(data, sum((Q - pbeta(q, alpha, beta))^2))
    }

    if (!force.optim) {
        fit <- try(suppressWarnings(
            nlm(
                f = min.RSS,
                p = init.pars,
                data = dat,
                hessian = hessian,
                ...
            )
        ),
        silent = TRUE
        )

        if (!inherits(fit, "try-error")) {
            nms <- names(fit)
            nms[nms == "estimate"] <- "estimate"
            names(fit) <- nms
            fit$opt.fun <- "nlm"
        }
        else
            force.optim <- TRUE
    }

    if (force.optim) {
        fit <- try(suppressWarnings(optim(
            par = init.pars,
            fn = min.RSS,
            gr = gr,
            data = dat,
            method = method,
            control = control,
            hessian = hessian,
            lower = lower,
            upper = upper
        )),
        silent = TRUE
        )

        if (!inherits(fit, "try-error")) {
            names(fit) <- c(
                "estimate", "value", "counts",
                "convergence", "message"
            )
        }
    }

    ## ========================= Cross validation =========================

    set.seed(seed)
    l <- length(q)
    cros.ind.1 <- sample.int(l, size = round(l / 2))
    cros.ind.2 <- setdiff(1:l, cros.ind.1)

    x1 <- q[cros.ind.1]
    x2 <- q[cros.ind.2]
    y1 <- Q(x1)
    y2 <- Q(x2)

    FIT1 <- try(suppressWarnings(
        nlm(
            f = min.RSS,
            p = init.pars,
            data = data.frame(Q = y1, q = x1),
            hessian = hessian,
            ...
        )
    ),
    silent = TRUE
    )

    FIT2 <- try(suppressWarnings(
        nlm(
            f = min.RSS,
            p = init.pars,
            data = data.frame(Q = y2, q = x2),
            hessian = hessian,
            ...
        )
    ),
    silent = TRUE
    )

    if (inherits(FIT1, "try-error") || inherits(FIT2, "try-error")) {
        R.cross.FIT <- 0
    } else {
        ## prediction using model 1
        p.FIT1 <- pbeta(
            q = x2,
            shape1 = FIT1$estimate[1],
            shape2 = FIT1$estimate[2]
        )
        R.FIT1 <- try(cor(p.FIT1, y2, use = "complete.obs"), silent = TRUE)
        if (inherits(R.FIT1, "try-error")) {
            R.FIT1 <- try(cor(p.FIT1, y2, use = "pairwise.complete.obs"),
                silent = TRUE
            )
        }
        ## prediction using model 2
        p.FIT2 <- pbeta(
            q = x1,
            shape1 = FIT2$estimate[1],
            shape2 = FIT2$estimate[2]
        )
        R.FIT2 <- try(cor(p.FIT2, y1, use = "complete.obs"), silent = TRUE)
        if (inherits(R.FIT2, "try-error")) {
            R.FIT2 <- try(cor(p.FIT2, y1, use = "pairwise.complete.obs"),
                silent = TRUE
            )
        }

        if (inherits(R.FIT1, "try-error") && inherits(R.FIT2, "try-error")) {
            R.cross.FIT <- 0
        } else {
            term1 <- length(p.FIT1) * R.FIT1
            term2 <- length(p.FIT2) * R.FIT2
            R.cross.FIT <- (term1 + term2) / (length(p.FIT1) + length(p.FIT2))
        }
    }


    ## ====================== Adjusted R^2 =============================
    fitted <- pbeta(
        q = q,
        shape1 = fit$estimate[1],
        shape2 = fit$estimate[2]
    )
    residuals <- dat$Q - fitted
    N <- length(fitted)
    aic <- AICmodel(residuals = residuals, np = 2)
    bic <- BICmodel(residuals = residuals, np = 2)

    X <- cbind(1, q)
    betaHat <- solve(t(X) %*% X) %*% t(X) %*% q
    var_betaHat <- var(q) * solve(t(X) %*% X)
    se_beta <- sqrt(diag(var_betaHat))
    t_value <- fit$estimate[1:2] / se_beta
    wald_test <- pt(q = t_value, df = N - 2, lower.tail = FALSE)
    wald_test[wald_test < 1e-16] <- "<1e-16"


    if (!inherits(fit, "try-error")) {
        ## **** R squares ****
        Adj.R.Square <- (1 - sum(residuals^2) / ((N - 2)) *
            var(q, use = "everything"))

        Adj.R.Square <- ifelse(is.na(Adj.R.Square) || Adj.R.Square < 0,
            0, Adj.R.Square
        )

        ## Stein adjusted R square
        rho <- (1 - ((N - 2) / (N - 3)) * ((N + 1) / (N)) *
            (1 - Adj.R.Square))
        rho <- ifelse(is.na(rho) | rho < 0, 0, rho)
    }

    stats <- data.frame(
        Estimate = c(
            shape1 = fit$estimate[1],
            shape2 = fit$estimate[2]
        ),
        Std.Error = se_beta,
        t_value = t_value,
        "Pr(>|t|)" = wald_test,
        Adj.R.Square = c(Adj.R.Square, NA),
        rho = c(rho, NA),
        R.Cross.val = c(R.cross.FIT, NA),
        AIC = c(aic, NA),
        BIC = c(bic, NA),
        n = c(N, NA)
    )
    names(stats) <- c(
        "Estimate", "Std.Error", "t_value", "Pr(>|t|)",
        "Adj.R.Square", "rho", "R.Cross.val", "AIC",
        "BIC", "n"
    )

    stats <- structure(stats, class = c("BetaModel", "data.frame"))
    return(stats)
}
