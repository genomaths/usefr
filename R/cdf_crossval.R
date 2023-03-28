## ########################################################################## #
##
## Copyright (C) 2019 Robersy Sanchez <https://genomaths.com/>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.
##
## ########################################################################## #
##

#' @rdname cdf_crossval
#' @title Cross validation of Cumulative Distribution Function model
#' @description This function returns A goodness-of-fit criteria for nonlinear
#' model selection, specifically, the cross-validation correlation coefficient R
#' (R.Cross.val).
#' @details The cross-validation correlation coefficient R (R.Cross.val) is
#' an estimator of the average cross-validation predictive power (1).
#' @param formula No required for when a model from class
#' \strong{\emph{CDFmodel}} or \code{\link[stats]{nls}} is provided. Otherwise,
#' it will be nonlinear model formula including variables and parameters,
#' which will be coerced to a formula if necessary. For example, for a Gamma
#' model the formula will be: "Y ~ pgamma(q, shape, scale)", where
#' \code{\link[stats]{pgamma}} function is available in 'stats' R package.
#' However, the \code{\link[minpack.lm]{nls.lm}} class model created by
#' function \code{\link{fitCDF}} has incorporated the formula information.
#' @param pars Estimated model parameters.
#' @param q Objective variable used to build the model, typically called a
#' vector of quantiles. The model's formula must be expressed in terms of
#' variable \eqn{'q'}.
#' @param min.val A number denoting the lower bound of the domain where CDF
#' is defined. For example, for Weibull and GGamma \strong{\emph{min.val = 0}}.
#' @param logx logical(1). If TRUE, then a logarithm transformation will be
#' applied: \eqn{q = log1p(q)}.
#' @param loss.fun Described in \code{\link{fitCDF}}.
#' @param maxiter,ptol,maxfev Arguments for function
#' \code{\link[minpack.lm]{nlsLM}} and/or \code{\link[minpack.lm]{nls.lm}}.
#' @param minFactor A positive numeric value specifying the minimum step-size
#' factor allowed on any step in the iteration. The increment is calculated
#' with a Gauss-Newton algorithm and successively halved until the residual sum
#' of squares has been decreased or until the step-size factor has been reduced
#' below this limit. Default value: 10^-6.
#' @importFrom minpack.lm nlsLM
#' @export
#' @seealso \code{\link{mcgoftest}} for Bootstrap test for Goodness of fit.
#' @references
#'   \enumerate{
#'     \item Stevens JP. Applied Multivariate Statistics for the Social
#'           Sciences. Fifth Edit. Routledge Academic; 2009.
#'  }
#' @aliases cdf_crossval
#' @examples
#' ## Let's simulate a sample from normal distribution
#' x1 = rnorm(10000, mean = 1.5, sd = 2) + runif(10^4)
#'
#' ## Let's build a model
#' cdfp <- fitCDF(x1, distNames = "Normal", plot = F)
#'
#' ## Next, we get an estimation of the cross-validation correlation
#' ## coefficient R (R.Cross.val)
#' cdf_crossval(model = cdfp$bestfit, q = x1)
setGeneric(
    "cdf_crossval",
    def = function(model,
    ...) {
        standardGeneric("cdf_crossval")
    }
)

setOldClass(c("CDFmodel", "CDFmodelList"))
setOldClass("nls")
setOldClass("nls.lm")
setOldClass("nlsModel")

setClassUnion("missingORNULL", c("missing", "NULL"))

#' @rdname cdf_crossval
#' @aliases cdf_crossval
#' @export
setMethod(
    "cdf_crossval", signature(model = "missingORNULL"),
    function(
        model,
        formula,
        pars,
        q,
        logx = FALSE,
        min.val = NULL,
        loss.fun = c(
            "linear", "huber", "smooth",
            "cauchy", "arctg"
        ),
        maxiter = 1024,
        maxfev = 1e5,
        ptol = 1e-12,
        minFactor = 1e-6) {
        if (!inherits(formula, "formula")) {
            stop("*** Argument for formula must be a 'formula' class object")
        }

        if (!is.null(min.val)) {
            q <- q[ which(q > min.val) ]
        }

        if (logx) {
            q <- q[ q > 0 ]
            q <- log(q)
        }
        n <- length(q)
        cros.ind.1 <- sample.int(n, size = round(n / 2))
        cros.ind.2 <- setdiff(seq_len(n), cros.ind.1)

        Fy <- ecdf(q)
        pX <- Fy(q)

        FIT1 <- try(nlsLM(
            formula,
            data = data.frame(
                q = q[cros.ind.1],
                Y = pX[cros.ind.1]
            ),
            start = as.list(pars),
            control = list(maxiter = maxiter, ptol = ptol)
        ),
        silent = TRUE
        )

        if (inherits(FIT1, "try-error")) {
            probfun <- match.fun(formula[[3L]][[1L]])
            FIT1 <- try(nls.lm(
                par = pars,
                fn = optFun,
                probfun = probfun,
                quantiles = q[cros.ind.1],
                prob = pX[cros.ind.1],
                loss.fun = loss.fun,
                control = nls.lm.control(
                    maxiter = maxiter,
                    maxfev = maxfev,
                    ptol = ptol
                )
            ),
            silent = TRUE
            )
        }

        FIT2 <- try(nlsLM(
            formula,
            data = data.frame(
                q = q[cros.ind.2],
                Y = pX[cros.ind.2]
            ),
            start = as.list(pars),
            control = list(maxiter = maxiter, ptol = ptol)
        ),
        silent = TRUE
        )

        if (inherits(FIT2, "try-error")) {
            probfun <- match.fun(formula[[3L]][[1L]])
            FIT2 <- try(nls.lm(
                par = pars,
                fn = optFun,
                probfun = probfun,
                quantiles = q[cros.ind.2],
                prob = pX[cros.ind.2],
                loss.fun = loss.fun,
                control = nls.lm.control(
                    maxiter = maxiter,
                    maxfev = maxfev,
                    ptol = ptol
                )
            ),
            silent = TRUE
            )
        }

        if (inherits(FIT1, "try-error") || inherits(FIT2, "try-error")) {
            R.Cross.val <- NA
            message(
                "\nError!\n",
                "*** Model fiting fails. R.Cross.val cannot be estimated "
            )
            if (inherits(FIT1, "try-error")) {
                message(FIT1)
            } else {
                message(FIT2)
            }
        } else {
            fun <- as.character(formula[[3L]][[1L]])

            ## prediction using model 1
            evalLIST <- as.list(coef(FIT1)[seq_along(pars)])
            evalLIST$q <- q[cros.ind.2]
            p.FIT1 <- do.call(fun, evalLIST)
            if (all(is.na(p.FIT1)) || all(is.na(pX[cros.ind.2]))) {
                R.FIT1 <- 0
            }
            else {
                R.FIT1 <- cor(p.FIT1, pX[cros.ind.2], use = "complete.obs")
            }


            ## prediction using model 2
            evalLIST <- as.list(coef(FIT2)[seq_along(pars)])
            evalLIST$q <- q[cros.ind.1]
            p.FIT2 <- do.call(fun, evalLIST)

            if (all(is.na(p.FIT2)) || all(is.na(pX[cros.ind.1]))) {
                R.FIT2 <- 0
            }
            else {
                R.FIT2 <- cor(p.FIT2, pX[cros.ind.1], use = "complete.obs")
            }

            R.Cross.val <- (length(p.FIT1) * R.FIT1 +
                length(p.FIT2) * R.FIT2) / (length(p.FIT1) +
                length(p.FIT2))
        }
        return(c(R.Cross.val = R.Cross.val))
    }
)


#' @rdname cdf_crossval
#' @aliases cdf_crossval
#' @export
setMethod(
    "cdf_crossval", signature(model = "nls"),
    function(model,
    q,
    logx = FALSE,
    min.val = NULL,
    maxiter = 1024,
    ptol = 1e-12,
    minFactor = 1e-6) {
        cdf_crossval(
            formula = model$m$formula(),
            pars = coef(model),
            q = q,
            logx = logx,
            min.val = min.val,
            maxiter = maxiter,
            ptol = ptol,
            minFactor = minFactor
        )
    }
)

#' @rdname cdf_crossval
#' @aliases cdf_crossval
#' @export
setMethod(
    "cdf_crossval", signature(model = "CDFmodel"),
    function(model,
    q,
    logx = FALSE,
    min.val = NULL,
    maxiter = 1024,
    ptol = 1e-12,
    minFactor = 1e-6) {
        formula <- model$formula

        cdf_crossval(
            formula = formula,
            pars = coef(model$bestfit),
            q = q,
            logx = logx,
            min.val = min.val,
            maxiter = maxiter,
            ptol = ptol,
            minFactor = minFactor
        )
    }
)

#' @rdname cdf_crossval
#' @aliases cdf_crossval
#' @export
setMethod(
    "cdf_crossval", signature(model = "nls.lm"),
    function(model,
    formula,
    q,
    logx = FALSE,
    min.val = NULL,
    maxiter = 1024,
    ptol = 1e-12,
    minFactor = 1e-6) {
        if (missing(formula) && is.null(model$formula))
            stop("*** An argument for 'formula' must be provided.")

        if (missing(formula) && !is.null(model$formula))
            formula <- model$formula

        cdf_crossval(
            formula = formula,
            pars = coef(model),
            q = q,
            logx = logx,
            min.val = min.val,
            maxiter = maxiter,
            ptol = ptol,
            minFactor = minFactor
        )
    }
)
