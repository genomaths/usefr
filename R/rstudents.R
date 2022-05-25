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

#' @rdname rstudents
#' @aliases rstudents
#' @title Studentized Residuals from a Nonlinear Model
#' @description Given a fitted model and additional information, this function
#' computes the corresponding internally studentized residuals.
#' @details If the \strong{model} argument is provided, then only the
#' argument \strong{varobj} are required. The internally studentized residuals
#' are computed as \eqn{t = residuals/(s * sqrt(1 - h))}, where \eqn{h} is the
#' diagonal of the hat matrix and \eqn{s} is the estimation of the residual
#' variation (\href{https://is.gd/2SINWt}{(see Wikipedia)}).
#'
#' If the hat matrix cannot be estimated, the standardized residual are
#' estimated as equal to the value of a residual, divided by an estimation of
#' its standard deviation, i.e., \eqn{t = residuals/s}.
#'
#' Studentized/Standardized residuals greater than 2 and less than -2 are
#' usually considered large.
#'
#' @param model An object from the ""nls" class. It can be derived. e.g., from
#' function \code{\link[stats]{nls}} or \code{\link[minpack.lm]{nlsLM}} from
#' the R packages \emph{stats} and \emph{minpack.lm}, respectively.
#' @param method Argument for \code{\link[numDeriv]{grad}} function from
#' \emph{numDeriv} package. One of "Richardson", "simple", or "complex"
#' indicating the method to use for the approximation.
#' @param varobj Objective variable used build the \strong{model}.
#' @param pars \strong{model} parameters.
#' @param fun The named expression (not a character string) of the
#' fitted model. For example, the normal distribution function is given by
#' pnorm (see ?pnorm).
#' @param residuals Residuals from the model fitting.
#' @importFrom numDeriv grad
#' @importFrom stats splinefun var
#' @export
#' @return A vector of studentized residuals when model is missing, NULL, or
#' \strong{\emph{"nls"}} class object. If \eqn{model="CDFmodel}, then it will
#' update the information carried on model$rstudent and the model will be
#' returned.

setGeneric(
    "rstudents",
    def = function(model,
    varobj,
    method = c("Richardson", "simple", "complex"),
    ...) {
        standardGeneric("rstudents")
    }
)

#' @rdname rstudents
#' @aliases rstudents
#' @export
setMethod(
    "rstudents", signature(model = "missingORNULL"),
    function(model,
    varobj,
    method = c("Richardson", "simple", "complex"),
    pars,
    fun,
    residuals) {
        method <- match.arg(method)

        #   # Derivative of the best fit function
        evalLST <- as.list(pars)
        evalLST$func <- fun
        evalLST$x <- varobj
        evalLST$method <- method

        gradient <- try(do.call(grad, evalLST),
            silent = TRUE
        )

        if (!inherits(gradient, "try-error")) {
            # NaN/missing values are replace by cubic spline interpolation
            grad.spline <- splinefun(varobj, gradient)
            ind <- which(is.na(gradient))
            gradient[ind] <- grad.spline(varobj[ind])

            rstudent <- try(
                nls.rstudent(gradient, residuals, length(pars)),
                silent = TRUE
            )

            if (inherits(gradient, "try-error")) {
                rstudent <- NA
            }
        }
        return(rstudent)
    }
)

#' @rdname rstudents
#' @aliases rstudents
#' @export
setMethod(
    "rstudents", signature(model = "nls"),
    function(model,
    varobj,
    method = c("Richardson", "simple", "complex"),
    residuals = NULL) {
        method <- match.arg(method)
        fun <- model$m$formula()[[3L]]
        if (is.null(residuals)) {
            residuals <- model$m$resid()
        }

        return(rstudents(
            varobj = varobj,
            method = method,
            pars = coef(model),
            fun = fun[[1]],
            residuals = residuals
        ))
    }
)

#' @rdname rstudents
#' @aliases rstudents
#' @importFrom stats coef
#' @export
setMethod(
    "rstudents", signature(model = "CDFmodel"),
    function(model,
    varobj,
    method = c("Richardson", "simple", "complex"),
    residuals = NULL) {
        method <- match.arg(method)
        if (is.null(residuals)) {
            residuals <- resid(model$bestfit)
        }

        cdf <- paste0("p", model$cdf)
        model$rstudent <- rstudents(
            varobj = varobj,
            method = method,
            pars = coef(model$bestfit),
            fun = cdf,
            residuals = residuals
        )
        return(model)
    }
)


#' @rdname rstudents
#' @aliases rstudents
#' @importFrom stats coef
#' @export
setMethod(
    "rstudents", signature(model = "nls.lm"),
    function(model,
    varobj,
    method = c("Richardson", "simple", "complex"),
    fun,
    residuals = NULL) {
        method <- match.arg(method)
        if (is.null(residuals)) {
            residuals <- resid(model)
        }
        return(rstudents(
            varobj = varobj,
            method = method,
            pars = coef(model),
            fun = match.fun(fun),
            residuals = residuals
        ))
    }
)


## ------------ *** Studentized residuals for nls *** #

nls.rstudent <- function(grad, residuals, num.par) {
    # grad: derivative of model as obtained from function deriv
    # R.S
    n <- length(residuals)
    s <- sqrt(n * var(residuals, na.rm = TRUE) / (n - num.par))
    h <- try(
        diag(grad %*% solve(crossprod(grad)) %*% t(grad)),
        silent = TRUE
    )
    if (!inherits(h, "try-error")) {
        res <- residuals / (s * sqrt(1 - h))
    } else {
        res <- residuals / s
        message(
            "*** Internally studentized residuals cannot be estimated\n",
            "Trying the estimation of standardized residuals: residuals/s"
        )
    }
    return(res)
}
