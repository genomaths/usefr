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
#' computes the corresponding studentized residuals.
#' @details If the \strong{model} argument is provided, then only the
#' argument \strong{varobj} are required.
#'
#' @param model An object from the ""nls" class. It can be derived. e.g., from
#' function \code{\link[stats]{nls}} or \code{\link[minpack.lm]{nlsLM}} from
#' the R packages \emph{stats} and \emph{minpack.lm}, respectively.
#' @param method Argument for \code{\link[numDeriv]{grad}} function from
#' \emph{numDeriv} package. One of "Richardson", "simple", or "complex"
#' indicating the method to use for the approximation.
#' @param varobj Objective variable used build the \strong{model}.
#' @param par \strong{model} parameters.
#' @param fun The named expression (not a character string) of the
#' fitted model. For example, the normal distribution function is given by
#' pnorm (see ?pnorm).
#' @param residuals Residuals from the model fitting.
#' @importFrom numDeriv grad
#' @importFrom stats splinefun var
#' @export
#' @return A vector of studentized residuals.

setGeneric(
    "rstudents",
    def = function(
                    model,
                    varobj,
                    method = c("Richardson",  "simple", "complex"),
                    ...) standardGeneric("rstudents"))

setOldClass(c("nls", "CDFmodel"))

setClassUnion("missingORNULL", c("missing", "NULL"))

#' @rdname rstudents
#' @aliases rstudents
#' @export
setMethod("rstudents", signature(model = "missingORNULL"),
    function(
            model,
            varobj,
            method = c("Richardson",  "simple", "complex"),
            par,
            fun,
            residuals) {

        method = match.arg(method)

        #   # Derivative of the best fit function
        evalLST <- as.list(par)
        evalLST$func <- fun
        evalLST$x <- varobj
        evalLST$method <- method

        gradient <- try(do.call(grad, evalLST),
                        silent = TRUE)

        if (!inherits(gradient, "try-error")) {
            # NaN/missing values are replace by cubic spline interpolation
            grad.spline = splinefun(varobj, gradient)
            ind = which(is.na(gradient))
            gradient[ind] <- grad.spline(varobj[ ind ])

            rstudent <- try(
                nls.rstudent(gradient, residuals, length(fit$par)),
                silent = TRUE)

            if (inherits(gradient, "try-error"))
                rstudent <- NA
        }
        return(rstudent)
    })

#' @rdname rstudents
#' @aliases rstudents
#' @export
setMethod("rstudents", signature(model = "nls"),
    function(
            model,
            varobj,
            method = c("Richardson",  "simple", "complex")) {

    method <- match.arg(method)
    fun = model$m$formula()[[ 3L ]]

    return(rstudents(
                    varobj = varobj,
                    par = coef(model),
                    fun = fun[[1]],
                    residuals = model$m$resid(),
                    method = method))
    }
)



## ------------ *** Studentized residuals for nls *** #

nls.rstudent <- function(grad, residuals, num.par) {
    # grad: derivative of model as obtained from function deriv
    # R.S
    n = length(residuals)
    h = diag(grad %*% solve(crossprod(grad)) %*% t(grad))
    s = sqrt(n * var(residuals, na.rm = TRUE) / (n - num.par))
    residuals / (s * sqrt(1 - h))
}
