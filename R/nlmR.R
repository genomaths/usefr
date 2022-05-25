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
#' @rdname nlmR
#' @title Nonlinear Goodness-of-fit and Determination Coefficients
#' @description This function returns goodness-of-fit criteria for nonlinear
#' model selection. Adjusted R squared (rho), the AIC, and BIC are estimated.
#' @details The Stein's formula for adjusted R squared (rho) is applied as an
#' estimator of the average cross-validation predictive power (1).
#' @param object An object which inherits from 'nls' class, typically returned
#' after performing a nonlinear regression fit with function
#' \code{\link[stats]{nls}} or \code{\link[minpack.lm]{nlsLM}}.
#' @export
#' @seealso \code{\link{mcgoftest}} for Bootstrap test for Goodness of fit.
#' @examples
#' ### Examples from 'nls' doc
#' DNase1 <- subset(DNase, Run == 1)
#'
#' ## using a selfStart model
#' fm1DNase1 <- nls(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
#' nlmR(fm1DNase1)
#'
#' @references
#' 1. Stevens JP. Applied Multivariate Statistics for the Social
#'     Sciences. Fifth Edit. Routledge Academic; 2009.
#' @author Robersy Sanchez - 10/22/2020

nlmR <- function(object) {
    if (!inherits(object, "nls")) {
        stop(
            "\n*** 'object' must inherits from class 'nls', i.e.,\n",
            "an object derived from nonlinear fit from fron functions: \n",
            "'nls' or 'nlsLM' from packages 'stats' or 'minpack.lm'."
        )
    }

    m <- length(coef(object))
    n <- length(object$m$lhs())
    ## **** R squares ****
    Adj.R.Square <- (1 - (deviance(object) / ((n - m) *
        var(object$m$lhs(),
            use = "everything"
        ))))
    Adj.R.Square <- ifelse(is.na(Adj.R.Square) ||
        Adj.R.Square < 0, 0, Adj.R.Square)
    ## Stein adjusted R square
    if (m > 2) {
        rho <- ((n - 1) / (n - 4)) * ((n - 2) / (n - 5)) * ((n + 1) / n)
    } else {
        rho <- ((n - 1) / (n - 3)) * ((n - 2) / (n - 4)) * ((n + 1) / n)
    }
    rho <- 1 - rho * (1 - Adj.R.Square)
    rho <- ifelse(is.na(rho) | rho < 0, 0, rho)
    res <- data.frame(
        Adj.R.Square = Adj.R.Square, rho = rho,
        AIC = AICmodel(object), BIC = BICmodel(object)
    )
    return(res)
}
