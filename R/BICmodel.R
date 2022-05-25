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

#' @rdname BICmodel
#'
#' @title Bayesian Information Criterion (BIC)
#' @description this function permits the estimation of the BIC for models for
#'     which the function 'BIC' from 'stats' packages does not work.
#' @details if for a given model 'm' BIC(m) works, then BICmodel(m) = BIC(m).
#'
#' @param model if provided, it is an R object from where the residuals and
#'     model parameters can be retrieved using resid(model) and coef(model),
#'     respectively.
#' @param residuals if provided, it is numerical vector with the residuals:
#'     residuals = observe.values - predicted.values, where predicted values
#'     are estimated from the model. If the parameter 'model' is not provided,
#'     then this parameter must be provided.
#' @param np number of model parameters. If the parameter 'model' is not
#'     provided, then 'np' and 'residuals' must be provided.
#'
#' @return BIC numerical value
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#'
#' @examples
#' set.seed(77)
#' x <- runif(100, 1, 5)
#' y <- 2 * exp(-0.5 * x) + runif(100, 0, 0.1)
#' plot(x, y)
#'
#' nlm <- nls(Y ~ a * exp(b * X),
#'     data = data.frame(X = x, Y = y),
#'     start = list(a = 1.5, b = -0.7),
#'     control = nls.control(maxiter = 10^4, tol = 1e-05),
#'     algorithm = "port"
#' )
#' ## The estimations of Akaike information criteria given by BIC' function
#' ## from stats' R package and from 'AICmodel' function are equals.
#' BICmodel(nlm) == BIC(nlm)
#'
#' ## Now, using residuals from the fitted model:
#' res <- y - coef(nlm)[1] * exp(coef(nlm)[2] * x)
#'
#' BICmodel(residuals = res, np = 2) == BIC(nlm)
#'
#' @importFrom stats resid coef
#' @seealso \code{\link{AICmodel}}
#' @export
BICmodel <- function(model = NULL, residuals = NULL, np = NULL) {
    if (is.null(model) && is.null(residuals)) {
        stop(paste(
            "At least one of the parameter 'model'",
            " or 'residual' must be provided"
        ))
    }
    if (!is.null(model)) {
        if (is.null(np)) np <- length(coef(model))
        RESID <- resid(model)
    }
    if (is.null(np)) {
        stop("The number of model parameter 'np' must be provided")
    }
    if (!is.null(residuals)) RESID <- residuals

    sse <- sum(RESID^2, na.rm = TRUE)
    N <- length(RESID)
    return(N * (1 + log(2 * pi) + log(sse / N)) + log(N) * (1L + np))
}
