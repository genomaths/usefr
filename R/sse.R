## ########################################################################## #
##
## Copyright (C) 2019-2025 Robersy Sanchez <https://genomaths.com/>
##
## Author: Robersy Sanchez
#
## This file is part of the R package "usefr".
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
##
## ########################################################################## #

#' @name sse
#' @aliases sse
#' @title Sum of Squares of Error Statistics
#' @description
#' Computes the residual sum of squares and the root of the mean of for squares
#' for models: \code{\link[stats]{nls}}, \code{\link[minpack.lm]{nls.lm}}, and
#' models CDFmodel and CDFmodelList obatained with [fitCDF].
#'
#' For a \code{\link[stats]{nls}} the RMSE is computed as:
#' \deqn{ y = resid(model)}
#' \deqn{RMSE = sqrt(mean(y^2))}
#'
#' @param model An object from one of the classes: \code{\link[stats]{nls}},
#' \code{\link[minpack.lm]{nls.lm}}, CDFmodel or CDFmodelList, which are
#' obtained after running [fitCDF].
#' @param stat Select which statistic to return: 1) The root mean of the square
#' of error or the square errors of residuals.
#' @returns The requested statistic for each model.
#' @examples
#' ## Simulate data from Normal distribution
#' ## and search for the best fitted model
#' x <-  rnorm(1000, mean = 5, sd = 2)
#' models <- fitCDF2(x, distNames = c(1, 7, 11, 12))
#'
#' ## Bases on the AIC criteria the best model is Weibull 3P
#' models
#'
#' ## In the current case Weibull is also the model with the lowest error
#' sse(model)
#'
#' @export
sse <- function(model, stat) {
    UseMethod("sse")
}

#' @rdname sse
#' @aliases sse
#' @export
sse.nls <- function(model, stat = c("rmse", "sse")) {

        stat <- match.arg(stat)
        y <- resid(model)

        switch(stat,
            "rmse" = sqrt(mean(y^2)),
            "sse" = sum(y^2)
        )
}


#' @rdname sse
#' @aliases sse
#' @export
sse.nls.lm <- function(model, stat = c("rmse", "sse")) {

    stat <- match.arg(stat)
    y <- resid(model)

    switch(stat,
           "rmse" = sqrt(mean(y^2)),
           "sse" = sum(y^2)
    )
}


#' @rdname sse
#' @aliases sse
#' @export
sse.CDFmodel <- function(model, stat = c("rmse", "sse")) {
    stat <- match.arg(stat)
    sapply(model$fit, function(y) {
        sse(model = y, stat = stat)
    })
}


#' @rdname sse
#' @aliases sse
#' @export
sse.CDFmodelList <- function(model, stat = c("rmse", "sse")) {
        stat <- match.arg(stat)
        n <- length(model)
        m <- model[-n]
        sapply(m, function(x) {
            sapply(x$fit, function(y) {
                sse(model = y, stat = stat)
            })
    })
}
