## Copyright (C) 2019 Robersy Sanchez <https://genomaths.com/>
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

#' @name weibull3p
#' @aliases weibull3p
#' @aliases rweibull3p
#' @aliases pweibull3p
#' @aliases dweibull3p
#' @aliases qweibull3p
#'
#' @title Weibull Distribution with Three Parameters
#' @description Probability density function (PDF), cummulative density function
#'     (CDF), quantile function and random generation for the Three-Parameter
#'     Weibull (weibull3p) distribution.
#' @details It is important to mention that Weibull distribution is defined for
#'     \eqn{x > 0} and, therefore, values \eqn{x < 0} will yield probability
#'     zero. While to obtain negative random values with function
#'     \strong{rweibull3p} is an indication that the location parameter
#'     \eqn{\mu} has not a valid value.
#' @param x,q numeric vector
#' @param n number of observations
#' @param shape shape parameter, or slope, defaulting to 1
#' @param scale scale parameter, or characteristic life,  defaulting to 1
#' @param mu location parameter (numerical, default 0).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X<=x],
#'     otherwise, P[X > x]
#' @param log.p logical; if TRUE, probabilities/densities p are returned as
#'     log(p).
#' @return Three-parameter PDF values for dweibull3p, three-parameter
#'     probability for pweibull3p, quantiles or three-parameter random
#'     generated values for rweibull3p.
#' @examples
#' set.seed(123) # set a seed for random generation
#' ##  Sampling from the specified Weibull distribution
#' r <- rweibull3p(n = 1e5, shape = 1.4, scale = 3.7, mu = 0.9)
#'
#' ## The histogram and density
#' hist(r, 100, freq = FALSE, xlim = c(0, 20))
#' r1 <- seq(0, 20, by = 0.1)
#' lines(r1, dweibull3p(r1, shape = 1.4, scale = 3.7, mu = 0.9),
#'     col = "red"
#' )
#'
#' @aliases dweibull3p
#' @rdname weibull3p
#' @title Weibull Distribution with Three Parameters
#' @export
dweibull3p <- function(x, shape = 1, scale = 1, mu = 0, log = FALSE) {
    d <- dweibull(x = x - mu, shape = shape, scale = scale, log = log)
    return(d)
}

#' @aliases pweibull3p
#' @rdname weibull3p
#' @title Weibull Distribution with Three Parameters
#' @export
pweibull3p <- function(q, shape = 1, scale = 1, mu = 0, lower.tail = TRUE,
    log.p = FALSE) {
    p <- pweibull(
        q = q - mu, shape = shape, scale = scale,
        lower.tail = lower.tail, log.p = log.p
    )
    return(p)
}

#' @aliases qweibull3p
#' @rdname weibull3p
#' @title Weibull Distribution with Three Parameters
#' @export
qweibull3p <- function(p, shape = 1, scale = 1, mu = 0, lower.tail = TRUE,
    log.p = FALSE) {
    q <- qweibull(
        p = p, shape = shape, scale = scale, lower.tail = lower.tail,
        log.p = log.p
    ) + mu
    return(q)
}

#' @aliases rweibull3p
#' @rdname weibull3p
#' @title Weibull Distribution with Three Parameters
#' @export
rweibull3p <- function(n, shape = 1, scale = 1, mu = 0) {
    r <- rweibull(n = n, shape = shape, scale = scale) + mu
    return(r)
}
