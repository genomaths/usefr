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

#' @name gamma3p
#' @aliases gamma3p
#' @aliases dgamma3p
#' @aliases pgamma3p
#' @aliases qgamma3p
#' @aliases rgamma3p
#'
#' @title Gamma Distribution with Three Parameters
#' @description Probability density function (PDF), cummulative density function
#'     (CDF), quantile function and random generation for the Three-Parameter
#'     Gamma (gamma3p) distribution.
#' @details It is important to mention that Gamma distribution is defined for
#'     \eqn{x > 0} and, therefore, values \eqn{x < 0} will yield probability
#'     zero. While to obtain negative random values with function
#'     \strong{rgamma3p} is an indication that the location parameter
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
#' @return Three-parameter PDF values for dgamma3p, three-parameter
#'     probability for pgamma3p, quantiles or three-parameter random
#'     generated values for rgamma3p.
#' @examples
#' set.seed(123) # set a seed for random generation
#' ##  Sampling from the specified Gamma distribution
#' r <- rgamma3p(n = 1e5, shape = 1.4, scale = 3.7, mu = 0.9)
#'
#' ## The histogram and density
#' hist(r, 100, freq = FALSE, xlim = c(0, 20))
#' r1 <- seq(0, 20, by = 0.1)
#' lines(r1, dgamma3p(r1, shape = 1.4, scale = 3.7, mu = 0.9),
#'     col = "red")
#'
#' @aliases dgamma3p
#' @rdname gamma3p
#' @title Gamma Distribution with Three Parameters
#' @export
dgamma3p <- function(x, shape = 1, rate = 1, scale = 1/rate, mu = 0, log = FALSE) {
  d <- dgamma(x = x - mu, shape = shape, scale = scale, log = log)
  return(d)
}

#' @aliases pgamma3p
#' @rdname gamma3p
#' @title Gamma Distribution with Three Parameters
#' @importFrom stats pgamma
#' @export
pgamma3p <- function(q, shape = 1,  rate = 1, scale = 1/rate, mu = 0, lower.tail = TRUE,
                      log.p = FALSE) {
   p <- pgamma(q = q - mu, shape = shape, scale = scale,
                   lower.tail = lower.tail, log.p = log.p)
   return(p)
}

#' @name qgamma3p
#' @rdname gamma3p
#' @title Gamma Distribution with Three Parameters
#' @export
qgamma3p <- function(p, shape = 1, rate = 1, scale = 1/rate, mu = 0, lower.tail = TRUE,
                      log.p = FALSE) {
   q <- qgamma(p = p, shape = shape, scale = scale, lower.tail = lower.tail,
                 log.p = log.p) + mu
   return(q)
}

#' @name rgamma3p
#' @rdname gamma3p
#' @title Gamma Distribution with Three Parameters
#' @importFrom stats qgamma
#' @export
rgamma3p <- function(n, shape = 1, rate = 1, scale = 1/rate, mu = 0) {
   r <- rgamma(n = n, shape = shape, scale = scale) + mu
   return(r)
}

