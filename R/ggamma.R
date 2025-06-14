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

#' @name ggamma
#' @aliases ggamma
#' @aliases rggamma
#' @aliases pggamma
#' @aliases dggamma
#' @aliases qggamma
#'
#' @title Generalized Gamma distribution
#' @description Probability density function (PDF), cumulative density function
#'     (CDF), quantile function and random generation for the Generalized Gamma
#'     (GG) distribution with 3 or 4 parameters: alpha, scale, mu, and psi. The
#'     function is reduced to GGamma distribution with 3 parameters by setting
#'     mu = 0.
#' @details Details about these function can be found in references 1 to 3. You
#' may also see section Note at ?pgamma or ?rgamma. Herein, we are using Stacy'
#' s formula (references 2 to 3) with the parametrization given in reference 4
#' (equation 6, page 12). As in the case of gamma distribution function, the
#' cumulative distribution function (as given in equation 12, page 13 from
#' reference 4) is expressed in terms of the lower incomplete gamma function
#' (see ?pgamma).
#'
#' The GG distribution with parameters \eqn{\alpha}, \eqn{\beta} (scale),
#' \eqn{\psi}, and \eqn{\mu} has density:
#'
#' \deqn{f(x | \alpha, \beta, \mu, \psi) = \alpha exp(-((x-\mu)/
#'     \beta)^\alpha) ((x-\mu)/\beta)^(\alpha * \psi - 1)/(\beta \Gamma(\psi))}
#'
#' The entropy, as given in reference 4, is computed using the closed form:
#' \deqn{S = {\ln \frac{\beta \,\Gamma \left( \delta \right)}
#' {\alpha} + \psi \left( \delta \right)\left( {\frac{1}{\alpha }
#' - \delta } \right) + \delta }},
#'
#' @param q numeric vector.
#' @param n number of observations.
#' @param p vector of probabilities.
#' @param alpha numerical parameter, strictly positive (default 1). The
#'     generalized gamma becomes the gamma distribution for alpha = 1.
#' @param scale,psi the same real positive parameters as is used for the Gamma
#'     distribution. These are numerical and strictly positives; default 1.
#'     (see ?pgamma).
#' @param mu location parameter (numerical, default 0).
#' @param lower.tail logical; if TRUE (default), probabilities are
#'     \eqn{P(X<=x)}, otherwise, \eqn{P(X > x)}
#' @param log.p logical; if TRUE, probabilities/densities p are returned as
#'     log(p).
#'
#' @return GG PDF values (3-parameters or 4-parameters) for dggamma,
#'     GG probability for pggamma, quantiles or GG random generated values for
#'     rggamma.
#'
#' @references 1. Handbook on  STATISTICAL DISTRIBUTIONS for experimentalists
#' (p. 73) by Christian Walck. Particle Physics Group Fysikum. University of
#' Stockholm (e-mail: walck@physto.se )
#'
#' 2. Stacy, E. W. A Generalization of the Gamma Distribution. Ann. Math. Stat.
#' 33, 1187–1192 (1962).
#'
#' 3. Stacy E, Mihram G (1965) Parameter estimation for a generalized gamma
#' distribution. Technometrics 7: 349-358.
#'
#' 4. Sanchez, R., & Mackenzie, S. A. (2023). On the thermodynamics of DNA
#' methylation process. Scientific Reports, 13(1), 8914.
#' https://doi.org/10.1038/s41598-023-35166-9.
#'
#' @examples
#' q <- (1:9) / 10
#' pggamma(q,
#'     alpha = 1, scale = 1, mu = 0,
#'     psi = 1, lower.tail = TRUE, log.p = FALSE
#' )
#'
#' ## To fit random generated numbers
#' set.seed(123)
#' x <- rggamma(2000, alpha = 1.03, psi = 0.75, scale = 2.1)
#' fitCDF(x, distNames = "Generalized 3P Gamma")
#'
#' @importFrom stats pgamma rgamma
#'
#' @aliases dggamma
#' @rdname ggamma
#' @title Generalized Gamma distribution
#' @export
dggamma <- function(q, alpha = 1, scale = 1, mu = 0,
    psi = 1, log.p = FALSE) {
    if (scale > 0 && alpha > 0 && psi > 0) {
        y <- (q - mu) / scale
        d <- exp(-y^alpha) * alpha * y^(alpha * psi -
            1) / (scale * gamma(psi))
        if (log.p) {
            return(log(d))
        } else {
            return(d)
        }
    } else {
        d <- NaN
    }
    return(d)
}
#'
#' @name pggamma
#' @rdname ggamma
#' @title Generalized Gamma distribution
#' @importFrom stats pgamma
#' @export
pggamma <- function(q, alpha = 1, scale = 1, mu = 0,
    psi = 1, lower.tail = TRUE, log.p = FALSE) {
    if (scale > 0 && alpha > 0 && psi > 0) {
        y <- ((q - mu) / scale)^alpha
        p <- pgamma(y,
            shape = psi, lower.tail = lower.tail,
            log.p = log.p
        )
    } else {
        p <- NaN
    }
    return(p)
}

#' @name qggamma
#' @rdname ggamma
#' @title Generalized Gamma distribution
#' @importFrom stats qgamma
#' @export
qggamma <- function(p, alpha = 1, scale = 1, mu = 0,
    psi = 1, lower.tail = TRUE, log.p = FALSE) {
    if (scale > 0 && alpha > 0 && psi > 0) {
        q <- qgamma(p,
            shape = psi, lower.tail = lower.tail,
            log.p = log.p
        )
        q <- scale * q^(1 / alpha) + mu
    } else {
        q <- NaN
    }
    return(q)
}

#' @name rggamma
#' @rdname ggamma
#' @title Generalized Gamma distribution
#' @importFrom stats rgamma
#' @export
rggamma <- function(n, alpha = 1, scale = 1, mu = 0, psi = 1) {
    if (scale <= 0) {
        stop("'scale' parameter must be positive")
    }
    if (alpha <= 0) {
        stop("'alpha' parameter must be positive")
    }
    if (psi <= 0) {
        stop("'psi' parameter must be positive")
    }

    r <- scale * (rgamma(n, psi * alpha)^(1 / alpha)) + mu
    return(r)
}


#' @name eggamma
#' @rdname ggamma
#' @title Generalized Gamma distribution
#' @export
eggamma <- function(alpha, scale, psi) {
    # Input validation
    if (alpha <= 0 || scale <= 0 || psi <= 0) {
        stop("Parameters alpha, scale, and psi must be positive")
    }

    # Compute entropy
    entropy <- -log(alpha) + log(scale) + lgamma(psi) +
        psi - ((alpha * psi - 1)/alpha) * digamma(psi)

    return(unname(entropy))
}
