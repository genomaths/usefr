## ######################################################################## #
## 
## Copyright (C) 2019 - 2024 Robersy Sanchez <https://genomaths.com/>
##     Author: Robersy Sanchez Rodriguez
## 
## This file is part of the R package "fitCDF".
## 
## Permission is hereby granted to a person or institution obtaining a copy of 
## this software (fitCDF) and associated documentation files (the "Software"), 
## to use it under the following conditions:
##     
## i) Any use of this software, commercial or non-commercial purposes including 
## teaching, academic and government research, public demonstrations, personal
## experimentation, and/or the evaluation of the Software requires a license
## from the author.
## 
## ii) Commercial use of the Software requires a commercial license. The author
## will negotiate commercial licenses upon request. These requests can
## be directed to https://genomaths.com/.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHOR, PERSON, INSTITUTION, ENTITY OR AFFILIATES, BE LIABLE FOR ANY CLAIM,
## DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
## OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
## USE OR OTHER DEALINGS IN THE SOFTWARE.
## 
## ######################################################################## #
## 

#' Beta-Geometric
#' 
#' @rdname bgeom
#' @name bgeom
#' @aliases bgeom
#' @aliases dbgeom
#' @aliases pbgeom
#' @aliases qbgeom
#' @aliases rbgeom
#' @aliases bgeom2
#' @aliases dbgeom2
#' @aliases pbgeom2
#' @aliases qbgeom2
#' @aliases rbgeom2
#'
#' @description 
#' Probability density function (PDF), cumulative density function (CDF),
#' quantile function, and random variate generation of the Beta-Geometric 
#' distribution (BGD). Herein, we consider the \strong{unshifted} sBGD PDF 
#' (dbgeom):
#' 
#' \deqn{f(x | \alpha,\beta) = 
#' \frac{B(\alpha + 1, \beta + x - 1)}{B(\alpha,\beta)}
#' }
#'
#' for \eqn{x = 1,2,...}, and \eqn{0 < p \leq 1} and the \strong{shifted} 
#' PDF (dbgeom2):
#' 
#' \deqn{f(x | \alpha,\beta) = 
#' \frac{B(\alpha + 1, \beta + x)}{B(\alpha,\beta)}
#' }

#' with parameters \eqn{\alpha} (shape1) and \eqn{\beta} (shape1) as given in
#' references (1-2), for \eqn{x = 0,1,2,...}, and \eqn{0 < p \leq 1}, where 
#' \eqn{B(.)} is the \code{\link[base]{beta}}(a,b) function:
#' 
#' \deqn{B(a,b)= \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}}
#' 
#' @details 
#' Details about the BGD can be found in reference (1-3). 
#' 
#' ## The Cumulative Distribution Function (CDF)
#' The CDF is computed as:
#' 
#' \deqn{F(x|\alpha, \beta) = \Sigma_{i=0}^{x} f(x |\alpha, \beta)}
#' 
#' ## The quantile function
#' 
#' The quantile function is the inverse of the CDF, i.e., \eqn{F^{-1}}. Since
#' there is not analytical function for the CDF, we compute the quantile
#' \eqn{\hat{q}_i} for the probability value \eqn{p_i} as the argument
#' \eqn{x_i} that minimizes the difference 
#' \eqn{1/2 (F(x_i|\alpha, \beta) - p_i)}, formally:
#' 
#' \deqn{\hat{q}_i = \underset{x \in \mathbb{Z}^+}{\operatorname{Min}}
#' \left \{ \frac{1}{2} (F(x_i|\alpha, \beta) - p_i)^2 \:
#' | \: F(x_i|\alpha, \beta) \leq p_i \right \}}
#' 
#' ## Random variate generation (random sampling)
#' 
#' To generate a random variate 
#' \eqn{X \sim \: F(x|\alpha, \beta)}, first, we
#' must generate a random variate \eqn{\pi} (probabilities): 
#' \eqn{\pi \sim \: B_x(x|\alpha, \beta)} (\code{\link[stats]{rbeta}}). Hence, 
#' we can generate \eqn{X \sim \: F(x|\alpha, \beta)} as:
#' 
#' \deqn{X = G^{-1}(R_B(n|\alpha, \beta) | \pi)}
#' 
#' where \eqn{R_B(n|\alpha, \beta)} stands for the function to generate random
#' variate from Beta distribution: \code{\link[stats]{rbeta}} and
#' \eqn{G^{-1}(p|\pi)} the quantile of the geometric distribution.
#' 
#' @param x,q A numeric vector.
#' @param p  A vector of probabilities.
#' @param n The number of observations
#' @param shape1,shape2 shape parameters to pass to Beta distribution,
#' \code{\link[stats]{Beta}}.
#' @param qmax The maximum expected quantile. An argument to be used with
#' \emph{qbgeom} function.
#' 
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P(X<=x)}, otherwise, \eqn{P(X > x)}.
#' @param log.p,log logical; if TRUE, probabilities/densities p are returned
#' as log(p).
#' @return The PDF, CDF, quantile and random generation functions for BGD.
#' 
#' @references 
#' 1. Weinberg, P. & Gladen, B.C. (1986). The Beta-Geometric distribution 
#'    applied to comparative fecundability studies. Biometrics, 42, 547-560.
#' 2. Gupta, A. K., Nadarajah, S (2004). The Beta-Geometric Distribution. In 
#'    Handbook of Beta Distribution and Its Applications Statistics. 
#' 3. Paul, Sudhir R. (2005) "Testing Goodness Of Fit Of The Geometric 
#'    Distribution: An Application To Human Fecundability Data,"
#'    Journal of Modern Applied Statistical Methods: Vol. 4 : Iss. 2 , 
#'    Article 8. DOI: 10.22237/jmasm/1130803620
#'    
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @import stats
#' @examples
#' ## Generate random variates from the unshifted Beta-Geometric
#' set.seed(122)
#' x1 <- rbgeom(1e3, shape1 = 2, shape2 = 40)
#' summary(x1)
#' 
#' ## The nonlinear fit 
#' cdf <- fit_cdf(x1, distNames = "Beta-Geometric I",
#'               plot = T)
#' cfs <- coefs(cdf)
#' cdf
#' 
#' ## The histogram with the density (PDF) curve.
#' hist(x1, 600, freq = FALSE, xlim = c(0, 50),
#'      panel.first={points(0, 0, pch=16, cex=1e6, col="grey95")
#'          grid(col="white", lty = 1)})
#' x <- round(seq(min(x1), max(x1), by = 1))
#' lines(x, dbgeom(x, shape1 = 2, shape2 = 40), col = "red")
#' lines(x, dbgeom(x, shape1 = cfs[1], shape2 = cfs[2]), col = "dodgerblue")
#' 
#' ## Kolmogorov-Smirnov tests can fail
#' set.seed(122)
#' x2 <- rbgeom(1e3, shape1 = 2, shape2 = 30)
#' summary(x2)
#' 
#' cdf <- fit_cdf(x2, distNames = "Beta-Geometric I",
#'               plot = TRUE)
#' cdf
#' 
#' ## Bootstrap test for Goodness of fit (GoF)
#' ## using Kolmogorov-Smirnov (for the sake of reducing compuational time,
#' ## the sample.size was set small, default is NULL).
#' mcgoftest(x2, cdf, stat = "ks", sample.size = 100)
#' 
#' 
#' ## Generate random variates from the shifted Beta-Geometric
#' set.seed(122)
#' x1 <- rbgeom2(1e3, shape1 = 2, shape2 = 40)
#' summary(x1)
#' 
#' ## The nonlinear fit 
#' cdf_sf <- fit_cdf(x1, distNames = "Beta-Geometric II",
#'             plot = T)
#' cdf_sf
#' 
#' @aliases dbgeom
#' @export
dbgeom <- function(
        x,
        shape1,
        shape2,
        log = FALSE) {
    
    if (shape1 > 0 && shape2 > 0) {
        
        if (log)
            x <- exp(x)
        
        x <- abs(x)
        x <- floor(x)
        
        x <- beta(1 + shape1, shape2 + x - 1) / beta(shape1, shape2)
        
        if (log)
            x <- log(x)
    }
    else 
        x <- NaN
    
    return(x)
}


#' @rdname bgeom
#' @aliases pbgeom
#' @title Beta-Geometric
#' @export
pbgeom <- function(
        q, 
        shape1,
        shape2,
        lower.tail = TRUE, 
        log.p = FALSE) {
    
    if (shape1 > 0 && shape2 > 0) {
        q <- abs(q)
        q <- floor(q)
        p <- 0 * q
        p[q == 0] <- dbgeom(1, 
                            shape1 = shape1, 
                            shape2 = shape2)
        unq <- sort(unique(q))
        idx <- which(unq > 0)
        if (length(idx) > 1) {
            unq <- unq[idx]
            sq <- seq(1, max(unq, na.rm = TRUE))
            ps <- cumsum(dbgeom(sq, 
                                shape1 = shape1, 
                                shape2 = shape2))[-1]
            sq <- sq[-1]
            
            for (k in seq_along(sq)) {
                idx <- which(q == sq[k])
                p[idx] <- ps[k]
            }
        }
        else {
            if (length(idx) == 1) {
                ps <- sum(dbgeom(seq(1, unq), 
                                shape1 = shape1, 
                                shape2 = shape2))
                idx <- which(q == unq)
                p[idx] <- ps
            }
        }
    }
    else 
        p <- NaN
    
    return(p)
}

#' @rdname bgeom
#' @aliases qbgeom
#' @title Beta-Geometric
#' @export
qbgeom <- function(
        p, 
        shape1,
        shape2,
        qmax = 1e4,
        lower.tail = TRUE, 
        log.p = FALSE) {
    
    if (any(p > 1 | p < 0))
        p[which(p > 1 | p < 0)] <- 0
    
    if (shape1 > 0 && shape2 > 0) {
        qs <- seq(1, qmax)
        ps <- cumsum(dbgeom(qs,
                             shape1 = shape1,
                             shape2 = shape2))
        
        if (length(p) > 1) {
            for (k in seq_along(p)) {
                idx <- which.min(1 / 2 * abs(ps - p[k])^2)
                p[k] <- qs[idx]
            }
        }
        else {
            idx <- which.min(1 / 2 * abs(ps - p)^2)
            p <- qs[idx]
        }
    }
    else 
        p <- NaN
    
    return(p)
}

#' @rdname bgeom
#' @aliases rbgeom
#' @title Beta-Geometric
#' @export
rbgeom <- function(
        n, 
        shape1,
        shape2) {
    
    if (shape1 > 0 && shape2 > 0) {
        
        n <- rgeom(
            n,
            prob = rbeta(
                n,
                shape1 = shape1,
                shape2 = shape2)
        )
        n <- floor(n)
        n <- n + 1
    }
    else 
        n <- NaN
    
    return(n)
}

 
 
#' @rdname bgeom
#' @aliases dbgeom2
#' @title Beta-Geometric
#' @export
dbgeom2 <- function(
        x,
        shape1,
        shape2,
        log = FALSE) {
    
    if (shape1 > 0 && shape2 > 0) {
        
        if (log)
            x <- exp(x)
        
        x <- abs(x)
        x <- floor(x)
        
        x <- beta(1 + shape1, shape2 + x) / beta(shape1, shape2)

        if (log)
            x <- log(x)
    }
    else 
        x <- NaN
    
    return(x)
}


#' @rdname bgeom
#' @aliases pbgeom2
#' @title Beta-Geometric
#' @export
pbgeom2 <- function(
        q, 
        shape1,
        shape2,
        lower.tail = TRUE, 
        log.p = FALSE) {
    
    if (shape1 > 0 && shape2 > 0) {
        q <- abs(q)
        q <- floor(q)
        p <- 0 * q
        p[q == 0] <- dbgeom2(0, 
                            shape1 = shape1, 
                            shape2 = shape2)
        unq <- sort(unique(q))
        idx <- which(unq > 0)
        if (length(idx) > 1) {
            unq <- unq[idx]
            sq <- seq(0, max(unq, na.rm = TRUE))
            ps <- cumsum(dbgeom2(sq, 
                                shape1 = shape1, 
                                shape2 = shape2))[-1]
            sq <- sq[-1]
            
            for (k in seq_along(sq)) {
                idx <- which(q == sq[k])
                p[idx] <- ps[k]
            }
        }
        else {
            if (length(idx) == 1) {
                ps <- sum(dbgeom2(seq(0, unq), 
                                shape1 = shape1, 
                                shape2 = shape2))
                idx <- which(q == unq)
                p[idx] <- ps
            }
        }
    }
    else 
        p <- NaN
    
    return(p)
}

#' @rdname bgeom
#' @aliases qbgeom2
#' @title Beta-Geometric
#' @export
qbgeom2 <- function(
        p, 
        shape1,
        shape2,
        qmax = 1e4,
        lower.tail = TRUE, 
        log.p = FALSE) {
    
    if (any(p > 1 | p < 0))
        p[which(p > 1 | p < 0)] <- 0
    
    if (shape1 > 0 && shape2 > 0) {
        qs <- seq(0, qmax)
        ps <- cumsum(dbgeom2(qs,
                            shape1 = shape1,
                            shape2 = shape2))
        
        if (length(p) > 1) {
            for (k in seq_along(p)) {
                idx <- which.min(1 / 2 * abs(ps - p[k])^2)
                p[k] <- qs[idx]
            }
        }
        else {
            idx <- which.min(1 / 2 * abs(ps - p)^2)
            p <- qs[idx]
        }
    }
    else 
        p <- NaN
    
    return(p)
}

#' @rdname bgeom
#' @aliases rbgeom2
#' @title Beta-Geometric
#' @export
rbgeom2 <- function(
        n, 
        shape1,
        shape2) {
    
    if (shape1 > 0 && shape2 > 0) {
        
        n <- rgeom(
                n,
                prob = rbeta(
                    n,
                    shape1 = shape1,
                    shape2 = shape2))
        
        n <- floor(n)
     }
    else 
        n <- NaN
    
    return(n)
}

