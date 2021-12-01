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

#' @name dirichlet
#' @rdname dirichlet
#' @aliases dirichlet
#' @aliases rdirichlet
#' @aliases pdirichlet
#' @aliases ddirichlet
#'
#' @title Dirichlet Distribution
#' @description Probability density function (PDF), cumulative density function
#' (CDF), and random generation for the Dirichlet distribution
#' with parameters \eqn{\lapha = \alpha_1, ... , \alpha_n}, all positive reals,
#' in short, \eqn{Dir(\alpha)}. Dirichlet distribution is a family of continuous
#' multivariate probability distributions, a multivariate generalization of the
#' Beta distribution. The PDF is a multidimensional generalization of the beta
#' distribution.
#' @details The computation of the function 'pdirichlet' is accomplished
#' using Monte Carlo integration of the density 'dirichlet'. The estimation are
#' based on the fact that if a variable \eqn{x = (x_1, x_2, ...x_n)} follows
#' Dirichlet Distribution with parameters \eqn{\alpha = \alpha_1, ... ,
#' \alpha_n} (all positive reals), in short, \eqn{x ~ Dir(\alpha)}, then
#' \eqn{x_i ~ Beta(\alpha_i, \alpha_0 - \alpha_i)}, where Beta(.) stands for the
#' Beta distribution and \eqn{\alpha_0 = \sum \alpha_i}.
#'
#' The density is computed directly from the logarithm of Dirichlet PDF
#' \eqn{log(Dir(\alpha))}.
#'
#' if \eqn{x = (x_1, x_2, ...x_n)} are observed counts, then they are
#' transformed estimated probability \eqn{p_i} of the observation \eqn{x_i} is
#' estimated as: \eqn{p_i = (x_i + \alpha_i)/(\sum x_i + sum \alpa_i)}.
#'
#' @param x A numerical matrix \eqn{x_ij >= 0} or vector \eqn{x_i >= 0} of
#' observed counts  or relative frequencies \eqn{\sum x_i = 1}.
#' @param q numeric vector.
#' @param n number of observations.
#' @param p vector of probabilities.
#' @param alpha Dirichlet distribution's parameters. A numerical parameter
#' vector, strictly positive (default 1):
#' \eqn{\lapha = \alpha_1, ... , \alpha_n}.
#' @param lower.tail logical; if TRUE (default), probabilities are
#'     \eqn{P(X<=x)}, otherwise, \eqn{P(X > x)}
#' @param log.p logical; if TRUE, probabilities/densities p are returned as
#'     log(p).
#' @param log.base (Optional) The logarithm's base if log.p = TRUE.
#' @param interval The interval whether to search for the quantile (see
#' Details).
#' @param tol The desired accuracy (convergence tolerance) for the quantile
#' estimation.
#' @param ... Additional named or unnamed arguments to be passed to function
#' \code{\link[cubature]{hcubature}} used in functions 'pdirichlet' and
#' 'qdirichlet'.
#' @return GG PDF values (3-parameters or 4-parameters) for ddirichlet,
#' GG probability for pdirichlet, quantiles or GG random generated values for
#' rdirichlet.
#'
#' @aliases ddirichlet
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @references
#' 1. Sjolander K, Karplus K, Brown M, Hughey R, Krogh A, Saira Mian I, et al.
#'    Dirichlet mixtures: A method for improved detection of weak but
#'    significant protein sequence homology. Bioinformatics. 1996;12: 327–345.
#'    doi:10.1093/bioinformatics/12.4.327.
#' @examples
#' ## A random generation of numerical vectors
#' x = rdirichlet(n = 1e3, alpha = rbind(c(2.1, 3.1, 1.2),
#'                                       c(2.5, 1.9, 1.6)))
#' head(x[[1]])
#' lapply(x, estimateDirichDist)
#'
#' ### Density computation
#' alpha = cbind(1:10, 5, 10:1)
#' x = rdirichlet(10, c(5, 5, 10))
#' ddirichlet(x = x, alpha = a.mat)
#'
#' ### Or just one vector alpha
#' x = rdirichlet(n = 10, alpha = c(2.1, 3.1, 1.2))
#' ddirichlet(x = x, alpha = c(2.1, 3.1, 1.2))
#'
#' ## Compute Dirichlet probability
#' set.seed(123)
#' q = rdirichlet(10, alpha = c(2.1, 3.2, 3.3))
#' pdirichlet(q, alpha = c(2.1, 3.2, 3.3))
#'

ddirichlet <- function( x,
                        alpha,
                        log.p = FALSE) {

    ### The logarithm of the beta function expressed in terms of
    ### gamma functions

    d1 <- is.null(dim(x))
    d2 <- is.null(dim(alpha))
    if (d2)
        l <- length(alpha)
    else
        l <- dim(alpha)[2]

    if (inherits(x, c("matrix", "data.frame"))) {
        if (is.data.frame(x))
            x <- as.matrix(x)
    }
    else
        if (!is.numeric(x) || !length(x) > 1 || length(x) != l)
            stop("\n*** 'x' must be a 'matrix' or a 'data.frame' or",
                " a numerical vector length(x) == length(alpha)")

    if (inherits(alpha, c("matrix", "data.frame"))) {
        if (is.data.frame(alpha))
            alpha <- as.matrix(alpha)
    }
    else {
        if (length(alpha) < 2)
            stop("\n*** 'alpha' must be a matrix or a numerical vector",
                 " with length(alpha) = ncol(x).")
    }


    if (any.greater(x)) {
        p <- (x + alpha)/(rsum(x) + rsum(alpha))
    }


    logConst <- function(a) {
        if (is.null(dim(a)))
            sum(lgamma(a)) - lgamma(sum(a))
        else
            apply(a, 1, function(x) sum(lgamma(x)) - lgamma(sum(x)))
    }

    if (inherits(x, c("matrix", "data.frame"))) {
        if (d2)
            s <- (log(x) %*% (alpha - 1))
        else
            s <- rowSums((alpha - 1) * log(x))
    }
    else {
        if (d2)
            s <- sum((alpha - 1) * log(x), na.rm = TRUE)
        else
            s <- ((alpha - 1) %*% log(x))
    }

    if (log.p)
        pdf <- (s - logConst(alpha))
    else
        pdf <- (exp(s - logConst(alpha)))

    if (d1) {
        pdf[ !all.equal(sum(x), 1) ] <- 0
        pdf[ !all.equal(sum(abs(x)), 1) ] <- 0
        pdf <- as.vector(pdf)
    }
    if (d2 && !d1) {
        pdf[ !is.equal(rowSums(x), 1) ] <- 0
        pdf[ !is.equal(rowSums(abs(x)), 1) ] <- 0
        pdf <- as.vector(pdf)
    }
    if (!d2 && !d1) {
        pdf[ !is.equal(rowSums(x), 1) ] <- 0
        pdf[ !is.equal(rowSums(abs(x)), 1) ] <- 0
    }
    return(pdf)
}

#' @name pdirichlet
#' @rdname dirichlet
#' @title Dirichlet Distribution
#' @importFrom cubature hcubature
#'
#' @export
pdirichlet <- function(
                        q,
                        alpha,
                        lower.tail = TRUE,
                        log.p = FALSE,
                        log.base = exp(1), ...) {

    d1 <- is.null(dim(q))
    if (d1) nc <- length(q)
    else nc <- ncol(q)

    d2 <- is.null(dim(alpha))
    if (d2 && (length(alpha) == nc))
        l <- length(alpha)
    else
        stop("\n*** 'alpha' must be a numerical vector",
             " with length(alpha) = ncol(q).")

    if (inherits(q, c("matrix", "data.frame")) && nrow(q) > 1) {
        if (is.data.frame(q))
            q <- as.matrix(q)
    }
    else {
        if (!d1 && nrow(q) == 1) {
            q <- as.numeric(q)
            d1 <- TRUE
        }
        if (!is.numeric(q) || !length(q) > 1 || length(q) != l)
            stop("\n*** 'x' must be a 'matrix' or a 'data.frame' or",
                 " a numerical vector length(x) == length(alpha)")
    }

    if (d1) {
        intgrl <- integralDir(
                               lower = rep(0, (l - 1)),
                               upper = q[ seq_len(l - 1) ],
                               alpha = alpha)
    }
    else {
        lowerLimit <- matrix(0, nrow(q), (l - 1))
        upperLimit <- q[, seq_len(l - 1)]

        intgrl <- lapply(seq_len(nrow(lowerLimit)), function(k, ...)
            integralDir(lower = lowerLimit[k,],
                        upper = upperLimit[k,],
                        alpha = alpha, ...))
        intgrl <- unlist(intgrl)
    }

    if (!lower.tail)
        intgrl <- ( 1 - intgrl )
    if (log.p)
        intgrl <- log(intgrl, base = log.base)

    return(intgrl)
}

#' @name rdirichlet
#' @rdname dirichlet
#' @title Dirichlet Distribution
#' @importFrom stats rgamma
#' @export
rdirichlet <- function(n, alpha) {

    d <- dim(alpha)

    if (inherits(alpha, c("matrix", "data.frame"))) {
        rn <- rownames(alpha)
        cn <- colnames(alpha)
        if (is.matrix(alpha))
            alpha <- as.list(data.frame(t(alpha)))
        else
            alpha <- as.list(t(alpha))
        r <- lapply(alpha, function(a) {
                r <- sapply(a, function(a) rgamma(n, shape = a))
                d2 <- is.null(dim(r))
                if (d2)
                    r <- r/sum(r, na.rm = TRUE)
                else
                    r <- r/rowSums(r, na.rm = TRUE)
                if (is.null(cn) && d2)
                    names(r) <- paste0("a", seq_len(length(r)))
                else
                    if (is.null(cn))
                        colnames(r) <- paste0("a", seq_len(d[2]))
                return(r)
            })
        if (is.null(rn)) names(r) <- paste0("alpha", seq_len(d[1]))
        else names(r) <- rn
    }
    else {
        if (length(alpha) < 2)
            stop("\n*** 'alpha' must be a matrix or a numerical vector",
                " with length(alpha) = ncol(x).")
        nms <- names(alpha)
        r <- sapply(alpha, function(a) rgamma(n, shape = a))
        d2 <- is.null(dim(r))
        if (d2)
            r <- r/sum(r, na.rm = TRUE)
        else
            r <- r/rowSums(r, na.rm = TRUE)

        if (is.null(nms) && d2)
            names(r) <- paste0("a", seq_len(length(r)))
        else
            colnames(r) <- paste0("a", seq_along(alpha))
    }
    return(r)
}

### ==================== Auxiliary function ========================

is.equal <- Vectorize(all.equal, vectorize.args = "target")

ddirich <- function(x, alpha) {
    x <- c(x, (1 - sum(x)))
    x <- ddirichlet(x, alpha = alpha)
    return(x)
}

integralDir <- function(lower, upper, alpha, ...) {
    hcubature(ddirich, lowerLimit = lower, upperLimit = upper,
              alpha = alpha, ...)$integral
}


mc_pdir <- function(q, alpha, n = 1e4, seed = 1) {
                set.seed(seed)
                l <- seq_len(length(alpha) - 1)
                x <- rdirichlet(n, alpha = alpha)
                sum(apply(x[,l], 1, function(x) prod(x <= q[l]))) / n
}

any.greater <- function(x, value, d = 1) {
    if (d > 1) {
        res <- any(colSums(x) > value)
    } else
        res <- any(rowSums(x) > value)
    return(res)
}
