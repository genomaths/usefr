## Copyright (C) 2019 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez
##
## This program is part 'usefr' R package.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
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

#' @rdname betaBinPost
#' @aliases betaBinPost
#' @title Beta binomial posteriors
#' @description In a Bayesian framework, the parameters from Beta
#' distribution, \eqn{\alpha} and \eqn{\beta}, estimated from count data are
#' used as pseudo-counts to estimate the posterior probabilities:
#' \deqn{p = (\alpha + success)/(\alpha + \beta + trials)}.
#'
#' @param x A numerical vector, a matrix or a data.frame object.
#' @param num.cores,tasks Parameters for parallel computation using
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param verbose if TRUE, prints the function log to stdout and a progress bar
#' @param ... Further arguments for \code{\link{betaDistEstimation}} function.
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @author Robersy Sanchez <https://genomaths.com>
#' @seealso \code{\link{betaDistEstimation}} and
#' \code{\link{estimateDirichDist}}


betaBinPost <- function(x, num.cores = 1L, tasks = 0L, verbose = TRUE, ...) {

    if (!((is.vector(x) && is.numeric(x)) || is.matrix(x) || is.data.frame(x)))
        stop("\n*** 'x' must be a numerical vector, a matrix or a data.frame")

    d <- dim(x)
    if (is.null(d)) {
        l <- length(x)
        if (l < 2) return(NA)
    } else l = d[2]

    if (is.null(d)) {
        p <- (x + 1)/(rsum(x) + l)
        pars <- try(betaDistEstimation(p, ...)$parameters,
                    silent = TRUE)
        if (!inherits(pars, "try-error")) {
            p <- betaBinPosteriors(
                x,
                rsum(x),
                a = pars[1],
                b = pars[2]
            )
        }
    } else {
        n <- rsum(x)
        p <- (x + 1)/(n + l)

        if (num.cores > 1) {
            cn <- colnames(p)
            # Set parallel computation
            progressbar = FALSE
            if (verbose) progressbar = TRUE
            if (Sys.info()['sysname'] == "Linux") {
                bpparam <- MulticoreParam(workers = num.cores,
                                          tasks = tasks,
                                          progressbar = progressbar)
            } else bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                                        progressbar = progressbar)

            p <- bplapply(seq_len(l), function(k) {

                pars <- try(betaDistEstimation(p[, k], ...)$parameters,
                            silent = TRUE)
                if (!inherits(pars, "try-error")) {
                    post <- betaBinPosteriors(
                                            x[, k],
                                            n[k],
                                            a = pars[1],
                                            b = pars[2]
                            )
                } else post <- p[, k]
                return(post)
            }, ..., BPPARAM = bpparam)
            p <- do.call(cbind, p)
            colnames(p) <- cn
        } else {
            for (k in seq_len(l)) {
                pars <- try(betaDistEstimation(p[, k], ...)$parameters,
                            silent = TRUE)
                if (!inherits(pars, "try-error")) {
                    p[, k] <- betaBinPosteriors(
                                                x[, k],
                                                n[k],
                                                a = pars[1],
                                                b = pars[2]
                                            )
                }
            }
        }
    }
    return(p)
}

## ====================== Auxiliary function =============================

betaBinPosteriors <- function(success, trials, a, b) {
    (a + success)/(a + b + trials)
}
