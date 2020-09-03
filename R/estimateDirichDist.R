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

#' @rdname estimateDirichDist
#' @aliases estimateDirichDist
#' @title Nonlinear Parameter Estimation for Dirichlet Distribution \eqn{}
#' @description The parameter estimation is accomplished using a count data
#' matrix. The estimation is based on the fact that if a variable \eqn{x = (x_1,
#' x_2, ...x_n)} follows Dirichlet Distribution with parameters \eqn{\lapha =
#' \alpha_1, ... , \alpha_n} (all positive reals), in short, \eqn{x ~
#' Dir(\alpha)}, then \eqn{x_i ~ Beta(\alpha_i, \alpha_0 - \alpha_i)}, where
#' Beta(.) stands for the Beta distribution and \eqn{\alpha_0 = \sum \alpha_i}.
#'
#' Dirichlet distribution is a family of continuous multivariate probability
#' distributions, a multivariate generalization of the Beta distribution.
#'
#' @param x A matrix or a data.frame object carrying count data.
#' @param num.cores,tasks Parameters for parallel computation using
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param verbose if TRUE, prints the function log to stdout and a progress bar
#' @param ... Further arguments for \code{\link{betaDistEstimation}} function.
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @author Robersy Sanchez <https://genomaths.com>
#' @return A vector of estimated parameter values
#' @seealso \code{\link{betaDistEstimation}} and \code{\link{betaBinPost}}
#' @export
#' @examples
#' library(DirichletReg)
#'
#' ## A random generation numerical vectors with
#' x = rdirichlet(n = 1000, alpha = c(2.1, 3.1, 1.2))
#' head(x)
#'
#' estimateDirichDist(x)

estimateDirichDist <- function(  x,
                                num.cores = 1L,
                                tasks = 0L,
                                verbose = TRUE,
                                 ...) {

    if (!((is.vector(x) && is.numeric(x)) || is.matrix(x) || is.data.frame(x)))
        stop("\n*** 'x' must be a numerical vector, a matrix or a data.frame")

    if (!inherits(x, c("matrix", "data.frame")))
        stop("\n*** Object 'x' must be a 'matrix' or a 'data.frame'.")

    cn <- colnames(x)
    d <- dim(x)
    p <- x/rsum(x)

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

        pars <- bplapply(seq_len(d[2]), function(k) {

            parm <- try(betaDistEstimation(q = p[, k], ...)$parameters[1],
                        silent = TRUE)
            if (!inherits(parm, "try-error")) {
                parm <- try(betaDistEstimation( q = p[, k],
                                                force.optim = TRUE,
                                                ...)$parameters[1],
                            silent = TRUE)
            }

            if (inherits(parm, "try-error"))
                stop("\n*** Parameter for marginal '", k, "' failed")

            return(parm)
        }, ..., BPPARAM = bpparam)
        pars <- unlist(pars)
    } else {
        pars <- vector(mode = "numeric", length = d[2])
        for (k in seq_len(d[2])) {
            # if (verbose)
            parm <- try(betaDistEstimation(q = p[, k], ...)$parameters[1],
                        silent = TRUE)
            if (!inherits(parm, "try-error")) {
                parm <- try(betaDistEstimation( q = p[, k],
                                                force.optim = TRUE,
                                                ...)$parameters[1],
                            silent = TRUE)
            }

            if (inherits(parm, "try-error"))
                stop("\n*** Parameter for marginal '", k, "' failed")
            pars[k] <- parm
        }
    }
    if (is.null(cn)) names(pars) <- paste0("a", seq_len(d[2]))
    else names(pars) <- cn
    return(pars)
}

### ===================== Auxiliary function =================

rsum <- function(x) {
    if (length(dim(x)) > 1)  return(rowSums(x, na.rm = TRUE))
    else  return(sum(x, na.rm = TRUE))
}
