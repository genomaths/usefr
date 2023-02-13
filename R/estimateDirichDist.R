## Copyright (C) 2019-23 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez
##
## This program is part 'usefr' R package.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the Free
## Software Foundation; either version 3 of the License, or (at your option)
## any later version.
##
## 'usefr' is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

#' @rdname estimateDirichDist
#' @aliases estimateDirichDist
#' @title Nonlinear Parameter Estimation for Dirichlet Distribution \eqn{}
#' @description The parameter estimation is accomplished using a count data
#' matrix. The estimation is based on the fact that if a variable
#' \eqn{x = (x_1, x_2, ...x_n)} follows Dirichlet Distribution with parameters
#' \eqn{\alpha = \alpha_1, ... , \alpha_n} (all positive reals), in short,
#' \eqn{x ~ Dir(\alpha)}, then \eqn{x_i ~ Beta(\alpha_i, \alpha_0 - \alpha_i)},
#' where Beta(.) stands for the Beta distribution and
#' \eqn{\alpha_0 = \sum \alpha_i}.
#'
#' Dirichlet distribution is a family of continuous multivariate probability
#' distributions, a multivariate generalization of the Beta distribution.
#'
#' @details As any non-linear fitting, results strongly depends on the start
#' parameter values.
#' @param x A matrix or a data.frame object carrying count data.
#' @param start Initial parameter values for \eqn{\alpha =
#' \alpha_1, ... , \alpha_n} (all positive reals). Defaults is NULL.
#' @param num.cores,tasks Parameters for parallel computation using
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#' use, i.e. at most how many child processes will be run simultaneously (see
#' \code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
#' for Linux OS).
#' @param verbose if TRUE, prints the function log to stdout and a progress bar
#' @param ... Further arguments for \code{\link{betaDistEstimation}} function.
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @author Robersy Sanchez <https://genomaths.com>
#' @return A vector of estimated parameter values
#' @seealso \code{\link{betaDistEstimation}} and \code{\link{betaBinPost}}
#' @export
#' @examples
#' #' ## A random generation numerical vectors with
#' x <- rdirichlet(n = 1000, alpha = c(2.1, 3.1, 1.2))
#' head(x)
#'
#' estimateDirichDist(x)
#'
estimateDirichDist <- function(x,
    start = NULL,
    num.cores = 1L,
    tasks = 0L,
    seed = 123,
    refit = TRUE,
    verbose = TRUE,
    ...) {
    if (!((is.vector(x) && is.numeric(x)) || is.matrix(x) ||
        is.data.frame(x))) {
        stop("\n*** 'x' must be a numerical vector, a matrix or a data.frame")
    }

    if (!inherits(x, c("matrix", "data.frame"))) {
        stop("\n*** Object 'x' must be a 'matrix' or a 'data.frame'.")
    }

    cn <- colnames(x)
    d <- dim(x)
    p <- x / rsum(x)

    if (is.null(start)) {
        start <- matrix(1, d[1], 2)
    } else {
        if (length(start) != ncol(x)) {
            warning(
                "*** Wrong 'start' parameter length.",
                " The length(start) == ncol(x) \n",
                "The 'start' values will be ignored"
            )
        }
        alfa <- sum(start)
        start <- cbind(
            shape1 = start,
            shape2 = sapply(start, function(x) alfa - x)
        )
    }

    FIT <- beta_fitting(
        p = p,
        start = start,
        num.cores = num.cores,
        tasks = tasks,
        seed = seed,
        verbose = verbose,
        ...
    )

    if (refit) {
        start <- coefs(FIT)
        alfa_0 <- sum(start)
        start <- cbind(
            shape1 = start,
            shape2 = sapply(start, function(x) alfa_0 - x)
        )

        FIT <- beta_fitting(
            p = p,
            start = start,
            num.cores = num.cores,
            tasks = tasks,
            seed = seed,
            verbose = verbose,
            ...
        )
    }

    return(FIT)
}


### ===================== Auxiliary function =================

rsum <- function(x) {
    if (length(dim(x)) > 1) {
        return(rowSums(x, na.rm = TRUE))
    } else {
        return(sum(x, na.rm = TRUE))
    }
}


## ================ Auxiliary function ======================

beta_fitting <- function(p,
    start = NULL,
    num.cores = 1L,
    tasks = 0L,
    seed = 123,
    verbose = TRUE,
    ...) {
    cn <- colnames(p)
    d <- dim(p)

    if (num.cores > 1) {
        cn <- colnames(p)
        # Set parallel computation
        progressbar <- FALSE
        if (verbose) progressbar <- TRUE
        if (Sys.info()["sysname"] == "Linux") {
            bpparam <- MulticoreParam(
                workers = num.cores,
                tasks = tasks,
                progressbar = progressbar
            )
        } else {
            bpparam <- SnowParam(
                workers = num.cores, type = "SOCK",
                progressbar = progressbar
            )
        }

        FIT <- bplapply(seq_len(d[2]), function(k) {
            fit <- try(betaDistEstimation(
                q = p[, k],
                init.pars = start[k, ],
                seed = seed,
                ...
            ),
            silent = TRUE
            )
            if (inherits(fit, "try-error")) {
                fit <- try(betaDistEstimation(
                    q = p[, k],
                    init.pars = start[k, ],
                    force.optim = TRUE,
                    seed = seed,
                    ...
                ),
                silent = TRUE
                )
            }

            if (inherits(fit, "try-error")) {
                stop("\n*** Model for marginal '", k, "' failed")
            }

            fit <- structure(fit, class = c("BetaModel", "data.frame"))

            return(fit)
        }, ..., BPPARAM = bpparam)
        ALFA <- sapply(FIT, function(x) x$Estimate[1])
    } else {
        pars <- vector(mode = "numeric", length = d[2])
        ALFA <- vector(mode = "numeric", length = d[2])
        FIT <- vector(mode = "list", length = d[2])

        for (k in seq_len(d[2])) {
            # if (verbose)
            fit <- try(betaDistEstimation(
                q = p[, k],
                init.pars = start[k, ],
                seed = seed
            ),
            silent = TRUE
            )
            if (inherits(fit, "try-error")) {
                fit <- try(betaDistEstimation(
                    q = p[, k],
                    init.pars = start[k, ],
                    force.optim = TRUE,
                    seed = seed,
                    ...
                ),
                silent = TRUE
                )
            }

            if (inherits(fit, "try-error")) {
                stop("\n*** Parameter for marginal '", k, "' failed")
            }

            fit <- structure(fit, class = c("BetaModel", "data.frame"))

            ALFA[k] <- fit$Estimate[1]
            FIT[[k]] <- fit
        }
    }
    if (is.null(cn)) {
        names(FIT) <- paste0("beta_", seq_len(d[2]))
    } else {
        names(FIT) <- cn
    }

    FIT <- list(
        alpha = ALFA,
        marginals = FIT
    )
    FIT <- structure(FIT, class = c("DirchModel", "list"))
    return(FIT)
}
