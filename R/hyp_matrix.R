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
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along
## with this program; if not, see <http://www.gnu.org/licenses/>.

#' @rdname hyp_matrix
#' @title Compute the Hypothesis Matrix for MANOVA testing.
#' @description Given the summarized data provided, this function compute
#' the hypothesis matrices used in several multivaraite statistics.
#' @details To perform the comparison of a model under different parameter
#' values (different curve fits with the same model), the covariance matrices
#' and means from each set of model parameters can be used in place of the
#' covariance matrices and means for the variables inside of each group.
#' @return A list carrying the following matrices and numerical data is
#' returned:
#'
#' \itemize{
#'   \item "Pooled" sums of squares and cross product (SSCP) matrix: 'E'
#'   \item Degree of freedom for matrix E: 'df.e'
#'   \item The hypothesis SSCP matrix: 'H'
#'   \item Degree of freedom for matrix H: 'df.h'
#'   \item 'H + E': 'HpE.inv'
#'   \item Error covariance matrix: 'Se'
#'   \item Inverse of Se: Se.inv
#'   \item Sum of mean of Squares: 'W'
#'   \item Inverse of W: 'W.inv'
#'   \item Variance-covariance matrix from each model: 'Cov'
#'   \item Parameter matrix for the groups: 'par'
#'   \item Number of variables or model parameters: 'p'
#'   \item Total sample size: 'N'
#'   \item Number of models/groups: 'n'
#'   \item Sample size: 'm'
#' }
#' @param obj A list of ANOVA R objects or linear/non-linear model objects.
#' If "obj" is not given, then the parmeters 'par.mat', 'cov.mat', and
#' 'sample.size' must be given.
#' @param par.mat A matrix where each column contains the model parameters of
#' each curve.
#' @param cov.mat List of parameter variance-covariance matrices.
#' @param sample.size Vector with the sample sizes used to estimate each
#' models.
#' @param groups list of sets of indices that identify models from a same
#' group.
#' @seealso \code{\link{MANOVA}}.
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @examples
#' ## Generate a dataset
#' set.seed(1230)
#' x1 <- rnorm(1e3, mean = 2.5, sd = 1)
#' x2 <- rnorm(1e3, mean = 2.5, sd = 1)
#'
#' line_1 <- 2.5 * x1 + 3 + runif(length(x1))/10
#' line_2 <- 2.5 * x2 + 3 + runif(length(x1))/10
#'
#' ## Fitting the linear models
#' lm_1 <- lm(line_1 ~ x1)
#' lm_2 <- lm(line_2 ~ x2)
#' res <- hyp_matrix(obj = list(model1 = lm_1, model2 = lm_2))

hyp_matrix <- function(
        obj,
        par.mat = NULL,
        cov.mat = NULL,
        sample.size = NULL,
        groups = NULL) {

    datos <- data_formating(obj = obj,
                            par.mat = par.mat,
                            cov.mat = cov.mat,
                            sample.size = sample.size,
                            groups = groups)
    m <- datos$m
    n <- datos$n
    N <- datos$N
    p <- datos$p
    df.h <- datos$df.h
    par <- datos$par
    Cov <- datos$Cov

    ## ========================================================== #
    ## "pooled" sums of squares and cross product matrix E (SSCP):
    ## ========================================================== #
    S <- mapply(function(x, y) (x - 1) * y, m, Cov, SIMPLIFY = F)
    E <- Reduce("+", S)
    df.e <- N - n # Degree of freedom for matrix E
    Se <- E / df.e # Error covariance matrix
    A <- mapply(function(x,y) x/y, S, m)

    ## ------------------------------------------------------ -#
    ## ----------- Estimation of hypothesis matrix -----------
    ## ------------------------------------------------------ -#
    ## vector of the grand means of all parameter in the models
    par.mean <- rowMeans(par)

    # SSCP for each group according to the hypothesis H:
    SSCPk <- list()
    for (k in seq_len(n)) {
        SSCPk[[k]] <- m[k] * (par[, k] - par.mean) %*% t(par[, k] - par.mean)
    }

    # The hypothesis SSCP matrix H:
    H <- Reduce("+", SSCPk)
    r <- min(p, df.h) # Necessary to compute dfs
    s <- max(p, df.h)
    #---------------------------------------- -

    W <- Reduce("+", A) # Sum of mean of Squares
    W.inv <- suppressWarnings(try(qr.solve(W), silent = TRUE))
    if (inherits(W.inv, "try-error")) {
        # error handling code, maybe just skip this iteration using
        Se.inv <- try(solve(qr(Se, LAPACK = TRUE)), silent = TRUE)
        HpE.inv <- try(solve(qr(H + E, LAPACK = TRUE)), silent = TRUE)
        E.inv <- try(solve(qr(E, LAPACK = TRUE)), silent = TRUE)
        W.inv <- try(solve(qr(W, LAPACK = TRUE)), silent = TRUE)

        matrix_inv <- c(
            "'Error-covariance-matrix'\n", "'SSCP'\n",
            "'SSCP + SSCP-hypothesis\n'",
            "'Sum-of-mean-of-Squares'\n"
        )

        try_error <- sapply(c(W.inv, Se.inv, HpE.inv, E.inv), inherits)
        if (any(try_error)) {
            stop(
                "*** It was possible to compute the inverse matrix for ",
                "some of the matrices:\n", matrix_inv[try_error],
                "\nCheck the covariance matrices of your models.\n",
                "'vcov(model)'"
            )
        }
    } else {
        Se.inv <- qr.solve(Se, tol = 1e-10)
        HpE.inv <- qr.solve(H + E, tol = 1e-10)
        E.inv <- qr.solve(E, tol = 1e-10)
    }

    invisible(
            list(
                E = E,
                df.e = df.e,
                E.inv = E.inv,
                H = H,
                A = A,
                df.h = df.h,
                HpE.inv = HpE.inv,
                Se = Se,
                Se.inv = Se.inv,
                W = W,
                W.inv = W.inv,
                par = par,
                Cov = Cov,
                p = p,
                N = N,
                df.h = df.h,
                m = m,
                n = n)
            )
}

## ================== Auxiliary function =========================

data_formating <- function(
                        obj,
                        par.mat = NULL,
                        cov.mat = NULL,
                        sample.size = NULL,
                        groups = NULL) {

    if (!missing(obj)) {
        n <- length(obj) ## number of models/groups
        Cov <- list()
        m <- c()
        par <- c()

        nls_lm_class <- all(sapply(obj, function(ob)
            inherits(ob, "nls.lm")))
        nls_class <- all(sapply(obj, function(ob)
            inherits(ob, "nls")))

        if (nls_class || nls_lm_class) {
            for (k in seq_len(n)) {
                ## To get the parameter variance-covariance matrix
                ## for each model
                Cov[[k]] <- unname(vcov(obj[[k]]))

                ## This generates a vector with the sample sizes used
                ## to estimate each models
                m <- c(m, sum(summary(obj[[k]])$df))

                ## This generated a matrix where each column contains
                ## the model parameters of each curve.
                par <- cbind(par, coef(obj[[k]]))
            }
        }

        lm_class <- all(sapply(obj, function(ob) inherits(ob, "lm")))
        aov_class <- all(sapply(obj, function(ob) inherits(ob, "aov")))

        if (lm_class || aov_class) {
            for (k in seq_len(n)) {
                ## To get the parameter variance-covariance matrix for each
                ##  model
                Cov[[k]] <- unname(vcov(obj[[k]]))

                ## It generates a vector with the sample sizes
                m <- c(m, length(obj[[k]]$fitted.values))

                ## It generates a matrix where each column contains group
                ## means or the model parameters of each group or linear
                ## model, respectively.
                par <- cbind(par, unname(coef(obj[[k]])))
            }
        }
    }

    if (is.null(par)) {
        stop(
            "*** Non-informative models provided. \n",
            "All models must belong to the same class.\n",
            "The supported classes are:\n",
            "'nls', 'nls.lm', 'lm', and 'aov'"
        )
    }

    if (missing(obj)) {
        Cov <- cov.mat
        par <- par.mat
        m <- sample.size
        n <- length(m) ## number of models/groups
    }

    if (!is.null(groups)) {
        M <- c()
        PAR <- c()
        COV <- list()
        for (k in seq_along(groups)) {
            ## Pooled covariance matrix for each group (joint covariance)
            COV[[k]] <- Reduce("+", mapply(function(x, y) (x - 1) * y,
                                           m[groups[[k]]],
                                           Cov[groups[[k]]],
                                           SIMPLIFY = F
            ))
            COV[[k]] <- COV[[k]] / sum(m[groups[[k]]] - 1)

            ## Sample size in group k
            M[k] <- sum(m[groups[[k]]])

            ## Parameter matrix for the groups
            PAR <- cbind(PAR, rowMeans(par[, groups[[k]]], na.rm = TRUE))
        }
        Cov <- COV
        par <- PAR
        m <- M
        n <- length(m)
    }

    p <- nrow(par) ## Number of variables or model parameters
    N <- sum(m) ## Total sample size
    df.h <- n - 1 ## Degree of freedom for matrix H
    return(list(par = par, Cov = Cov, p = p,
                N = N, df.h = df.h, m = m, n = n))
}

