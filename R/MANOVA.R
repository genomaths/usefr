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

##*===================================================================== #
##
## -------- Multivariate Analysis of Variance (MANOVA) -------
## -------------------- FOR SUMMARIZED DATA -----------------#
##
##====================================================================== #
#' @rdname MANOVA
#' @title Multivariate Analysis of Variance for Summarized Data
#' @description This is a MANOVA rest for summarized data. This means that if
#' we want to perform the multivariate comparison of, for example, three
#' groups, then we can use the covariance matrices and means for the variables
#' inside of each group to compute the statistics for a MANOVA.
#' @details To perform the comparison of a model under different parameter
#' values (different curve fits with the same model), the covariance matrices
#' and means from each set of model parameters can be used in place of the
#' covariance matrices and means for the variables inside of each group.
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
#' @param boxM Logic. Whether to perform the Box M test for covaraince matrix
#' equality. There are two limitations associated with the Box test ([1] pag
#' 42). First, it can be sensitive to multivariate nonnormality. Second, the
#' degrees of freedom for this test can often be very large, resulting in an
#' extremely sensitive/powerful test of the null hypothesis. Both of these
#' limitations could lead a researcher to question the validity of the
#' multivariate test on the vectors of means when the multivariate test is
#' valid. As a result, researchers often do not rely heavily on the results
#' of the test but rather rely on the robustness of the multivariate test on
#' the equality of mean vectors when sample sizes are equal.
#' @return If boxM = TRUE and the p-value of the test is greater than 0.05,
#' the function will return a list with data frame as the first element
#' carrying the statistics 'Hotelling T^2', 'Pillai's trace', 'Wilks lambda',
#' and 'Hotelling-Lawley' trace' and the corresponding p-values. Second
#' element will be the 'Box M' test result. If Box M' test result is not
#' significant then the data frame with the statistics will include, in
#' addition "Yao's F", "Johansen's F", and the "Nel-Merwe's F".
#'
#' If boxM = FALSE, then the function will return a data frame with the
#' statistics "Hotelling T^2", "Wilks lambda" and their corresponding
#' p-values.
#' @importFrom stats pchisq pf vcov coef
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @examples
#' ## =========== Example 1 ============
#' library(minpack.lm)
#' ## Build an initial random data set
#' set.seed(1230)
#' x1 <- rnorm(10000, mean = 0.5, sd = 1)
#' x2 <- rnorm(10000, mean = 0.6, sd = 1.2)
#'
#' cdfp1 <- fitCDF(x1, distNames = 1, plot = FALSE)
#' cdfp2 <- fitCDF(x2, distNames = 1, plot = FALSE)
#'
#' MANOVA(list(model1 = cdfp1$fit[[1]], model2 = cdfp2$fit[[1]]))
#'
#' ## =========== Example 2 ============
#' ## Define a non-linear function
#' mutr <- function(x, alpha, beta, lambda) {
#'     alpha / (1 + beta * exp(-lambda * x))
#' }
#'
#' ## Build two data set of points from specific curves
#' x1 <- sort(x1[x1 > 0]) * 10
#' y1 <- mutr(x1, alpha = 25, beta = 4.5, lambda = 2.1) + runif(length(x1))
#' y2 <- mutr(x1, alpha = 25.5, beta = 4.5,
#'             lambda = 2.05) + runif(length(x1))
#'
#' ## Non-linear fit
#' fit1 <- nlsLM(Y ~ alpha / (1 + beta * exp(-lambda * X)),
#'     data = data.frame(X = x1, Y = y1),
#'     start = list(alpha = 25.5, beta = 4.5, lambda = 2.05),
#'     control = list(
#'         maxiter = 1024, tol = 1e-12,
#'         minFactor = 10^-6
#'     )
#' )
#'
#' fit2 <- nlsLM(Y ~ alpha / (1 + beta * exp(-lambda * X)),
#'     data = data.frame(X = x1, Y = y2),
#'     start = list(alpha = 25, beta = 4.5, lambda = 2.1),
#'     control = list(
#'         maxiter = 1024, tol = 1e-12,
#'         minFactor = 10^-6
#'     )
#' )
#'
#' ## Manova result
#' MANOVA(list(model1 = fit1, model2 = fit2))
#'
#' @references
#' \enumerate{
#'     \item CARL J. HUBERTY, STEPHEN OLEJNIK.2006. Applied MANOVA and
#'         Discriminant Analysis. Second Edition. John Wiley & Sons, Inc.,
#'         Hoboken, New Jersey.
#'     \item Muller KE, Lavange LM, Ramey SL, Ramey CT. Power Calculations
#'         for General Linear Multivariate Models Including Repeated Measures
#'         Applications. J Am Stat Assoc. 1992;87: 1209 - 1226.
#'         doi:10.1080/01621459.1992.10476281.
#'     \item Finch H, French B (2013) A Monte Carlo Comparison of Robust
#'        MANOVA Test Statistics. J Mod Appl Stat Methods 12. Available:
#'        http://digitalcommons.wayne.edu/jmasm/vol12/iss2/4.
#'     \item Fouladi RT (1998) Type I Error Control of Two-Group Multivariate
#'        Tests on Means under Conditions of Heterogeneous Correlation
#'        Structure. Annual Meeting of the American Educational Research
#'        Association. San Diego, CA,. pp. 1 - 28. Available:
#'        http://files.eric.ed.gov/fulltext/ED420716.pdf.
#' }
#'

MANOVA <- function(obj, par.mat = NULL, cov.mat = NULL, sample.size = NULL,
    groups = NULL, boxM = TRUE) {

    hmat <- hyp_matrix(
            obj = obj,
            par.mat = par.mat,
            cov.mat = cov.mat,
            sample.size = sample.size,
            groups = groups)

    E = hmat$E
    df.e = hmat$df.e
    E.inv = hmat$E.inv
    H = hmat$H
    df.h = hmat$df.h
    HpE.inv = hmat$HpE.inv
    Se = hmat$Se
    Se.inv = hmat$Se.inv
    W = hmat$W
    W.inv = hmat$W.inv
    par = hmat$par
    Cov = hmat$Cov
    p = hmat$p
    N = hmat$N
    df.h = hmat$df.h
    m = hmat$m
    n = hmat$n
    rm(hmat); gc()

    r <- min(p, df.h) ## Necessary to compute dfs
    s <- max(p, df.h)

    if (inherits(W.inv, "try-error")) {
        ## error handling code, maybe just skip this iteration using
        war <- warning(W.inv)
    }

    ## The trace of a (square) matrix
    tr <- function(x)
                sum(diag(x), na.rm = TRUE)

    ##- ----------------------------------------------------- -#
    ### ------------------ BOX M TEST ----------------
    ##- ----------------------------------------------------- -#

    if (boxM) # TEST FOR COVARIANCE MATRIX EQUALITY
    {
        ## The Box M statistic
        M <- Reduce("+", mapply(function(x, y) (x - 1) * log(det(y)),
                m, Cov,
                SIMPLIFY = F
        ))
        M <- df.e * log(det(Se)) - M
        ## A transformation of M results in a statistic having a central
        ## chi-squared distribution with degrees of freedom,
        ## v = (J - 1)(p + 1)p/2, under the assumptions that the unit
        ## scores are independent and multivariate normal in each
        ## population ([1] pag 41)
        if (mean(m) == m[1]) {
            # If sample sizes are equal
            C <- 1 - ((2 * p^2 + 3 * p - 1) * (n + 1)) /
                            (6 * (p + 1) * df.e)
        } else {
            C <- 1 - ((2 * p^2 + 3 * p - 1) /
                    (6 * (p + 1) * (n - 1))) * (sum(1 /(m - 1)) - 1/df.e)
        }

        chisq <- C * M
        df <- (n - 1) * (p + 1) * p / 2

        chi.pvalue <- pchisq(chisq, df, lower.tail = FALSE)
        ## There are two limitations associated with the Box test.
        ## (pag 42) First, it can be sensitive to multivariate
        ## non-normality. Second, the degrees of freedom for this test
        ## can often be very large, resulting in an extremely
        ## sensitive/powerful test of the null hypothesis. Both of these
        ## limitations could lead a researcher to question the validity
        ## of the multivariate test on the vectors of means when the
        ## multivariate test is valid. As a result, researchers often do
        ## not rely heavily on the results of the test but rather rely on
        ## the robustness of the multivariate test on the equality of
        ## mean vectors when sample sizes are equal.
    }

    ##- ---------------------------------------------- -#
    ## --------- Estimation of Wilks' Lambda -----------
    ##- ---------------------------------------------- -#

    # The Wilks lambda criterion ([1], pag 49)
    WLambda <- det(E) / det(E + H)

    if (n == 2) {
        WL.F <- ((1 - WLambda) / WLambda) * (df.e - p + 1) / p
        WL.Fpval <- 1 - pf(WL.F, p, df.e - p + 1) # ([1], 3.17)
    }

    if (n == 3) {
        WL.F <- ((1 - sqrt(WLambda)) / sqrt(WLambda)) * (df.e - p + 1) / p
        WL.Fpval <- 1 - pf(WL.F, 2 * p, 2 * (df.e - p + 1)) # ([1] 3.18 )
    }

    ################################################################## #
    #### ------------ UNEQUAL COVARIANCE MATRICES -------------------
    ################################################################## #
    ## [1] pag 43. [3] pag 78
    if (!inherits(W.inv, "try-error")) {
        ##- ----------------------------------------------------- -#
        #### --------------- HOTELLING STATISTIC ---------------
        ##- ----------------------------------------------------- -#
        if (n == 2) {
            ## Hotelling T^2 (pag 39 -40):
            T2.H <- t(par[, 1] - par[, 2]) %*% Se.inv %*%
                        (par[, 1] - par[, 2])
            T2.H <- m[1] * (m[2] / N) * T2.H
            ## Se is is the p x p error covariance matrix

            ## The F-statistic:
            T2.F <- ((N - p - 1) * T2.H) / (df.e * p)

            ## has a central F distribution with degrees of freedom v1 = p and
            ## v2 = dfe - p + 1 = N - p - 1. This is true assuming independent
            ## units, bivariate normality, equal population covariance
            ## matrices, and null hypothesis, H0: par1 = par2, being true ([1]
            ## pag 40), where par1 & par2 are the parameter vectors from model
            ## 1 & 2, respectivaly. "Hotelling T^2"

            ## The F p-value:
            T2.Fpval <- 1 - pf(T2.F, p, N - p - 1)
        }

        if (boxM & n < 3 & chi.pvalue < 0.05) {
            ### ---------------------------------------------------------###
            #### -------------------------- YAO's FY -------------------
            ### ---------------------------------------------------------###
            P <- par[, 1] - par[, 2]
            b.1 <- t(P) %*% W.inv %*% A[[1]] %*% W.inv %*% P
            b.2 <- t(P) %*% W.inv %*% A[[2]] %*% W.inv %*% P
            T2.u <- T2(W, par)
            f <- 1 / (((b.1 / T2.u)^2) / (m[1] - 1) +
                            ((b.2 / T2.u)^2) / (m[2] - 1))
            v2 <- f - p + 1
            Y.F <- v2 * T2.u / (p * f)
            Y.Fpval <- 1 - pf(Y.F, p, v2)

            ### --------------------------------------------------------###
            ### -------------------- JOHANSEN's FJ ------------------
            ### --------------------------------------------------------###
            ## [3] pag 78, [4]
            W.1 <- W.inv %*% A[[1]]
            W.2 <- W.inv %*% A[[2]]
            s1 <- (tr(W.1^2) + tr(W.1)^2) / (m[1] - 1) # [4] pag 4
            s2 <- (tr(W.2^2) + tr(W.2)^2) / (m[2] - 1)
            C <- 0.5 * (s1 + s2)
            c1 <- p + 2 * C - 6 * C / (p + 1)
            J.F <- T2.u / c1 ## JOHANSEN's FJ
            v2 <- p * (p + 2) / (3 * C)
            J.Fpval <- 1 - pf(J.F, p, v2)

            ### -------------------------------------------------------###
            #### -------------------- NEL-MERWE's FNM -------------------
            ### -------------------------------------------------------###
            s0 <- tr(W^2) + tr(W)^2
            s1 <- (tr(A[[1]]^2) + tr(A[[1]])^2) / (m[1] - 1) # [4] pag 4
            s2 <- (tr(A[[2]]^2) + tr(A[[2]])^2) / (m[2] - 1) # [4] pag 4
            f <- s0 * (s1 + s2)
            v2 <- f - p + 1
            NM.F <- v2 * T2.u / (p * f) ## NEL-MERWE's FNM
            NM.Fpval <- 1 - pf(NM.F, p, v2)
        }

        ##- --------------------------------------------- -#
        #### ------------ PILLAI' S TRACE  -------------
        ##- --------------------------------------------- -#
        ## Pillai's trace. Also called Bartlett-Pillai-Hotelling trace
        ## criterion and Bartlett - Pillai Criterion.

        q <- (abs(df.h - p) - 1) / 2
        t <- (df.e - p - 1) / 2

        U <- tr(H %*% HpE.inv) # Equal to Eq 3.23 [1], pag 51

        ## This statistic can be transformed to a statistic having an F
        ## distribution with
        v1 <- r * (2 * q + r + 1)
        v2 <- r * (2 * t + r + 1) ## degrees of freedom.
        ## SPSS Statistics 17.0 Algorithms (pag 401):

        U.F1 <- (U / (r - U)) * (v2 / v1)
        ## approach at reference [1] to estimate the degrees of freedom
        v2 <- r * (df.e - p + r)
        v1 <- r * max(p, df.h)
        U.F2 <- (U / (r - U)) * (v2 / v1)
        U.F <- min(U.F1, U.F2) # the lesser

        ## The p-value for Pillai's trace:
        U.Fpval <- 1 - pf(U.F, v1, v2)

        ## ----------------------------------------------------- -#
        ## ------ Comparison of two or more than two models ------
        ## ----------------------------------------------------- -#
        ## --------------------------------------------- -#
        #### -------------- WILKS'S LAMBDA -------------
        ## --------------------------------------------- -#
        ##  F approximation
        if (n > 3) { # ([1] & [2])
            a <- df.e - (p - df.h + 1) / 2 # (3.20 [1] & [2] pag 7-8 )
            cond <- (p^2 + df.h^2 - 5)
            if (cond > 0) {
                b <- sqrt((p^2 * df.h^2 - 4) / cond)
            } else {
                b <- 1
            } # (3.21, [1] & [2] pag 7-8 )
            df.r <- (a * b - (p * df.h - 2) / 2) / (p * df.h)

            WL.F <- ((1 - WLambda^(1 / b)) / (WLambda^(1 / b))) * df.r
            WL.Fpval <- 1 - pf(WL.F, p * df.h, a * b - p * df.h / 2 + 1)
        }

        ##- ------------------------------------------------------ -#
        #### ------------- HOTELLING-LAWLEY'S TRACE  ----------
        ##- ------------------------------------------------------ -#
        ## This criterion, uses the eigenvalues of the E^-1*H matrix
        ## in a different statistic ([1], pag 51).

        V <- tr(E.inv %*% H) # Equal to Eq 3.27, pag 52 [1]
        ## eig = Re( eigen( E.inv %*% H , symmetric = FALSE)$values )
        ## V = sum( eig )
        ## This statistic can be transformed to a statistic having an F
        ## distribution with
        v1 <- r * (2 * q + r + 1)
        v2 <- 2 * (r * t + 1) ## degrees of freedom

        V.F <- V * v2 / (r * v1)
        ## The p-value for Hotelling-Lawley's trace:
        V.Fpval <- 1 - pf(V.F, v1, v2)

        ## Bartlett - Pillai criterion is about the same as that obtained
        ## using the Wilks criterion.
    } else {
        T2.u <- NA
        Y.F <- NA
        Y.Fpval <- NA
        J.F <- NA
        J.Fpval <- NA
        NM.F <- NA
        NM.Fpval <- NA
        T2 <- NA
        T2.F <- NA
        T2.Fpval <- NA
        U <- NA
        U.F <- NA
        U.Fpval <- NA
        WLambda <- NA
        WL.F <- NA
        WL.Fpval <- NA
        V <- NA
        V.F <- NA
        V.Fpval <- NA
        ## warning( "The problem probably arises because of a high degree
        ## of correlation between model parameters" )
    }

    ## - --------------------------------------------------- -#
    ## -------------- Report for two models --------------
    ## - -------------------------------------------------- - #
    if (n == 2) {
        if (boxM) {
            if (chi.pvalue > 0.05) {
                report <- list()
                report$stats <- rbind(
                    c(T2.H, T2.F, T2.Fpval),
                    c(WLambda, WL.F, WL.Fpval),
                    c(U, U.F, U.Fpval),
                    c(V, V.F, V.Fpval)
                )
                colnames(report$stats) <- c("Stat", "Fstat", "p.value")
                rownames(report$stats) <- c(
                    "Hotelling T^2",
                    "Pillai's trace",
                    "Wilks lambda",
                    "Hotelling-Lawley' trace"
                )
                report$Box.M <- data.frame(M, chisq, chi.pvalue)
                colnames(report$Box.M) <- c("Box M", "Chi^2", "p.value")
            } else {
                report <- list()
                report$stats <- rbind(
                    c(T2.H, T2.F, T2.Fpval),
                    c(WLambda, WL.F, WL.Fpval),
                    c(T2.u, Y.F, Y.Fpval),
                    c(T2.u, J.F, J.Fpval),
                    c(T2.u, NM.F, NM.Fpval),
                    c(U, U.F, U.Fpval),
                    c(V, V.F, V.Fpval)
                )
                colnames(report$stats) <- c("Stat", "Fstat", "p.value")
                rownames(report$stats) <- c(
                    "Hotelling T^2",
                    "Wilks lambda",
                    "Yao's F",
                    "Johansen's F",
                    "Nel-Merwe's F",
                    "Pillai's trace",
                    "Hotelling-Lawley' trace"
                )
                ind <- order(report$stats[, 3], decreasing = TRUE)
                report$stats <- report$stats[ind, ]

                report$Box.M <- data.frame(M, chisq, chi.pvalue)
                colnames(report$Box.M) <- c("Box M", "Chi^2", "p.value")
            }
        } else {
            report <- rbind(
                c(T2.H, T2.F, T2.Fpval),
                c(WLambda, WL.F, WL.Fpval)
            )
            colnames(report) <- c("Stat", "Fstat", "pvalue")
            rownames(report) <- c("Hotelling T^2", "Wilks lambda")
        }
    }

    ##- --------------------------------------------------------- -#
    #### --------- Report for three or more models ----------
    ##- --------------------------------------------------------- -#
    if (n >= 3) {
        if (boxM) {
            report <- list()
            report$stats <- rbind(
                c(U, U.F, U.Fpval),
                c(WLambda, WL.F, WL.Fpval),
                c(V, V.F, V.Fpval)
            )
            colnames(report$stats) <- c("Stat", "Fstat", "pvalue")
            rownames(report$stats) <- c(
                "Pillai's trace", "Wilks lambda",
                "Hotelling-Lawley' trace"
            )

            report$Box.M <- data.frame(M, chisq, chi.pvalue)
            colnames(report$Box.M) <- c("Box M", "Chi^2", "chi.pvalue")
        } else {
            report <- rbind(
                c(U, U.F, U.Fpval),
                c(WLambda, WL.F, WL.Fpval),
                c(V, V.F, V.Fpval)
            )
            colnames(report) <- c("Stat", "Fstat", "pvalue")
            rownames(report) <- c(
                "Pillai's trace", "Wilks lambda",
                "Hotelling-Lawley' trace"
            )
        }
    }
    return(report)
}

### ====================== Auxiliary functions ========================

## ------------------------------------------------------------ ##
## --------------------------- Hotelling T^2-------------------
## ------------------------------------------------------------ ##
T2 <- function(S.e, pars) {
    ## S.e is an error matrix
    S.e.inv <- suppressWarnings(try(qr.solve(S.e, tol = 1e-10),
                                    silent = TRUE
    ))
    if (inherits(S.e.inv, "try-error")) {
        S.e.inv <- solve(qr(S.e, LAPACK = TRUE))
    }
    t(pars[, 1] - pars[, 2]) %*% S.e.inv %*% (pars[, 1] - pars[, 2])
}



