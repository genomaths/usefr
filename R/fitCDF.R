## ########################################################################### #
##
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
##
## ########################################################################### #
#
#' @rdname fitCDF
#' @title Nonlinear fit of a commulative distribution function
#' @description Usually the parameter estimation of a cumulative distribution
#' function (*CDF*) are accomplished using the corresponding probability density
#' function (*PDF*). Different optimization algorithms can be used to
#' accomplished this task and different algorithms can yield different estimated
#' parameters. Hence, why not try to fit the CDF directly?
#' @param varobj A a vector, a named list, a matrix or a data.frame, containing
#' the observations from the variable for which the CDF parameters will be
#' estimated. When the argument is a matrix or a data.frame, the columns must be
#' named, carrying the objective variables.
#' @param distNames a vector of distribution numbers to select from the listed
#' below in details section, e.g. c(1:10, 15). If 'distNames' is not any of
#' current 20 named distributions, then it can be any arbitrary character
#' string, but the argument 'distf' must be given (see below).
#' @param start A named numerical vector giving the parameters to be optimized
#' with initial values or a list of numerical vectors (only when
#' \emph{\strong{varobj}} is a list, a matrix or a data.frame). This can be
#' omitted for some of the named distributions (see Details). This argument will
#' be used if provided for only one distribution. The default parameter values
#' are:
#' \enumerate{
#'     \item norm = c( mean = MEAN, sd = SD )
#'     \item lnorm = c(meanlog = mean( log1p(X), na.rm = TRUE),
#'             sdlog = sd( log1p( X ), na.rm = TRUE))
#'     \item hnorm = c(theta = sqrt(pi)/(SD*sqrt(2)))
#'     \item gnorm = c( mean = MEAN, sigma = SD, beta = 2)
#'     \item tgnorm = c( mean = MEAN, sigma = SD, beta = 2)
#'     \item laplace = c( mean = MEAN, sigma = sqrt( VAR))
#'     \item gamma = c(shape_scale(X, gg = FALSE))
#'     \item gamma3p = c(shape_scale(X, gg = FALSE), mu = 0),
#'     \item ggamma = c(shape_scale(X, gg = TRUE), mu = MIN, psi = 1)
#'     \item ggamma = c(shape_scale(X, gg = TRUE), psi = 1)
#'     \item weibull = c( shape = log( 2 ), scale = Q)
#'     \item weibull3p = c( mu = MIN, shape = log( 2 ), scale = Q)
#'     \item beta = c(shape1 = 1, shape2 = 2)
#'     \item beta3 = c(shape1 = 1, shape2 = 2, a = MIN)
#'     \item beta4 = c(shape1 = 2, shape2 = 3, a=0.9 * MIN, b=1.1 * MAX)
#'     \item bweibull = c(alpha=1, beta=2, shape = log( 2 ), scale = Q)
#'     \item gbeta = c( shape1 = 1, shape2 = 2, lambda = 1)
#'     \item rayleigh = c( sigma = SD )
#'     \item exp = c( rate = 1)
#'     \item exp2 = c( rate = 1, mu = 0)
#'     \item geom = c(prob = ifelse(MEAN > 0, 1/(1 + MEAN), 1))
#'     \item lgamma = shape_scale(log1p(X), gg = FALSE)
#'     \item lpgamma3p = c(shape_scale(log1p(X), gg = FALSE), mu = 0)
#' }
#' @param loss.fun Loss function(s) used in the regression (see
#' \href{https://en.wikipedia.org/wiki/Loss_function}{(Loss function)}). After
#' \eqn{z =  1/2 * sum((f(x) - y)^2)} we have:
#' \enumerate{
#'     \item "linear": linear function which gives a standard least squares:
#'           \eqn{loss(z) = z}.
#'     \item "huber": Huber loss, \eqn{loss(z) = ifelse(z <= 1, z, sqrt(z) -1)}.
#'     \item "smooth": Smooth approximation to the sum of residues absolute
#'           values: \eqn{loss(z) = 2*(sqrt(z + 1) - 1)}.
#'     \item "cauchy": Cauchy loss: \eqn{loss(z) = log(z + 1)}.
#'     \item "arctg": arc-tangent loss function: \eqn{loss(x) = atan(z)}.
#' }
#' @param min.val A number denoting the lower bound of the domain where CDF
#' is defined. For example, for Weibull and GGamma \strong{\emph{min.val = 0}}.
#' @param plot Logical. Default FALSE Whether to produce the plots for the best
#' fitted CDF.
#' @param plot.num The number of distributions to be plotted.
#' @param distf A character string naming a cumulative distribution function(s)
#' (CDF) present in the R session environment . For example, \strong{gamma} or
#' \strong{norm}, etc, from where, internally, we can get: density, distribution
#' function, quantile function and random generation as: \strong{dnorm},
#' \strong{pnorm}, \strong{qnorm}, and \strong{rnorm}, respectively. If the
#' function is not present in the environment, then an error will be returned.
#' It must given only if 'distNames' is not any of current 20 named
#' distributions (see  details below). Default is NULL.
#' @param only.info Logic. Default TRUE. If true, only information about the
#' parameter estimation is returned.
#' @param maxiter,maxfev,ptol Parameters to control of various aspects of the
#' Levenberg-Marquardt algorithm through function
#' \code{\link[minpack.lm]{nls.lm.control}} from *minpack.lm* package.
#' @param nls.model Logical. Whether to return the best fitted model as an
#' object from \strong{nlsModel} class. Default is FALSE. If TRUE, then
#' the estimated parameters are used new fitting with \code{\link[stats]{nls}}
#' function.
#' @param algorithm Only if \strong{nls.model} = TRUE. The same as for
#' \code{\link[stats]{nls}} function.
#' @param xlabel (Optional) Label for variable \strong{\emph{varobj}}.
#' Default is \emph{xlabel = "x"}.
#' @param mar,mgp,las,cex.main (Optional) Graphical parameters (see
#' \code{\link[graphics]{par}}).
#' @param cex.text,cex.point Numerical value to scale text and points.
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to use,
#' i.e. at most how many child processes will be run simultaneously (see
#' \code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
#' for Linux OS).
#' @param ... (Optional) Further graphical parameters (see
#' \code{\link[graphics]{par}}). Graphical parameter will simultaneously affect
#' all the plots.
#' @param verbose Logic. If TRUE, prints the function log to stdout
#' @details The nonlinear fit (NLF) problem for CDFs is addressed with
#' Levenberg-Marquardt algorithm implemented in function
#' \code{\link[minpack.lm]{nls.lm}} from package *minpack.lm*. The Stein's rho
#' for adjusted R squared (rho) is applied as an estimator of the average
#' cross-validation predictive power [1]. This function is inspired in a script
#' for the function \code{\link[propagate]{fitDistr}} from the package propagate
#' [2]. Some parts or script ideas from function
#' \code{\link[propagate]{fitDistr}} are used, but here we to estimate CDF and
#' not the PDF as in the case of "\code{\link[propagate]{fitDistr}}. More
#' informative results are given now. The studentized residuals are provided as
#' well. The list (so far) of possible CDFs is:
#'     \enumerate{
#'         \item Normal \href{https://goo.gl/xaEAdT}{(Wikipedia)}
#'         \item Log-normal \href{https://goo.gl/a7MtYq}{(Wikipedia)}. This
#'         This function is set to fit \eqn{log(1+x)}. Users can transform their
#'         variable by themself and then try the fitting to Normal
#'         distribution.
#'         \item Half-normal \href{https://goo.gl/yxMF6T}{(Wikipedia)}. An
#'             Alternatively using a scaled precision (inverse of the variance)
#'             parametrization (to avoid issues if \eqn{\sigma} is near zero),
#'             obtained by setting \eqn{\theta=sqrt(\pi)/\sigma*sqrt(2)}.
#'         \item Generalized Normal \href{https://goo.gl/EPk8mH}{(Wikipedia)}
#'         \item T-Generalized Normal [3].
#'         \item Laplace \href{https://goo.gl/fCykV9}{(Wikipedia)}
#'         \item Gamma \href{https://goo.gl/cYkvar}{(Wikipedia)}
#'         \item 3P Gamma [4].
#'         \item Generalized 4P Gamma [4]
#'                 \href{https://goo.gl/1n4kpW.}{(Wikipedia)}
#'         \item Generalized 3P Gamma [4].
#'         \item Weibull \href{https://goo.gl/WMXmQP}{(Wikipedia)}
#'         \item 3P Weibull \href{https://goo.gl/WMXmQP}{(Wikipedia)}
#'         \item Beta \href{https://goo.gl/893wzR}{(Wikipedia)}
#'         \item 3P Beta \href{https://goo.gl/893wzR}{(Wikipedia)}
#'         \item 4P Beta \href{https://goo.gl/893wzR}{(Wikipedia)}
#'         \item Beta-Weibull \href{https://goo.gl/dpaG8h}{ReliaWiki}
#'         \item Generalized Beta \href{https://goo.gl/UcVsct}{(Wikipedia)}
#'         \item Rayleigh \href{https://goo.gl/d9b3zv}{(Wikipedia)}
#'         \item Exponential \href{https://goo.gl/stVsi7}{(Wikipedia)}
#'         \item 2P Exponential \href{https://goo.gl/stVsi7}{(Wikipedia)}
#'         \item Geometric \href{https://is.gd/94HW4w}{(Wikipedia)}
#'         \item Log-Gamma \href{https://is.gd/kMQVxX}{(Mathematica)}
#'         \item Log-Gamma 3P \href{https://is.gd/kMQVxX}{(Mathematica)}
#'     }
#'
#' Where, shape_scale function is an internal function that can be
#' retrieve by typing: usefr:::shape_scale.
#' @return After return the plots, a list with following values is provided:
#'     \itemize{
#'         \item aic: Akaike information creterion
#'         \item fit: list of results of fitted distribution, with parameter
#'             values
#'         \item bestfit: the best fitted distribution according to AIC
#'         \item fitted: fitted values from the best fit
#'         \item rstudent: studentized residuals
#'         \item residuals: residuals
#'     }
#'     After cdf = fitCDF( varobj, ...), attributes( cdf$bestfit ) shows the
#'     list of objects carry on cdf$bestfit:
#'      \itemize{
#'         \item names: "par" "hessian" "fvec" "info" "message" "diag" "niter"
#'                      "rsstrace"  "deviance"
#'         \item class: "nls.lm"
#'      }
#'
#' And fitting details can be retrieved with summary(cdf$bestfit)
#'
#' @importFrom numDeriv grad
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @importFrom utils flush.console
#' @importFrom mixdist weibullpar
#' @importFrom stats var sd quantile ecdf pgamma pnorm pbeta pexp pweibull
#' @importFrom stats plnorm na.omit splinefun qqnorm qqline pgeom
#' @importFrom graphics par grid lines mtext abline
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @seealso \code{\link[MASS]{fitdistr}} and \code{\link{fitMixDist}} and
#'     for goodness-of-fit: \code{\link{mcgoftest}}.
#' @references
#'   \enumerate{
#'     \item Stevens JP. Applied Multivariate Statistics for the Social
#'           Sciences. Fifth Edit. Routledge Academic; 2009.
#'     \item Andrej-Nikolai Spiess (2014). propagate: Propagation of
#'           Uncertainty. R package version 1.0-4.
#'           http://CRAN.R-project.org/package=propagate
#'     \item Abramowitz, M. and Stegun, I. A. (1972) Handbook of Mathematical
#'           Functions. New York: Dover. Chapter 6: Gamma and Related Functions.
#'     \item Hand-book on STATISTICAL DISTRIBUTIONS for experimentalists
#'           (pag 73) by Christian Walck. Particle Physics Group Fysikum.
#'           University of Stockholm (e-mail: walck@physto.se).
#'   }
#' @examples
#' set.seed(1230)
#' x1 <- rnorm(10000, mean = 0.5, sd = 1)
#' cdfp <- fitCDF(x1, distNames = "Normal", plot = FALSE)
#' summary(cdfp$bestfit)
#'
#' ## Add some cosmetics to the plots
#' cdfp <- fitCDF(x1,
#'     distNames = "Normal", xlabel = "My Nice Variable Label",
#'     plot = T, font.lab = 3, font = 2, font.axis = 2, family = "serif",
#'     cex.lab = 1.3, cex.axis = 1.3
#' )
#'
#' ## Fitting a Weibull distribution with 3 paramaters
#' x1 <- rweibull3p(1000, shape = 0.5, scale = 1, mu = 0.1)
#' cdfp <- fitCDF(x1,
#'     distNames = "3P Weibull",
#'     xlabel = "My Nice Variable Label",
#'     plot = T, font.lab = 3, font = 2, font.axis = 2, family = "serif",
#'     cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.1,
#'     mgp = c(2.5, 1, 0)
#' )
#' @aliases fitCDF
setGeneric(
    "fitCDF",
    def = function(varobj,
    ...) {
        standardGeneric("fitCDF")
    }
)

#' @rdname fitCDF
#' @aliases fitCDF
#' @export
setMethod(
    "fitCDF", signature(varobj = "numeric"),
    function(varobj,
    distNames,
    plot = FALSE,
    plot.num = 1L,
    distf = NULL,
    start = NULL,
    loss.fun = c(
        "linear", "huber", "smooth",
        "cauchy", "arctg"
    ),
    min.val = NULL,
    only.info = FALSE,
    maxiter = 1024,
    maxfev = 1e+5,
    ptol = 1e-12,
    nls.model = FALSE,
    algorithm = "default",
    xlabel = "x",
    mar = c(4, 4, 3, 1),
    mgp = c(2.5, 0.6, 0),
    las = 1,
    cex.main = 1,
    cex.text = 0.8,
    cex.point = 0.5,
    verbose = TRUE,
    ...) {
        if (is.numeric(distNames)) {
            distNames <- as.integer(distNames)
            if (any(distNames > 23)) {
                stop("*** 'distNames' must be a string or a number <= 23")
            }
        }

        loss.fun <- match.arg(loss.fun)

        # "count" function from "propagate"
        counter <- function(i) {
            if (i %% 10 == 0) {
                cat(i)
            } else {
                cat(".")
            }
            if (i %% 50 == 0) {
                cat("\n")
            }
            flush.console()
        }

        options(warn = -1)
        if (is.vector(varobj)) {
            X <- sort(varobj)
        } else {
            stop("varobj must be a numeric vector!")
        }

        if (!is.null(min.val)) {
            X <- X[which(X > min.val)]
        }

        MEAN <- mean(X, na.rm = TRUE)
        VAR <- var(X, na.rm = TRUE)
        SD <- sd(X, na.rm = TRUE)
        MIN <- min(X, na.rm = TRUE)
        MAX <- max(X, na.rm = TRUE)
        Q <- unname(quantile(X, 0.632, na.rm = TRUE))
        Fy <- ecdf(X)
        pX <- Fy(X)

        distNAMES <- c(
            "Normal",
            "Log-normal",
            "Half-Normal",
            "Generalized Normal",
            "T-Generalized Normal",
            "Laplace",
            "Gamma",
            "3P Gamma",
            "Generalized 4P Gamma",
            "Generalized 3P Gamma",
            "Weibull",
            "3P Weibull",
            "Beta",
            "3P Beta",
            "4P Beta",
            "Beta-Weibull",
            "Generalized Beta",
            "Rayleigh",
            "Exponential",
            "2P Exponential",
            "Geometric",
            "Log-Gamma",
            "Log-Gamma 3P"
        )

        funLIST <- list(
            pnorm,
            plnorm,
            phnorm,
            pgnorm,
            ptgnorm,
            plaplace,
            pgamma,
            pgamma3p,
            pggamma,
            pggamma,
            pweibull,
            pweibull3p,
            pbeta,
            pbeta3,
            pbeta4,
            pbweibull,
            pgbeta,
            prayleigh,
            pexp,
            pexp2,
            pgeom,
            plgamma,
            plgamma3p
        )

        ## ------------------------- parLIST ---------------------------
        parLIST <- list(
            norm = c(mean = MEAN, sd = SD),
            lnorm = c(
                meanlog = mean(log1p(X), na.rm = TRUE),
                sdlog = sd(log1p(X), na.rm = TRUE)
            ),
            hnorm = c(theta = sqrt(pi) / (SD * sqrt(2)), mu = 0),
            gnorm = c(
                mean = MEAN,
                sigma = SD,
                beta = 2
            ),
            tgnorm = c(
                mean = MEAN,
                sigma = SD,
                beta = 2
            ),
            laplace = c(mean = MEAN, sigma = sqrt(VAR)),
            gamma = shape_scale(X, gg = FALSE),
            gamma3p = c(shape_scale(X, gg = FALSE), mu = 0),
            ggamma = c(shape_scale(X, gg = TRUE), mu = MIN, psi = 1),
            ggamma = c(shape_scale(X, gg = TRUE), psi = 1),
            weibull = weibullpars(mu = MEAN, sigma = SD),
            weibull3p = c(weibullpars(mu = MEAN, sigma = SD), mu = MIN),
            beta = beta_start_par(MEAN, VAR),
            beta3 = c(
                beta_start_par(MEAN, VAR),
                mu = MIN
            ),
            beta4 = c(
                beta_start_par(MEAN, VAR),
                mu = 0.9 * MIN,
                b = 1.1 * MAX
            ),
            bweibull = c(
                alpha = 1,
                beta = 2,
                shape_scale(X, gg = FALSE)
            ),
            gbeta = c(
                shape1 = 1,
                shape2 = 2,
                lambda = 1
            ),
            rayleigh = c(sigma = SD),
            exp = c(rate = 1),
            exp2 = c(rate = 1, mu = 0),
            geom = c(prob = ifelse(MEAN > 0, 1 / (1 + MEAN), 1)),
            lgamma = c(shape = 1, scale = 1),
            lgamma3p = c(shape = 1, scale = 1, mu = 0)
        )
        ## ------------- Distribution names ----------------------
        if (is.character(distNames)) {
            elemt <- all(is.element(distNames, distNAMES))
        } else
        if (is.numeric(distNames)) {
            elemt <- all(is.element(distNames, seq_along(distNAMES)))
            distNames <- as.integer(distNames)
        }
        if (missing(distNames)) {
            distNames <- seq_along(distNAMES)
        }

        if (length(distNames) < 24 && elemt) {
            if (is.character(distNames)) {
                distNames <- as.integer(na.omit(match(distNames, distNAMES)))
            }

            if (sum(!is.element(distNames, 1:length(distNAMES))) > 0) {
                stop(
                    "At least one CDF is not found between the",
                    " possible selections \n"
                )
            }

            if (is.integer(distNames)) {
                distnms <- distNAMES
                distNAMES <- distNAMES[distNames]
                funLIST <- funLIST[distNames]
                distr <- distr[distNames]

                if (!is.null(start)) {
                    parLIST <- list(start)
                } else {
                    parLIST <- parLIST[distNames]
                }
            }
        }

        if (is.character(distNames) && !elemt) {
            if (is.null(distf)) {
                stop("*** A user defined distribution function must be given")
            }
            distp <- try(sapply(distf, function(s) match.fun(paste0("p", s))),
                silent = TRUE
            )
            if (inherits(distp, "try-error")) {
                stop(
                    "*** 'distf' must be a symbol to call a commulative ",
                    "distribution function e.g, pnorm, pgamma"
                )
            }

            funLIST <- distp
            if (is.null(start)) {
                stop("*** 'start' parameter values must be provided")
            } else {
                parLIST <- if (is.list(start)) start else list(start)
            }

            if (length(distf) != length(distNames)) {
                stop(
                    "*** The lengths of 'distNames' and 'distf' arguments",
                    " must be equal"
                )
            }

            distNAMES <- distNames
            distnms <- distNAMES
            rfunLIST <- paste0("r", distf)
            qfunLIST <- paste0("q", distf)
        }

        if (!is.null(distf) && is.character(distf)) {
            funName <- paste0("p", distf)
        }

        if (elemt) {
            funName <- paste0("p", distr)
        }

        ## ------------------ Fitting ------------------------

        fitLIST <- vector("list", length = length(distNAMES))
        AICS <- rep(NA, length(distNAMES))

        for (i in seq_along(distNAMES)) {
            if (verbose) {
                message("\n*** Fitting ", distNAMES[i],
                    " distribution ... ",
                    sep = ""
                )
            }

            PARS <- matrix(parLIST[[i]], nrow = 1)
            colnames(PARS) <- names(parLIST[[i]])
            PARS <- PARS[1, ]

            if (verbose) {
                counter(i)
            }
            FIT <- try(nls.lm(
                par = PARS,
                fn = optFun,
                probfun = funLIST[[i]],
                quantiles = X,
                prob = pX,
                loss.fun = loss.fun,
                control = nls.lm.control(
                    maxiter = maxiter,
                    maxfev = maxfev,
                    ptol = ptol
                )
            ),
            silent = TRUE
            )

            if (inherits(FIT, "try-error")) {
                message("\n", FIT, "Function model: ", distNAMES[i])
            }

            if (nls.model && !inherits(FIT, "try-error")) {
                pars <- names(coef(FIT))
                pars <- paste(c("q", pars), collapse = ",")
                formula <- as.formula(
                    paste0(
                        "Y ~ ", funName[i],
                        "(", pars, ")"
                    )
                )

                FIT1 <- try(nls(
                    formula,
                    data = data.frame(q = X, Y = pX),
                    start = as.list(coef(FIT)),
                    control = list(maxiter = maxiter, tol = ptol)
                ),
                silent = TRUE
                )

                if (inherits(FIT1, "try-error")) {
                    FIT1 <- try(nlsLM(
                        formula,
                        data = data.frame(q = X, Y = pX),
                        start = as.list(coef(FIT)),
                        control = list(maxiter = maxiter, ptol = ptol)
                    ),
                    silent = TRUE
                    )
                }

                if (!inherits(FIT1, "try-error")) {
                    FIT <- FIT1
                }
            }
            if (!inherits(FIT, "try-error")) {
                fitLIST[[i]] <- FIT

                evalLIST <- as.list(coef(FIT))
                evalLIST$q <- X
                fitted <- do.call(funLIST[[i]], evalLIST)
                RESIDUAL <- pX - fitted

                AICS[i] <- tryCatch(
                    AICmodel(FIT, residuals = RESIDUAL),
                    error = function(e) {
                        NA
                    }
                )
                if (verbose) {
                    cat("Fitting Done.\n")
                }
            } else {
                fitLIST[[i]] <- "Error"
                if (verbose) {
                    cat("Error!\n")
                }
            }
        }

        ## -------------- Ordering fitting results -----------------

        ORDER <- order(AICS)
        aicDAT <- data.frame(Distribution = distNAMES, AIC = AICS)
        aicDAT <- aicDAT[ORDER, ]

        distNAMES <- distNAMES[ORDER]
        funLIST <- funLIST[ORDER]
        fitLIST <- fitLIST[ORDER]
        funName <- funName[ORDER]
        bestFIT <- fitLIST[[1]]

        names(fitLIST) <- funName

        evalLIST <- as.list(coef(fitLIST[[1]]))
        evalLIST$q <- X
        fitted <- do.call(funLIST[[1]], evalLIST)
        RESIDUAL <- pX - fitted

        if (!nls.model) {
            bestFIT$fvec <- RESIDUAL
        }

        rfunLIST <- rfunLIST[match(distNAMES, distnms)]
        qfunLIST <- qfunLIST[match(distNAMES, distnms)]

        if (only.info) {
            rho <- Stein_rho(fit = bestFIT, varobj = pX)
            res <- list(bestfit = bestFIT, gof = rho, AIC = aicDAT)
            cdf <- substring(funName[1], 2, nchar(funName[1]))

            res <- structure(
                list(
                    aic = aicDAT,
                    bestfit = bestFIT,
                    fit = fitLIST,
                    fitted = fitted,
                    info = distNAMES[1],
                    rstudent = NA,
                    gof = c(res$gof, R.Cross.val = NA, AIC = aicDAT[1, 2]),
                    cdf = cdf
                ),
                class = "CDFmodel"
            )
        } else {
            ## ------------------------ Start Graphics --------------------- ##
            if (plot) {
                opar <- par()
                for (k in 1:min(plot.num, length(distNames))) {
                    cat(
                        " * Estimating Studentized residuals for",
                        distNAMES[k],
                        "distribution\n"
                    )
                    # SEL <- ORDER[ k ]
                    FITs <- fitLIST[[k]]
                    evalLIST <- as.list(FITs$par)
                    evalLIST$q <- X
                    evalY <- do.call(funLIST[[k]], evalLIST)
                    resids <- pX - evalY

                    rstudent <- try(rstudents(
                        model = FITs,
                        varobj = X,
                        fun = funName[k],
                        residuals = resids
                    ),
                    silent = TRUE
                    )
                    if (inherits(rstudent, "try-error")) {
                        rstudent <- NA
                    }
                    outliers <- sum(abs(rstudent) > 2, na.rm = TRUE)

                    cat(" * Plots for", distNAMES[k], "distribution...\n")
                    par(
                        mfrow = c(2, 2),
                        mar = mar,
                        mgp = mgp,
                        las = las,
                        ...
                    )
                    plot(
                        Fy,
                        verticals = TRUE,
                        panel.first = {
                            points(0,
                                0,
                                pch = 16,
                                cex = 1e6,
                                col = "grey95"
                            )
                            grid(col = "white", lty = 1)
                        },
                        col = "blue",
                        pch = 20,
                        bty = "n",
                        xlab = xlabel,
                        ylab = "CDF",
                        cex = cex.point,
                        main = paste(distNAMES[k], "Distribution"),
                        cex.main = cex.main
                    )
                    lines(
                        evalLIST$q,
                        evalY,
                        col = 2,
                        lty = 2,
                        lwd = 2
                    )
                    mtext(
                        text = paste("AIC =", round(aicDAT$AIC[k], 3)),
                        cex = cex.text
                    )

                    ## PP-plot
                    # par(mar = c( 5, 2, 2.2, 1 ) + 0.2, mgp = c(1.2, 0.4,0),
                    #     las = 1)
                    plot(
                        pX,
                        evalY,
                        panel.first = {
                            points(0,
                                0,
                                pch = 16,
                                cex = 1e6,
                                col = "grey95"
                            )
                            grid(col = "white", lty = 1)
                        },
                        col = "blue",
                        bty = "n",
                        main = "P-P Plot",
                        pch = 20,
                        cex = cex.point,
                        xlab = "Empirical CDF",
                        ylab = "Theoretical CDF",
                        cex.main = cex.main
                    )
                    abline(0, 1, col = "red", lwd = 2) # Reference line y = x
                    pars <- FITs$par
                    par_n <- names(pars)
                    pars <- unlist(pars)
                    pars <-
                        t(cbind(par_n, format(round(
                            pars, 3
                        ), 3)))
                    mtext(
                        text = paste(pars, collapse = " "),
                        cex = cex.text
                    )

                    # par(mar = c( 5, 3, 1, 1 ) + 0.07, mgp = c(1.2, 0.4, 0))
                    plot(
                        X,
                        rstudent,
                        panel.first = {
                            points(0,
                                0,
                                pch = 16,
                                cex = 1e6,
                                col = "grey95"
                            )
                            grid(col = "white", lty = 1)
                        },
                        bty = "n",
                        pch = 20,
                        xlab = xlabel,
                        ylab = "Studentized residuals",
                        main = "Outliers",
                        col = "blue",
                        cex = cex.point,
                        cex.main = cex.main
                    )
                    abline(
                        h = 2,
                        col = "red",
                        lwd = 2,
                        lty = 2
                    ) # Reference line y = x
                    abline(
                        h = -2,
                        col = "red",
                        lwd = 2,
                        lty = 2
                    ) # Reference line y = x
                    mtext(
                        text = paste(
                            "Number of outliers (|st| > 2): ",
                            outliers,
                            collapse = " "
                        ),
                        cex = cex.text
                    )

                    rdistr <- try(get(rfunLIST[[k]],
                        mode = "function",
                        envir = parent.frame()
                    ),
                    silent = TRUE
                    )
                    qdistr <- try(get(qfunLIST[[k]],
                        mode = "function",
                        envir = parent.frame()
                    ),
                    silent = TRUE
                    )
                    ## --- Q-Q plot
                    # par( mar = c( 5, 2, 1, 1) + 0.07, mgp = c(1.2,0.4,0))
                    if (!inherits(rdistr, "try-error") &&
                        !inherits(qdistr, "try-error")) {
                        r <- rValues(rfunLIST[[k]], FITs, 1e4)
                        xl <- c(min(r, na.rm = TRUE), max(r, na.rm = TRUE))
                        qqplot(
                            x = r,
                            y = X,
                            panel.first = {
                                points(
                                    0,
                                    0,
                                    pch = 16,
                                    cex = 1e6,
                                    col = "grey95"
                                )
                                grid(col = "white", lty = 1)
                            },
                            bty = "n",
                            xlim = xl,
                            ylim = xl,
                            col = "blue",
                            pch = 20,
                            cex = cex.point,
                            cex.main = cex.main,
                            qtype = 6,
                            main = "Q-Q Plot",
                            xlab = "Theoretical Quantiles",
                            ylab = "Empirical Quantiles"
                        )
                        # Reference line y = x
                        qqline(
                            y = X,
                            distribution = function(p) {
                                qValues(p, qfunLIST[[k]], FITs)
                            },
                            xlim = xl,
                            ylim = xl,
                            cex = cex.point,
                            col = "red",
                            lwd = 2,
                            qtype = 6
                        )
                    } else {
                        qqnorm(
                            rstudent,
                            panel.first = {
                                points(
                                    0,
                                    0,
                                    pch = 16,
                                    cex = 1e6,
                                    col = "grey95"
                                )
                                grid(col = "white", lty = 1)
                            },
                            bty = "n",
                            main = "Q-Q Plot",
                            col = "blue",
                            pch = 20,
                            cex = cex.point,
                            cex.main = cex.main,
                            qtype = 6
                        )
                        # Reference line y = x
                        qqline(
                            rstudent,
                            col = "red",
                            lwd = 2,
                            qtype = 6,
                            cex.main = cex.main,
                            cex = cex.point
                        )
                    }
                }
                par(opar)
                names(fitLIST) <- distNAMES
                fitLIST <- fitLIST[as.character(aicDAT$Distribution)]

                rho <- Stein_rho(fit = fitLIST[[1]], varobj = pX)
                cdf <- substring(funName[1], 2, nchar(funName[1]))

                res <- structure(
                    list(
                        aic = aicDAT,
                        bestfit = bestFIT,
                        fit = fitLIST,
                        fitted = fitted,
                        info = distNAMES[1],
                        rstudent = NA,
                        cdf = cdf,
                        gof = rho,
                        formula = NA
                    ),
                    class = "CDFmodel"
                )

                res <- rstudents(
                    model = res, varobj = X,
                    residuals = RESIDUAL
                )

                pars <- names(coef(res$bestfit))
                pars <- paste(c("q", pars), collapse = ",")
                res$formula <- as.formula(
                    paste0(
                        "Y ~ ", funName[1],
                        "(", pars, ")"
                    )
                )

                ## ================ Cross FIT =============
                if (identical(cdf, "lnorm") ||
                    identical(cdf, "lgamma") ||
                    identical(cdf, "lgamma3p")) {
                    logx <- TRUE
                }
                else
                    logx <- FALSE

                cross_val <- cdf_crossval(
                    model = res,
                    q = X,
                    logx = logx,
                    min.val = min.val,
                    maxiter = maxiter,
                    ptol = ptol,
                    minFactor = 1e-6
                )

                res$gof <- c(res$gof, cross_val, AIC = aicDAT[1, 2])

                ## --------------------- End Graphics ------------------ ##
            } else {
                names(fitLIST) <- distNAMES

                rho <- Stein_rho(fit = bestFIT, varobj = pX)

                fitLIST <- fitLIST[as.character(aicDAT$Distribution)]
                cdf <- substring(funName[1], 2, nchar(funName[1]))

                res <- structure(
                    list(
                        aic = aicDAT,
                        bestfit = bestFIT,
                        fit = fitLIST,
                        fitted = fitted,
                        info = distNAMES[1],
                        rstudent = NA,
                        gof = rho,
                        cdf = cdf,
                        formula = NA
                    ),
                    class = "CDFmodel"
                )

                pars <- names(coef(res$bestfit))
                pars <- paste(c("q", pars), collapse = ",")
                res$formula <- as.formula(
                    paste0(
                        "Y ~ ", funName[1],
                        "(", pars, ")"
                    )
                )

                ## ================ Cross FIT =============
                if (identical(cdf, "lnorm") ||
                    identical(cdf, "lgamma") ||
                    identical(cdf, "lgamma3p")) {
                    logx <- TRUE
                }
                else
                    logx <- FALSE

                cross_val <- cdf_crossval(
                    model = res,
                    q = X,
                    logx = logx,
                    min.val = min.val,
                    maxiter = maxiter,
                    ptol = ptol,
                    minFactor = 1e-6
                )

                res$gof <- c(res$gof, cross_val, AIC = aicDAT[1, 2])

                res <- rstudents(
                    model = res, varobj = X,
                    residuals = RESIDUAL
                )

                if (verbose) {
                    cat("** Done ***\n")
                }
            }
        }
        return(res)
    }
)

setClassUnion(
    "list_OR_matrix_OR_dataframe",
    c("list", "matrix", "data.frame")
)

#' @rdname fitCDF
#' @aliases fitCDF
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @export
setMethod(
    "fitCDF", signature(varobj = "list_OR_matrix_OR_dataframe"),
    function(varobj,
    distNames,
    plot = FALSE,
    plot.num = 1,
    distf = NULL,
    start = NULL,
    loss.fun = c(
        "linear", "huber", "smooth",
        "cauchy", "arctg"
    ),
    only.info = FALSE,
    maxiter = 1024,
    maxfev = 1e+5,
    ptol = 1e-12,
    xlabel = "x",
    mar = c(4, 4, 3, 1),
    mgp = c(2.5, 0.6, 0),
    las = 1,
    cex.main = 1,
    cex.text = 0.8,
    cex.point = 0.5,
    num.cores = 1L,
    tasks = 0L,
    verbose = TRUE,
    ...) {

        ## -------------- Setting parallel computation ----------------- #
        progressbar <- FALSE
        if (verbose) progressbar <- TRUE
        if (Sys.info()["sysname"] == "Linux") {
            bpparam <- MulticoreParam(
                workers = num.cores, tasks = tasks,
                progressbar = progressbar
            )
        } else {
            bpparam <- SnowParam(
                workers = num.cores, type = "SOCK",
                progressbar = progressbar
            )
        }
        ## -------------------------------------------------------------- #

        var_nms <- names(varobj)
        if (is.list(varobj)) {
            num_samp <- length(varobj)
        }
        if (inherits(varobj, c("matrix", "data.frame"))) {
            num_samp <- ncol(varobj)
        }

        if (!is.numeric(distNames) && !is.character(distNames)) {
            stop(
                "*** 'distNames' argument must be a 'character'",
                " or numeric."
            )
        }

        if (is.character(distNames)) {
            distNames <- list(distNames)
            if (length(distNames) < num_samp) {
                distNames <- rep(distNames, num_samp)
            }
        }

        if (is.numeric(distNames)) {
            distNames <- as.integer(distNames)
            if (any(distNames > 21)) {
                stop("*** 'distNames' must be a string or a number < 21")
            }
            distNames <- list(distNames)
            if (length(distNames) < num_samp) {
                distNames <- rep(distNames, num_samp)
            }
        }

        if (!is.character(distf) && !is.null(distf)) {
            stop(
                "*** 'distf' argument must be a 'character' class",
                "object."
            )
        }

        if (is.character(distf) && length(distf) < num_samp) {
            distf <- list(distf)
            distf <- rep(distf, num_samp)
        }
        if (is.list(start)) {
            if (!all(sapply(start, is.numeric))) {
                stop("*** 'start' must be a list of numerical vectors")
            }
        }

        if (is.numeric(start)) {
            start <- list(start)
            if (length(start) < num_samp) {
                start <- rep(start, num_samp)
            }
        }

        if (is.list(start) && length(start) < num_samp) {
            start <- list(start)[rep(1, num_samp)]
        }

        nms <- (length(loss.fun) < num_samp)
        if (num_samp > 1 && nms) {
            loss.fun <- rep(loss.fun, num_samp)
        }

        if (inherits(varobj, "list")) {
            res <- bplapply(seq_along(varobj),
                function(k) {
                    fitCDF(
                        varobj = varobj[[k]],
                        distNames = distNames[[k]],
                        plot = plot,
                        plot.num = plot.num,
                        distf = distf[[k]],
                        start = start[[k]],
                        loss.fun = loss.fun[k],
                        only.info = only.info,
                        maxiter = maxiter,
                        maxfev = maxfev,
                        ptol = ptol,
                        xlabel = xlabel,
                        mar = mar,
                        mgp = mgp,
                        las = las,
                        cex.main = cex.main,
                        cex.text = cex.text,
                        cex.point = cex.point,
                        num.cores = num.cores,
                        tasks = tasks,
                        verbose = verbose,
                        ...
                    )
                },
                BPPARAM = bpparam
            )
        }

        if (inherits(varobj, c("matrix", "data.frame"))) {
            res <- bplapply(seq_len(ncol(varobj)),
                function(k) {
                    fitCDF(
                        varobj = varobj[, k],
                        distNames = distNames[[k]],
                        plot = plot,
                        plot.num = plot.num,
                        distf = distf[[k]],
                        start = start[[k]],
                        loss.fun = loss.fun[k],
                        only.info = only.info,
                        maxiter = maxiter,
                        maxfev = maxfev,
                        ptol = ptol,
                        xlabel = xlabel,
                        mar = mar,
                        mgp = mgp,
                        las = las,
                        cex.main = cex.main,
                        cex.text = cex.text,
                        cex.point = cex.point,
                        num.cores = num.cores,
                        tasks = tasks,
                        verbose = FALSE,
                        ...
                    )
                },
                BPPARAM = bpparam
            )
        }

        names(res) <- var_nms
        res$AICs <- summary_aic(res)
        res <- structure(res, class = "CDFmodelList")
        return(res)
    }
)


## ============================= Auxiliary functions ========================= #

optFun <- function(par, probfun, quantiles, prob, eval = FALSE, loss.fun) {
    if (identical(probfun, plnorm) || identical(probfun, plgamma) ||
        identical(probfun, plgamma3p)) {
        quantiles <- quantiles[ quantiles > 0 ]
        quantiles <- log(quantiles)
        Fy <- ecdf(quantiles)
        prob <- Fy(quantiles)
    }

    START <- as.list(par)
    START$q <- quantiles
    EVAL <- try(do.call(probfun, START), silent = TRUE)
    if (inherits(EVAL, "try-error")) {
        return(NA)
    }
    EVAL[is.nan(EVAL)] <- 0
    if (eval) {
        return(EVAL)
    } else {
        # nls.lm will minimize the sum of squares of vector RSS
        RSS <- 1 / 2 * (prob - EVAL)^2
        RSS <- loss(z = RSS, fun = loss.fun)
        return(RSS)
    }
}

loss <- function(z,
    fun = c(
        "linear", "huber",
        "smooth", "cauchy",
        "arctg"
    )) {
    fun <- match.arg(fun)
    fun <- switch(fun,
        linear = z,
        huber = ifelse(z <= 1, z, sqrt(z) - 1),
        smooth = 2 * (sqrt(z + 1) - 1),
        cauchy = log(z + 1),
        arctg = atan(z)
    )
    return(fun)
}

## ---- parameter estimation of Gamma Dist
shape_scale <- function(x, gg = TRUE) {
    n <- length(x)
    s1 <- n * sum(x * log(x)) - sum(log(x)) * sum(x)
    alpha <- n * sum(x) / s1
    scale <- s1 / n^2

    if (gg) {
        return(c(alpha = alpha, scale = scale))
    } else {
        return(c(shape = alpha, scale = scale))
    }
}

weibullpars <- function(mu, sigma) {
    return(weibullpar(mu = mu, sigma = sigma)[-3])
}

distr <- c(
    "norm",
    "lnorm",
    "hnorm",
    "gnorm",
    "tgnorm",
    "laplace",
    "gamma",
    "gamma3p",
    "ggamma",
    "ggamma",
    "weibull",
    "weibull3p",
    "beta",
    "beta3",
    "beta4",
    "bweibull",
    "gbeta",
    "rayleigh",
    "exp",
    "exp2",
    "geom",
    "lgamma",
    "lgamma3p"
)


rfunLIST <- paste0("r", distr)
qfunLIST <- paste0("q", distr)

rValues <- function(distn, fit, n) {
    evalLIST <- as.list(fit$par)
    evalLIST$n <- n
    return(do.call(distn, evalLIST))
}

qValues <- function(p, distn, fit) {
    evalLIST <- as.list(fit$par)
    evalLIST$p <- p
    return(do.call(distn, evalLIST))
}


## -------------- Stein adjusted R square ---------------

Stein_rho <- function(fit, varobj) {
    m <- length(fit$par)
    n <- length(varobj)
    Adj.R.Square <- 1 - (deviance(fit) / ((n - m) *
        var(varobj, use = "everything")))
    if (m > 2) {
        rho <- ((n - 1) / (n - 4)) * ((n - 2) / (n - 5)) * ((n + 1) / n)
    } else {
        rho <- ((n - 1) / (n - 3)) * ((n - 2) / (n - 4)) * ((n + 1) / n)
    }
    rho <- 1 - rho * (1 - Adj.R.Square)
    rho <- ifelse(is.na(rho) || rho < 0, 0, rho)
    return(c(Adj.R.Square = Adj.R.Square, rho = rho))
}

## ---------------------------------------------------------------------- #
## *** Generalized normal CDF *** #
## https://en.(Wikipedia).org/wiki/Generalized_normal_distribution

pgnorm <- function(q,
    mean = 0,
    sigma = 1,
    beta = 1) {
    y <- (abs(q - mean) / sigma)^beta
    # 1/2 + sign( q - mean ) * pgamma( y, 1/beta )/( 2*gamma( 1/beta ) )
    1 / 2 + sign(q - mean) * pgamma(y, 1 / beta) / 2
}

## *** Thermodynamic based Generalized normal CDF *** #

ptgnorm <- function(q,
    mean = 0,
    sigma = 1,
    beta = 1) {
    # pgamma is closely related to the incomplete gamma function.
    # As defined by Abramowitz and Stegun 6.5.1 (and by 'Numerical Recipes')
    # this is: P(a,x) = 1/Gamma(a) integral_0^x t^(a-1) exp(-t) dt
    # P(a, x) is pgamma(x, a).
    # Abramowitz, M. and Stegun, I. A. (1972) Handbook of Mathematical
    # Functions. New York: Dover. Chapter 6: Gamma and Related Functions.
    y <- (abs(q - mean) / (sigma))^beta
    1 / 2 + sign(q - mean) * pgamma(y, 1 / beta) / 2
}

## *** Definition of Laplace CDF *** #

plaplace <- function(q, mean = 0, sigma = 1) {
    1 / sqrt(2) + 1 / sqrt(2) * (1 - exp(-abs(q - mean) / sigma))
}

## *** Generalized beta CDF *** #

pgbeta <-
    function(q,
    shape1,
    shape2,
    lambda = 1,
    lower.tail = TRUE,
    log.p = FALSE) {
        ans <- pbeta(
            q = 1 / (1 + (1 / q - 1) / lambda),
            shape1 = shape1,
            shape2 = shape2,
            lower.tail = lower.tail,
            log.p = log.p
        )
        ans[lambda <= 0] <- NaN
        return(ans)
    }

## *** 3P beta CDF *** #

pbeta3 <-
    function(q,
    shape1 = 2,
    shape2 = 3,
    mu = 0,
    lower.tail = TRUE,
    log.p = FALSE) {
        pbeta(
            q - mu,
            shape1,
            shape2,
            lower.tail = lower.tail,
            log.p = log.p
        )
    }

## *** 4P beta CDF *** #

pbeta4 <- function(q,
    shape1 = 2,
    shape2 = 3,
    mu = 0,
    b = 1,
    lower.tail = TRUE,
    log.p = FALSE) {
    pbeta((q - mu) / (b - mu),
        shape1,
        shape2,
        lower.tail = lower.tail,
        log.p = log.p
    )
}

## *** Definition of Rayleigh distribution *** #

prayleigh <- function(q, sigma) {
    1 - exp(-q^2 / (2 * sigma^2))
}

## *** Definition of 2P Exponential *** #
pexp2 <- function(q, rate, mu) {
    pexp(
        q - mu,
        rate = 1,
        lower.tail = TRUE,
        log.p = FALSE
    )
}

## *** Definition of Beta-Weibull distribution *** #
## It is taken from R package "Newdistns"

pbweibull <- function(q, alpha, beta, shape, scale) {
    pbeta(
        pweibull(
            q,
            shape = shape,
            scale = scale
        ),
        shape1 = alpha,
        shape2 = beta
    )
}

## *** Definition of Log-Gamma 3P distribution *** #
## Similar to the one defined in package VGAM
plgamma3p <- function(
    q,
    shape = 1,
    scale = 1,
    mu = 0,
    lower.tail = TRUE,
    log.p = FALSE) {

    y <- (q - mu)/scale
    res <- pgamma(exp(y), shape, lower.tail = lower.tail, log.p = log.p)
    res[scale < 0] <- NaN
    res
}

plgamma <- function(
        q,
        shape = 1,
        scale = 1,
        lower.tail = TRUE,
        log.p = FALSE) {

    y <- q/scale
    res <- pgamma(exp(y), shape, lower.tail = lower.tail, log.p = log.p)
    res[scale < 0] <- NaN
    res
}

## Initial parameters estimation by Method of Moments

beta_start_par <- function(MEAN, VAR) {
    shape1 <- MEAN * (MEAN * (1 - MEAN) / VAR - 1)
    shape2 <- shape1 * (1 - MEAN) / MEAN
    return(c(shape1 = shape1, shape2 = shape2))
}

## AIC
summary_aic <- function(m) {
    nms <- m[[1]]$aic[, 1]
    res <- do.call(cbind, lapply(m, function(x) {
        idx <- match(nms, x$aic$Distribution)
        x$aic[idx, 2]
    }))
    rownames(res) <- nms
    return(res)
}
