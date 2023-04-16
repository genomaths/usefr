## ########################################################################## #
##
## Copyright (C) 2019-2023 Robersy Sanchez <https://genomaths.com/>
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
## ########################################################################## #
#
#' @rdname fitCDF2
#' @title Nonlinear fit of a commutative distribution function
#' @description Usually the parameter estimation of a cumulative distribution
#' function (*CDF*) are accomplished using the corresponding probability
#' density function (*PDF*). Different optimization algorithms can be used to
#' accomplished this task and different algorithms can yield different
#' estimated parameters. Hence, why not try to fit the CDF directly?
#' @param varobj A a vector, a named list, a matrix or a data.frame, containing
#' the observations from the variable for which the CDF parameters will be
#' estimated. When the argument is a matrix or a data.frame, the columns must
#' be named, carrying the objective variables.
#' @param distNames a vector of distribution numbers to select from the listed
#' below in details section, e.g. c(1:10, 15). If 'distNames' is not any of
#' current 20 named distributions, then it can be any arbitrary character
#' string, but the argument 'distf' must be given (see below).
#' @param start A named numerical vector giving the parameters to be optimized
#' with initial values or a list of numerical vectors (only when
#' \emph{\strong{varobj}} is a list, a matrix or a data.frame). This can be
#' omitted for some of the named distributions (see Details). This argument
#' will be used if provided for only one distribution. The default parameter
#' values are:
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
#' \strong{norm}, etc, from where, internally, we can get: density,
#' distribution function, quantile function and random generation as:
#' \strong{dnorm}, \strong{pnorm}, \strong{qnorm}, and \strong{rnorm},
#' respectively. If the function is not present in the environment, then an
#' error will be returned. It must given only if 'distNames' is not any of
#' current 20 named distributions (see  details below). Default is NULL.
#' @param only.info Logic. Default TRUE. If true, only information about the
#'   parameter estimation is returned.
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
#' @param num.cores The number of cores to use, i.e. at most how many child
#' processes will be run simultaneously. This argument will be passed to
#' \code{\link[parallel]{clusterApply}}.
#' @param ... (Optional) Further graphical parameters (see
#' \code{\link[graphics]{par}}). Graphical parameter will simultaneously affect
#' all the plots.
#' @param verbose Logic. If TRUE, prints the function log to stdout
#' @details This function works as function \code{\link{fitCDF}}, except for
#' the parallel computation, which is done applying function
#' \code{\link[parallel]{clusterApply}} from the 'parallel' R package.
#'
#' The nonlinear fit (NLF) problem for CDFs is addressed with
#' Levenberg-Marquardt algorithm implemented in function
#' \code{\link[minpack.lm]{nls.lm}} from package *minpack.lm*. The Stein's rho
#' for adjusted R squared (rho) is applied as an estimator of the average
#' cross-validation predictive power [1]. This function is inspired in a script
#' for the function \code{\link[propagate]{fitDistr}} from the package
#' propagate [2]. Some parts or script ideas from function
#' \code{\link[propagate]{fitDistr}} are used, but here we to estimate CDF and
#' not the PDF as in the case of "\code{\link[propagate]{fitDistr}}. More
#' informative results are given now. The studentized residuals are provided as
#' well. The list (so far) of possible CDFs is:
#'     \enumerate{
#'         \item Normal \href{https://goo.gl/xaEAdT}{(Wikipedia)}
#'         \item Log-normal \href{https://goo.gl/a7MtYq}{(Wikipedia)}. This
#'         This function is set to fit \eqn{log(1+x)}. Users can transform
#'         their variable by themself and then try the fitting to Normal
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
#'     After cdf = fitCDF2( varobj, ...), attributes( cdf$bestfit ) shows the
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
#' @import stats
#' @importFrom graphics par grid lines mtext abline
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @seealso \code{\link{fitCDF}}, \code{\link[MASS]{fitdistr}},
#' \code{\link{fitMixDist}} and for goodness-of-fit: \code{\link{mcgoftest}}.
#' @references
#'   \enumerate{
#'     \item Stevens JP. Applied Multivariate Statistics for the Social
#'           Sciences. Fifth Edit. Routledge Academic; 2009.
#'     \item Andrej-Nikolai Spiess (2014). propagate: Propagation of
#'           Uncertainty. R package version 1.0-4.
#'           http://CRAN.R-project.org/package=propagate
#'     \item Abramowitz, M. and Stegun, I. A. (1972) Handbook of Mathematical
#'           Functions. New York: Dover. Chapter 6: Gamma and Related
#'           Functions.
#'     \item Hand-book on STATISTICAL DISTRIBUTIONS for experimentalists
#'           (pag 73) by Christian Walck. Particle Physics Group Fysikum.
#'           University of Stockholm (e-mail: walck@physto.se).
#'   }
#' @examples
#' set.seed(1230)
#' x1 <- rnorm(10000, mean = 0.5, sd = 1)
#' cdfp <- fitCDF2(x1, distNames = "Normal", plot = FALSE)
#' summary(cdfp$bestfit)
#'
#' ## Add some cosmetics to the plots
#' cdfp <- fitCDF2(x1,
#'     distNames = "Normal", xlabel = "My Nice Variable Label",
#'     plot = T, font.lab = 3, font = 2, font.axis = 2, family = "serif",
#'     cex.lab = 1.3, cex.axis = 1.3
#' )
#'
#' ## Fitting a Weibull distribution with 3 paramaters
#' x1 <- rweibull3p(1000, shape = 0.5, scale = 1, mu = 0.1)
#' cdfp <- fitCDF2(x1,
#'     distNames = "3P Weibull",
#'     xlabel = "My Nice Variable Label",
#'     plot = T, font.lab = 3, font = 2, font.axis = 2, family = "serif",
#'     cex.lab = 1.3, cex.axis = 1.3, cex.main = 1.1,
#'     mgp = c(2.5, 1, 0)
#' )
#' @aliases fitCDF2
setGeneric(
    "fitCDF2",
    def = function(varobj,
    ...) {
        standardGeneric("fitCDF2")
    }
)

#' @rdname fitCDF2
#' @aliases fitCDF2
#' @export
setMethod(
    "fitCDF2", signature(varobj = "numeric"),
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
        res <- fitCDF(
            varobj = varobj,
            distNames = distNames,
            plot = plot,
            plot.num = plot.num,
            distf = distf,
            start = start,
            loss.fun = loss.fun,
            min.val = min.val,
            only.info = only.info,
            maxiter = maxiter,
            maxfev = maxfev,
            ptol = ptol,
            nls.model = nls.model,
            algorithm = algorithm,
            xlabel = xlabel,
            mar = mar,
            mgp = mgp,
            las = las,
            cex.main = cex.main,
            cex.text = cex.text,
            cex.point = cex.point,
            verbose = verbose,
            ...
        )
        return(res)
    }
)

setClassUnion(
    "list_OR_matrix_OR_dataframe",
    c("list", "matrix", "data.frame")
)

#' @rdname fitCDF2
#' @aliases fitCDF2
#' @importFrom parallel detectCores makeCluster stopCluster clusterApply
#' @export
setMethod(
    "fitCDF2", signature(varobj = "list_OR_matrix_OR_dataframe"),
    function(
        varobj,
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
        num.cores = 1L,
        verbose = TRUE,
        ...) {

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

        if (inherits(varobj, c("matrix", "data.frame"))) {
            varobj <- as.matrix(varobj)
            varobj <- split(varobj,
                            rep(1:ncol(varobj), each = nrow(varobj)))
        }

        ## ------------ Setting up parallel computation ------------ #

        vlen <- length(varobj)
        nc <- detectCores()
        nc <- min(vlen, nc)

        if (num.cores > nc)
            num.cores <- nc

        if (Sys.info()["sysname"] == "Linux")
            cl <- makeCluster(num.cores, type = "FORK")
        else
            cl <- makeCluster(num.cores, type = "SOCK")

        # registerDoParallel(cl)

        ## -------------------------------------------------------- #


        res <- clusterApply(cl, seq_along(varobj), function(k) {
                    fitCDF(
                    varobj = varobj[[k]],
                    distNames = distNames[[k]],
                    plot = plot,
                    plot.num = plot.num,
                    distf = distf[[k]],
                    start = start[[k]],
                    loss.fun = loss.fun[k],
                    min.val = min.val,
                    only.info = only.info,
                    maxiter = maxiter,
                    maxfev = maxfev,
                    ptol = ptol,
                    nls.model = nls.model,
                    algorithm = algorithm,
                    xlabel = xlabel,
                    mar = mar,
                    mgp = mgp,
                    las = las,
                    cex.main = cex.main,
                    cex.text = cex.text,
                    cex.point = cex.point,
                    verbose = FALSE)})

        stopCluster(cl)

        names(res) <- var_nms
        res$AICs <- summary_aic(res)
        res <- structure(res, class = "CDFmodelList")
        return(res)
    }
)

