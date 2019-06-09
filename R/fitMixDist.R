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

#' @rdname fitMixDist
#'
#' @title Nonlinear fit of Mixture distributions
#' @description This function performs the nonlinear fit of mixture
#'     distributions exploiting a firth approach on parameterized finite
#'     Gaussian mixture models obtained through the function
#'     \code{\link[mclust]{Mclust}} from package \emph{mclust}.
#' @details The approch tries to fit the proposed mixture distributions using a
#'     modification of Levenberg-Marquardt algorithm implemented in function
#'     \code{\link[minpack.lm]{nls.lm}} from \emph{minpack.lm} package that is
#'     used to perform the nonlinear fit. Cross-validations for the nonlinear
#'     regressions (R.Cross.val) are performed as described in reference [1]. In
#'     addition, Stein's formula for adjusted R squared (rho) was used as an
#'     estimator of the average cross-validation predictive power [1].
#' @param X numerical vector
#' @param arg A list of named vectors with the corresponding named distribution
#'     parameters values. The names of the vector of parameters and the
#'     parameter names must correspond to defined functions. For example, if
#'     one of the involved distributions is the gamma density (see
#'     \code{\link[stats]{GammaDist}}), then the corresponding vector of
#'     parameters must be gamma = c(shape = 'some value', scale = 'some value').
#'     For the following distributions, the arguments can be provided with NULL
#'     values:
#'     \itemize{
#'         \item "norm" \href{https://goo.gl/xaEAdT}{(Wikipedia)}
#'         \item "halfnorm" \href{https://goo.gl/yxMF6T}{(Wikipedia)}.
#'         \item "gnorm" \href{https://goo.gl/EPk8mH}{(Wikipedia)}
#'         \item "gamma" \href{https://goo.gl/cYkvar}{(Wikipedia)}
#'         \item "beta" \href{https://goo.gl/893wzR}{(Wikipedia)}
#'         \item "laplace" \href{https://goo.gl/fCykV9}{(Wikipedia)}
#'         \item "weibull" \href{https://goo.gl/WMXmQP}{(Wikipedia)}
#'         \item "rayleigh" \href{https://goo.gl/d9b3zv}{(Wikipedia)}
#'         \item "exp" \href{https://goo.gl/stVsi7}{(Wikipedia)}
#'     }
#'     Notice that the distribution given names correspond to the root-names as
#'     given for R functions. For example, 'gamma' is the root-name for
#'     functions \code{\link[stats]{GammaDist}}. See example, for more details.
#' @param npoints number of points used in the fit of the density function or
#'     NULL. These are used as histogram break points to estimate the empirical
#'     density values. If \emph{npoints} = NULL and \emph{dens} = TRUE, then.
#'     Kernel Density Estimation function \code{\link[stats]{density}} from
#'     \emph{stats} package is used to estimate the empirical density. Default
#'     value is 100.
#' @param maxiter positive integer. Termination occurs when the number of
#'     iterations reaches maxiter. Default value: 1024.
#' @param prior Same as in \code{\link[mclust]{Mclust}} function.
#' @param ftol non-negative numeric. Termination occurs when both the actual
#'     and predicted relative reductions in the sum of squares are at most ftol.
#'     Therefore, ftol measures the relative error desired in the sum of
#'     squares. Default value: 1e-12
#' @param ptol non-negative numeric. Termination occurs when the relative error
#'     between two consecutive iterates is at most ptol. Therefore, ptol
#'     measures the relative error desired in the approximate solution.
#'     Default value: 1e-12.
#' @param maxfev Integer; termination occurs when the number of calls to fn has
#'   reached maxfev. Note that nls.lm sets the value of maxfev to
#'   100*(length(par) + 1) if maxfev = integer(), where par is the list or
#'   vector of parameters to be optimized.
#' @param equalPro An argument to pass to \code{\link[mclust]{emControl}}
#'     function. Logical variable indicating whether or not the mixing
#'     proportions are equal in the model. Default: equalPro = FALSE.
#' @param eps,tol Arguments to pass to \code{\link[mclust]{emControl}} function.
#' @param usepoints Integer. Computation by function
#'     \code{\link[mclust]{Mclust}} could take long time when the sample size is
#'     about >= 10000. This number can be used to extract a random sample of
#'     size 'usepoints' and to do the estimation with it.
#' @param seed Seed for random number generation.
#' @param dens Logic. Whether to use fit the 'PDF' or 'CDF'. Default is TRUE.
#' @param kmean Logic. Whether to use \code{\link[stats]{kmeans}} algorithm to
#'     perform the estimation in place of \code{\link[mclust]{Mclust}}. Deafult
#'     is FALSE.
#' @param verbose if TRUE, prints the function log to stdout
#' @param ... Further arguments to pass to other functions like
#'     \code{\link[mclust]{Mclust}} and \code{\link[stats]{density}}.
#' @return A list with the model table with coefficients and goodness-of-fit
#'     results, the fitted model returned by function
#'     \code{\link[minpack.lm]{nls.lm}}, and a named list of fitted arguments.
#'
#' @references 1. Stevens JP. Applied Multivariate Statistics for the Social
#'     Sciences. Fifth Edit. Routledge Academic; 2009.
#' @author Robersy Sanchez - 05/13/2019
#'
#' @importFrom stats var stepfun as.formula coef AIC pweibull BIC vcov knots
#' @importFrom stats deviance BIC cor quantile
#' @importFrom nls2 nls2
#' @importFrom utils head
#' @importFrom stats ecdf nls.control sd
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @importFrom mclust Mclust mclustBIC priorControl emControl
#' @importFrom mixdist weibullpar
#' @seealso \code{\link[MASS]{fitdistr}} and \code{\link{fitCDF}}
#' @examples
#' set.seed(123) # set a seed for random generation
#' ## ========= A mixture of three distributions =========
#' phi = c(3/10, 7/10) #' Mixture proportions
#' ## ---------------------------------------------------------
#'
#' ## === Named vector of the corresponding distribution function parameters
#' ## must be provided
#' args <- list(gamma = c(shape = 2, scale = 0.1), weibull = c(shape = 3, scale = 0.5))
#' ## ------------------------------------------------------------
#'
#' ## ===== Sampling from the specified mixture distribution ====
#' X <- rmixtdistr(n = 1e5, phi = phi , arg = args)
#' ## ------------------------------------------------------------
#'
#' ## ===== Nonlinear fit of the specified mixture distribution ====
#' FIT <- fitMixDist(X, args = list(gamma = c(shape = NA, scale = NA),
#'                                  weibull = c(shape = NA, scale = NA)),
#'                   npoints = 200, usepoints = 1000)
#'
#' ## === The graphics for the simulated dataset and the corresponding theoretical
#' ## mixture distribution
#' par(bg = "gray98",  mar = c(3, 4, 2, 1) )
#' hist(X, 90, freq = FALSE, las = 1, family = "serif", col = rgb(0, 0, 1, 0.),
#'      border = "deepskyblue", main = "Histogram of Mixture Distribution")
#' x1 <- seq(-4, 10, by = 0.001)
#' lines(x1, dmixtdistr(x1, phi = phi, arg = args), col = "red")
#' lines(x1, dmixtdistr(x1, phi = FIT$phi, arg = FIT$args), col = "blue")
#' legend(1, 1.5, legend=c("Theoretical Mixture PDF", "Estimated Mixture PDF"),
#'        col=c("red", "blue"), lty=1, cex=0.8)
#' @export
fitMixDist <- function(X, args = list(norm = c(mean = NA, sd = NA),
                                       weibull = c(shape = NA, scale = NA)),
                       npoints = 100, maxiter=1024, prior = priorControl(),
                       ftol=1e-14, ptol=1e-14, maxfev = 1e+5, equalPro = FALSE,
                       eps, tol, usepoints, seed = 123, dens = TRUE,
                       kmean = FALSE, verbose=TRUE, ...) {

   dfns <- names(args)
   par = unname(unlist(args))
   numpar <- unlist(lapply(args, length))
   nampar <- lapply(args, names)

   if (any(is.element(dfns, c("weibull", "gamma", "lnorm", "hnorm","beta"))))
       X <- X[X > 0]

   ## ========== Auxiliary functions & Starting parameter values ============= #

   weibullpars <- function(mu, sigma) {
     return(weibullpar(mu = mu, sigma = sigma)[-3])
   }

   gammapars <- function(mu, sigma) {
     return(c(shape = (mu/sigma)^2, scale = sigma^2/mu))
   }

   betapars <- function(mu, var) {
     alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
     return(c(shape1 = alpha, shape2 = alpha * (1 / mu - 1)))
   }

   distNames <- c("norm", "hnorm", "gnorm", "gamma", "beta", "laplace",
                  "weibull", "rayleigh", "exp")

   par_vector2list <- function(par) {
       i = 1
       arg <- list()
       for(j in 1:length(dfns)) {
           parm <- c()
           for(k in 1:numpar[j]) {
               parm <- c(parm, par[i])
               i <- i + 1
           }
           names(parm) <- unlist(nampar[j])
           arg[[j]] <- parm
       }
       names(arg) <- dfns
       return(arg)
   }

   parLIST <- function(dfn, MEAN = NULL, VAR = NULL, SD = NULL) {
     return(switch(dfn,
                   norm = c(mean = MEAN, sd = SD),
                   hnorm = c(theta = 1/MEAN),
                   gnorm = c( mean = MEAN, sigma = SD, beta = 2 ),
                   laplace = c( mean = MEAN, sigma = SD ),
                   gamma = gammapars(mu = MEAN, sigma = SD),
                   weibull = weibullpars(mu = MEAN, sigma = SD),
                   beta = betapars(mu = MEAN, var = VAR),
                   rayleigh = c( sigma = SD ),
                   exp = c( rate = 1 )
           )
       )
   }

   # ==================== Starting parameter values ==================== #
   N <- length(X)
   if (!is.null(npoints)) {
       DENS <- hist(X, breaks = npoints, plot = FALSE)
       x <- DENS$mids
       n <- length(x)
   } else n <- N

   if (dens & is.null(npoints)) {
       DENS <- density(X, ...)
       x <- DENS$x
       y <- DENS$y
   }

   if (dens & !is.null(npoints)) y <- DENS$density

   if (!dens) {
       Fn <- ecdf(X)
       if (!is.null(npoints)) {
           y <- Fn(x)
       } else {
           y <- Fn(X)
           x <- X
       }
   }

   if (kmean && any(is.na(unlist(args))) && all(is.element(dfns, distNames)))
   {
       cl <- mixture_stat(u = X, m = length(dfns))
       mu <- cl$mu
       sigma <- cl$sigma
       phi <- cl$prop
       args <- mapply(parLIST, dfns, MEAN = mu, SD = sigma, SIMPLIFY = FALSE)
   }

   if (any(is.na(unlist(args))) && all(is.element(dfns, distNames)) &&
       !kmean) {
       if (missing(eps)) eps = .Machine$double.eps
       if (missing(tol)) tol = c(1e-5,sqrt(.Machine$double.eps))
       if (!missing(usepoints)) {
          set.seed(seed)
          Z <- sample(X, usepoints)
       } else Z <- X
       fit <- Mclust(Z, G = length(dfns), model="V", prior = prior,
                     control = emControl(eps = eps, tol =  tol,
                                         equalPro = equalPro))
       rm(Z); gc()
       mu <- fit$parameters$mean
       sigma <- sqrt(fit$parameters$variance$sigmasq)
       phi <- fit$parameters$pro

       args <- mapply(parLIST, dfns, MEAN = mu, SD = sigma, SIMPLIFY = FALSE)
   } else {
       if (any(is.na(unlist(args))))
           stop("*** The list of functions and starting arguments",
               " must be given")
       if (!kmean) {
           if (missing(eps)) eps = .Machine$double.eps
           if (missing(tol)) tol = c(1e-5,sqrt(.Machine$double.eps))
           if (!missing(usepoints)) {
               set.seed(seed)
               Z <- sample(X, usepoints)
           } else Z <- X
           fit <- Mclust(Z, G = length(dfns), model="V", prior = prior,
                        control = emControl(eps = eps, tol =  tol,
                                            equalPro = equalPro))
           phi <- fit$parameters$pro

       } else {
           cl <- mixture_stat(u = X, m = length(dfns))
           phi <- cl$prop
       }
   }

   ## -------------------- END starting parameter values --------------------- #


   ## =========================== Fitting models ============================= #
   optFun <- function(par, objFun, quantiles, obsVals, eval = FALSE) {
      if (dens) {
         START <- list(phi = phi, arg = par_vector2list(par),
                       x = quantiles)
      } else {
         START <- list(phi = phi, arg = par_vector2list(par),
                       q = quantiles)
      }
      EVAL <- try(do.call(objFun, START), silent = TRUE)
      if (inherits(EVAL, "try-error")) return(NA)
      EVAL[is.nan(EVAL)] <- 0
      RSS <- (obsVals - EVAL)
      if (eval) {
         return(EVAL)
      } else return(RSS)
   }

   if (dens) {
      FIT <- try(nls.lm(par = unlist(args), fn = optFun, objFun = dmixtdistr,
                        quantiles = x, obsVals = y,
                        control = nls.lm.control(maxiter = maxiter, ftol = ftol,
                                                 maxfev = maxfev, ptol = ptol)),
                 silent = TRUE)
   } else {
      # optFun(par = unlist(args), objFun = pmixtdistr, quantiles = x,
      #        obsVals = y)

      FIT <- try(nls.lm(par = unlist(args), fn = optFun, objFun = pmixtdistr,
                        quantiles = x, obsVals = y,
                        control = nls.lm.control(maxiter = maxiter, ftol = ftol,
                                                 maxfev = maxfev, ptol = ptol)),
                   silent = TRUE)
   }

   if (inherits( FIT, "try-error" ) && !kmean) {
       cl <- mixture_stat(u = X, m = length(dfns))
       mu <- cl$mu
       sigma <- cl$sigma
       phi <- cl$prop
       args <- mapply(parLIST, dfns, MEAN = mu, SD = sigma, SIMPLIFY = FALSE)

       if (dens) {
           FIT <- try(nls.lm(par = unlist(args), fn = optFun,
                       objFun = dmixtdistr, quantiles = x, obsVals = y,
                       control = nls.lm.control(maxiter = maxiter, ftol = ftol,
                                               maxfev = maxfev, ptol = ptol)),
                   silent = TRUE)
       } else {
           FIT <- try(nls.lm(par = unlist(args), fn = optFun,
                       objFun = pmixtdistr, quantiles = x, obsVals = y,
                       control = nls.lm.control(maxiter = maxiter, ftol = ftol,
                                               maxfev = maxfev, ptol = ptol)),
                   silent = TRUE)
       }
   }

   if (!inherits( FIT, "try-error" )) {
       ## **** R squares ****
       Adj.R.Square <- (1 - (deviance(FIT) / ((N - length(coef(FIT))) *
                                             var(y, use="everything"))))
       Adj.R.Square <- ifelse(is.na(Adj.R.Square) || Adj.R.Square < 0,
                           0, Adj.R.Square)

       ## Stein adjusted R square
       rho = (1 - ((N - 2) / (N - 3)) * ((N + 1) / (N)) * (1 - Adj.R.Square))
       rho = ifelse( is.na( rho ) | rho < 0, 0, rho )

       ## =========================== Cross Validation ======================= #
       ##--- Crossvalidation standard model for Nonlinear regression: x versus r

       if (verbose) {
           cat(paste("*** Performing nonlinear regression model ",
                "crossvalidation...\n" ))
       }
       arg1 <- unlist(coef(FIT))

       if (dens && is.null(npoints)) {
          set.seed(seed)
          l <- length(x)
          cros.ind.1 <- sample.int(l, size=round(l / 2))
          cros.ind.2 <- setdiff(1:l, cros.ind.1)
          x1 <- x[ cros.ind.1 ]
          x2 <- x[ cros.ind.2 ]
          y1 <- y[ cros.ind.1 ]
          y2 <- y[ cros.ind.2 ]
       } else {
          if (dens && !is.null(npoints)) {
             set.seed(seed)
             cros.ind.1 <- sample.int(N, size=round(N / 2))
             cros.ind.2 <- setdiff(1:N, cros.ind.1)
             if (npoints > length(cros.ind.1)) {
                breaks <- round(0.6 * length(cros.ind.1))
             } else breaks = npoints

             DENS <- hist(X[ cros.ind.1 ], breaks = breaks, plot = FALSE)
             x1 <- DENS$breaks[-1]
             y1 <- DENS$density; rm(DENS)

             DENS <- hist(X[ cros.ind.2 ], breaks = breaks, plot = FALSE)
             x2 <- DENS$breaks[-1]
             y2 <- DENS$density; rm(DENS); gc()
          }
       }

       if (dens) {
           FIT1 <- try(nls.lm(par = arg1, fn = optFun, objFun = dmixtdistr,
                            quantiles = x1, obsVals = y1,
                            control = nls.lm.control(maxiter = maxiter,
                                                   ftol = ftol, maxfev = maxfev,
                                                   ptol = ptol)),
                       silent = TRUE)

           FIT2 <- try(nls.lm(par = arg1, fn = optFun, objFun = dmixtdistr,
                           quantiles = x2, obsVals = y2,
                           control = nls.lm.control(maxiter = maxiter,
                                                   ftol = ftol, maxfev = maxfev,
                                                   ptol = ptol)),
                       silent = TRUE)
       } else {
           set.seed(seed)
           l <- length(x)
           cros.ind.1 <- sample.int(l, size=round(l / 2))
           cros.ind.2 <- setdiff(1:l, cros.ind.1)

           x1 <- X[ cros.ind.1 ]
           x2 <- X[ cros.ind.2 ]
           y1 <- Fn(x1)
           y2 <- Fn(x2)

           FIT1 <- try(nls.lm(par = arg1, fn = optFun, objFun = pmixtdistr,
                           quantiles = x1, obsVals = y1,
                           control = nls.lm.control(maxiter = maxiter,
                                                   ftol = ftol, maxfev = maxfev,
                                                   ptol = ptol)),
                       silent = TRUE)

           FIT2 <- try(nls.lm(par = arg1, fn = optFun, objFun = pmixtdistr,
                           quantiles = x2, obsVals = y2,
                           control = nls.lm.control(maxiter = maxiter,
                                                   ftol = ftol, maxfev = maxfev,
                                                   ptol = ptol)),
                       silent = TRUE)
       }

       if (inherits(FIT1, "try-error") && inherits(FIT2, "try-error"))
           R.cross.FIT <- 0
       else {
           if (dens) getPred <- dmixtdistr else getPred <- pmixtdistr
           ## prediction using model 1
           p.FIT1 <- getPred(x2, phi=phi, arg = par_vector2list(coef(FIT1)))
           R.FIT1 <- try(cor(p.FIT1, y2, use = "complete.obs"), silent = TRUE)
           if (inherits(R.FIT1, "try-error")) {
              R.FIT1 <- try(cor(p.FIT1, y2, use = "pairwise.complete.obs"),
                            silent = TRUE)
           }
           ## prediction using model 2
           p.FIT2 <- getPred(x1, phi = phi, arg = par_vector2list(coef(FIT2)))
           R.FIT2 <- try(cor(p.FIT2, y1, use = "complete.obs"), silent = TRUE)
           if (inherits(R.FIT2, "try-error")) {
              R.FIT2 <- try(cor(p.FIT2, y1, use = "pairwise.complete.obs"),
                            silent = TRUE)
           }

           if (inherits(R.FIT1, "try-error") && inherits(R.FIT2, "try-error")) {
               R.cross.FIT <- NA
           } else {
               term1 <- length(p.FIT1) * R.FIT1
               term2 <- length(p.FIT2) * R.FIT2
               R.cross.FIT <- (term1 + term2)/(length(p.FIT1) + length(p.FIT2))
           }
       }
       gaps <- rep("", length(coef(FIT)) - 1)
       sumario <- try(summary(FIT)$coefficients, silent = TRUE)
       if (!kmean) bic <- fit$bic else bic <- NA

       if (inherits(sumario, "try-error")) {
           stats <- data.frame(Estimate = c("Choleski Decomposition fail", NA),
                               Adj.R.Square=c(Adj.R.Square, ""),
                               rho=c(rho, ""),
                               R.Cross.val=c(R.cross.FIT, ""),
                               DEV=c(deviance(FIT), ""),
                               BIC=c(bic, ""),
                               n=c(N , n))
       } else {
           stats <- data.frame(sumario,
                               Adj.R.Square=c(Adj.R.Square, gaps),
                               rho=c(rho, gaps),
                               R.Cross.val=c(R.cross.FIT, gaps),
                               DEV=c(deviance(FIT), gaps),
                               BIC=c(bic, gaps),
                               n=c(N , n, gaps[-1] ))
       }



   } else {
       warning(paste("Data did not fit to the model.",
                       "Returning result from Mclust"))
       stats <- data.frame(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                           NA, NA, NA)
   }

   if (!inherits(sumario, "try-error")) {
      colnames(stats) <- c( "Estimate", "Std. Error", "t value", "Pr(>|t|))",
                            "Adj.R.Square", "rho", "R.Cross.val", "DEV",
                            "BIC", "size")
   }

   return(list(stats = stats, fit = FIT, args = par_vector2list(coef(FIT)),
               phi = phi))
}

######################################################################
#
# =========Suitable starting values for mixture distributions ====== #
#
######################################################################

# This function uses a k-means algorithm to heuristically select suitable
# starting values for mixture distributions

mixture_stat <- function( u, m,...) {
   km  <- kmeans( u, centers = m, ...)
   prop <- km$size/sum( km$size )  #  Estimate proportion
   mu  <- sapply( 1:m, function(i) km$centers[i, ] )
   # the covariance matrix for each group
   if( class( u ) == "numeric" ) {
      get.var <- function(i) sd( u[ km$cluster == i ] )
      sigma <- unlist( lapply( 1:m, get.var ) )
      idx = order( mu )
      sigma = sigma[ idx ]
      prop = prop[ idx ]
      mu = mu[ idx ]
   }else{
      get.var <- function(i) var( u[ km$cluster == i, ] )
      sigma <- lapply( 1:m, get.var )
   }
   list( num.cluster = m, dim = ncol( u ), prop = prop, mu = mu, sigma = sigma )
}

