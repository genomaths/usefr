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

#' @rdname mcgoftest
#' @title Bootstrap test for Goodness of fit (GoF)
#' @description To accomplish the nonlinear fit of a probability distribution
#'     function (*PDF*), dIfferent optimization algorithms can be used. Each
#'     algorithm will return a different set of estimated parameter values. AIC
#'     and BIC are not useful (in this case) to decide which parameter set of
#'     values is the best. The goodness-of-fit tests (GOF) can help in this
#'     case. Please, see below the examples on how to use this function.
#' @details The test is intended for continuos distributions. If sampling size
#'     is lesser the size of the sample, then the test becomes a Monte Carlo
#'     test. The thes is based on the use of measures of goodness of fit,
#'     statistics. The following statistics are availible:
#'
#'     \itemize{
#'         \item Kolmogorov- Smirnov statistic (ks). Limitations: sensitive to
#'               ties [1]. Only the parametric Monte Carlo resampling
#'               (provided that there is not ties in the data) can be used.
#'
#'         \item Anderson–Darling statistic (ad) [2]. Limitation: by
#'               construction, it depends on the sample size. So, the size of
#'               the sampling must be close to the sample size if Monte Carlo
#'               resampling is used, which could be a limitation if the sample
#'               size is too large [2]. In particular, could be an issue in some
#'               genomic applications. It is worth highlighting that, for the
#'               current application, the Anderson–Darling statistic is not
#'               standardized as typically done in testing GoF for normal
#'               distribution with Anderson–Darling test. It is not required
#'               since, the statistic is not compared with a corresponding
#'               theoretical value. In addition, since the computation of this
#'               statistic requires for the data to be put in order [2], it does
#'               not make sense to perform a permutation test. That is, the
#'               maximum sampling size is the sample size less 1.
#'
#'         \item Pearson's Chi-squared statistic (chisq). Limitation: the sample
#'               must be discretized (partitioned into bins), which is could be
#'               a source of bias that leads to the rejection of the null
#'               hypothesis. Here, the discretization is done using function
#'               the resources from function \code{\link[graphics]{hist}}.
#'
#'         \item Root Mean Square statistic (rmse). Limitation: the same
#'               as 'chisq'.
#'     }
#' @param varobj A a vector containing observations, the variable for which the
#'     CDF parameters was estimated.
#' @param distr The name of the cummulative distribution function (CDF) or a
#'     concrete CDF from where estimate the cummulative probabilities.
#'     Distribution \emph{distr} must be defined in environment-namespace from
#'     any package or environment defined by user.
#' @param pars CDF model parameters. A list of parameters to evaluate the CDF.
#' @param stat One string denoting the statistic to used in the testing: "ks":
#'     Kolmogorov–Smirnov, "ad": Anderson–Darling statistic, "chisq: Pearson's
#'     Chi-squared, and "rmse": Root Mean Square of the error.
#' @param breaks Default is NULL. Basically, the it is same as in function
#'     \code{\link[graphics]{hist}}. If \emph{breaks} = NULL, then function
#'     'nclass.FD' (see \code{\link[grDevices]{nclass}} is applied to estimate
#'     the breaks.
#' @param parametric Logical object. If TRUE, then samples are drawn from the
#'     theoretical population described by \emph{distr}. Default: TRUE.
#' @param num.sampl Number of resamplings.
#' @param sample.size Size of the samples used for each sampling.
#' @param seed An integer used to set a 'seed' for random number generation.
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param verbose if verbose, comments and progress bar will be printed.
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @importFrom stats ks.test
#' @return A numeric vector with the following data:
#'     \enumerate{
#'
#'         \item Statistic value.
#'
#'         \item mc_p.value: the probability of finding the observed, or more
#'               extreme, results when the null hypothesis \eqn{H_0} of a study
#'               question is true obtained Monte Carlo resampling approach.
#'     }
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#' @references
#'     \enumerate{
#'         \item Feller, W. On the Kolmogorov-Smirnov Limit Theorems for
#'             Empirical Distributions. Ann. Math. Stat. 19, 177–189 (1948).
#'         \item Anderson, T. . & Darling, D. A. A Test Of Goodness Of Fit. J.
#'             Am. Stat. Assoc. 49, 765–769 (1954).
#'         \item Watson, G. S. On Chi-Square Goodness-Of-Fit Tests for
#'             Continuous Distributions. J. R. Stat. Soc. Ser. B Stat.
#'             Methodol. 20, 44–72 (1958).
#'     }
#' @seealso Distribution fitting: \code{\link{fitMixDist}},
#'     \code{\link[MASS]{fitdistr}}, \code{\link{fitCDF}}.
#' @examples
#' ## ======== Example 1 =======
#' # Let us generate a random sample a from a specified Weibull distribution:
#' # Set a seed
#' set.seed( 1 )
#' # Random sample from Weibull( x | shape = 0.5, scale = 1.2 )
#' x = rweibull(10000, shape = 0.5, scale = 1.2)
#'
#' # MC KS test accept the null hypothesis that variable x comes
#' # from Weibull(x | shape = 0.5, scale = 1.2), while the standard
#' # Kolmogorov-Smirnov test reject the Null Hypothesis.
#' mcgoftest(x, distr = "weibull", pars = c( 0.5, 1.2 ), num.sampl = 500,
#'         sample.size = 1000, num.cores = 4)
#'
#' ## ========= Example 2 ======
#' # Let us generate a random sample a random sample from a specified Normal
#' # distribution:
#' # Set a seed
#' set.seed( 1 )
#' x = rnorm(10000, mean = 1.5, sd = 2)
#'
#' # MC KS test accept the null hypothesis that variable x comes
#' # from N(x | mean = 0.5, sd = 1.2), while the standard
#' # Kolmogorov-Smirnov test reject the Null Hypothesis.
#' mcgoftest(x, distr = "norm", pars = c(1.5, 2), num.sampl = 500,
#'           sample.size = 1000, num.cores = 1)
#'
#' ## ========= Example 3 ======
#' ## Define a Weibull 3-parameter distribution function
#' pwdist <- function(x, pars) pweibull(x - pars[1], shape = pars[2],
#'                    scale = pars[3])
#' rwdist <- function(n, pars) rweibull(n, shape = pars[2],
#'                     scale = pars[3]) + pars[1]
#'
#' ## A random generation from Weibull-3P
#' set.seed(123)
#' pars <- c(mu = 0.9, shape = 1.4, scale = 3.7)
#' w <- rwdist(200, pars = pars)
#'
#' ## Testing GoF
#' mcgoftest(varobj = w, distr = "wdist", pars = list(pars), num.sampl = 100,
#'           sample.size =  199, stat = "chisq", num.cores = 4, breaks = 100,
#'           seed = 123)
#'
#' ## ========= Example 4 ======
#' ## ----- Testing GoF of a mixture distribution. ----
#' ## Define a mixure distribution to be avaluated with functions 'mixtdistr'
#' ## (see ?mixtdistr). In the current case, it will be mixture of a Log-Normal
#' ## and a Weibull distributions:
#'
#' phi = c( 0.37, 0.63) # Mixture proportions
#' args <- list(lnorm = c(meanlog = 0.837, sdlog = 0.385),
#'              weibull = c(shape = 2.7, scale = 5.8))
#' ##  Sampling from the specified mixture distribution
#' set.seed(123)
#' x <- rmixtdistr(n = 1e5, phi = phi , arg = args)
#' hist(x, 100, freq = FALSE)
#' x1 <- seq(0, 10, by = 0.001)
#' lines(x1, dmixtdistr(x1, phi = phi, arg = args), col = "red")
#'
#' ## The GoF for the simulated sample
#' pars <- c(list(phi = phi), arg = list(args))
#' mcgoftest(varobj = x, distr = "mixtdistr", pars = pars, num.sampl = 999,
#'         sample.size =  999, stat = "chisq", num.cores = 4, breaks = 200,
#'         seed = 123)
#'
mcgoftest <- function(varobj, distr, pars, num.sampl = 999, sample.size,
                   stat = c("ks", "ad", "rmse", "chisq"), breaks = NULL,
                   parametric = TRUE, seed = 1, num.cores = 1, tasks = 0,
                   verbose = TRUE) {
   ## A starting permutation script for permutation test used the idea
   ## published on \href{https://goo.gl/hfnNy8}{Alastair Sanderson}
   ## An idea presented in his post: \emph{Using R to analyse data statistical
   ## and numerical data analysis with R}. Herein, it is modified and extended.
   set.seed( seed )
   stat <- match.arg(stat)
   if (!is.character(distr))
       stop("*** 'distr' must be a characterstring naming a distribution ",
            "function e.g, 'norm' will permit accessing ",
            "functions 'pnorm' and 'rnorm'")

   if (is.character(distr)) {
       pdistr <- paste0("p", distr)
       pdistr <- get(pdistr, mode = "function", envir = parent.frame())
       if (!is.function(pdistr))
          stop("*** 'distr' must be a characterstring naming a distribution ",
               "function e.g, 'norm' will permit accessing ",
               "functions 'pnorm' and 'rnorm'")

       if (parametric) {
           rdistr <- paste0("r", distr)
           rdistr <- get(rdistr, mode = "function", envir = parent.frame())
           if (!is.function(rdistr))
              stop("*** 'distr' must be a characterstring naming a ",
                   "distribution function e.g, 'norm' will permit accessing ",
                   "functions 'pnorm' and 'rnorm'")
       }
   }


   if (missing( sample.size ) || sample.size > round(length(varobj)))
       sample.size = round(length(varobj)/3)

   GoFtest <- function(x, R, szise, distr, pars, stat, breaks, parametric,
                       mc.cores) {

       myfun <- function(a, distr, pars, stat, breaks, parametric) {
           switch(stat,
               ks = ks_stat(x = a, distr = distr, pars = pars)$stat,
               ad = ad_stat(x = a, distr = distr, pars = pars),
               rmse = rmse(x = a, distr = distr, pars = pars, breaks = breaks),
               chisq = chisq(x = a, distr = distr, pars = pars, breaks = breaks)
          )
       }

       DoIt <- function(r, distr, pars, stat, breaks, parametric) {
           if (parametric) a <- distfn(x = szise, dfn = distr,
                                       type = "r", arg = pars)
           else {
               i <- sample( length( x ), szise )
               a = x[ i ]
           }

           # to test empirical versus theoretical values
           myfun(a = a, distr = distr, pars = pars, stat = stat,
                   breaks = breaks, parametric = parametric)
       }

       # DoIt(1, distr=distr, pars=pars, stat=stat, breaks=breaks,
       #      parametric=parametric)

       if (verbose) progressbar = TRUE else progressbar = FALSE
       if (Sys.info()['sysname'] == "Linux") {
           bpparam <- MulticoreParam(workers=num.cores, tasks=tasks,
                                       progressbar = progressbar)
       } else bpparam <- SnowParam(workers=num.cores, type = "SOCK",
                                   progressbar = progressbar)

       pstats <- unlist(bplapply(1:R, DoIt, distr = distr, pars = pars,
                               stat = stat, breaks = breaks,
                               parametric = parametric, BPPARAM=bpparam))

       res <- switch(stat,
                   ks = {
                           statis = ks_stat(x = x, distr = distr, pars = pars)
                           p.value = statis$p.value
                           statis = statis$stat
                           c(KS.stat.D = statis,
                               mc_p.value = mean(c(statis, pstats) >= statis,
                                               na.rm = TRUE),
                               KS.stat.p.value = p.value,
                               sample.size = sample.size, num.sampl = num.sampl)
                        },
                   ad = {
                           statis = ad_stat(x = x, distr = distr, pars = pars)
                           c(AD.stat = statis,
                               mc_p.value=mean(c(statis, pstats) >= statis,
                                               na.rm = TRUE),
                               sample.size = sample.size, num.sampl = num.sampl)
                        },
                   rmse = {
                           statis = rmse(x = x, distr = distr, pars = pars ,
                                       breaks = breaks)
                           c(rmse = statis,
                             mc_p.value=mean(c(statis, pstats) >= statis,
                                             na.rm = TRUE),
                             sample.size = sample.size, num.sampl = num.sampl)
                          },
                   chisq = {
                               statis = chisq(x = x, distr = distr, pars=pars,
                                            breaks = breaks)
                               c(Chisq = statis,
                                   mc_p.value=mean(c(statis, pstats) >= statis,
                                                   na.rm = TRUE),
                                   sample.size = sample.size,
                                   num.sampl = num.sampl)
                           }
               )
       return(res)
   }

   statis <- switch(stat,
          ks = "Kolmogorov-Smirnov",
          ad = "Anderson–Darling",
          rmse = "Root Mean Square",
          chisq = "Pearson's Chi-squared"
   )

   if (verbose) {
       method <- ifelse(sample.size < length(varobj), "Monte Carlo GoF testing",
                       "Permutation GoF testing")

       if (sample.size == length(varobj)) parametric <- FALSE
       approach <- ifelse(parametric, "parametric approach",
                           "non-parametric approach")

       cat( "***", method,"based on" , statis,"statistic", "(", approach, ")",
           " ...\n")
   }

   if (stat == "ks" && !parametric && sample.size < length(varobj)) {
      warning("*** ", method," based on ", statis," statistic ",
              "(", approach, ") ", "is not reliable. \n",
              "Use the parametric approach or permutation instead.")
   }

   GoFtest(x = varobj, R = num.sampl, szise = sample.size, distr = distr,
           pars=pars, stat = stat, breaks = breaks, parametric = parametric,
           mc.cores = num.cores)
}

# ======================= Anderson–Darling statistic ========================= #
ad_stat <- function(x, distr, pars = NULL) {
   x <- sort(x[complete.cases(x)])
   if (!missing(distr))  x <- distfn(x = x, dfn = distr, type = "p", arg = pars)

   n <- length(x)
   if (n < 8)
      stop("AD statistic the sample size must be greater than 7")

   h <- x * (1 - rev(x))
   h <- (2 * seq(x) - 1) * log(h)
   return(-n - mean(h))
}

# ================ Auxiliary funciton to compute the frequencies ============= #
freqs <- function(x, distr, pars, breaks = NULL) {
   n <- length(x)
   if (is.null(breaks)) breaks <- nclass.FD(x)
   DENS <- hist(x, breaks = breaks, plot = FALSE)
   freq0 <- DENS$counts
   bounds <- data.frame(lower = DENS$breaks[ -length(DENS$breaks)],
                        upper = DENS$breaks[ -1 ])

   up_val <- distfn(x = bounds$upper, dfn = distr, type = "p", arg = pars)
   low_val <- distfn(x = bounds$lower, dfn = distr, type = "p", arg = pars)

   freq <- round(n * (up_val - low_val))
   return(data.frame(obsf = freq0, expf = freq))
}

# ===================== Root-Mean-Square statistic ========================== #
rmse <- function(x, distr, pars, breaks = NULL) {
   freq <- freqs(x = x, distr = distr, pars = pars, breaks = breaks)
   return(sqrt(sum((freq$obsf - freq$expf)^2, na.rm = TRUE))/length(freq$obsf))
}

# ================== Pearson's Chi-squared  statistic ======================= #
chisq <- function(x, distr, pars, breaks = NULL) {
   freq <- freqs(x = x, distr = distr, pars = pars, breaks = breaks)
   if (any(freq$expf == 0)) {
      freq$expf <- freq$expf + 1
      freq$obsf <- freq$obsf + 1
   }
   return(sum((freq$obsf - freq$expf)^2/freq$expf))
}

# ================== Kolmogorov-Smirnov  statistic ======================= #

ks_stat <- function(x, distr, pars) {
   cdf <- function(x, pars) distfn(x = x, dfn = distr, type = "p", arg = pars)
   res <- suppressWarnings(ks.test(x = x, y = cdf, pars))
   return(list(stat = unname(res$statistic), p.value = res$p.value))
}


# --------------------- Auxiliary function to get distribution --------------- #
distfn <- function(x, dfn, type = "r", arg) {
   switch(type,
          p = do.call(paste0(type, dfn),
                      c(list(x), arg)),
          r = do.call(paste0(type, dfn), c(list(x), arg))
   )
}
# -------------------------- End auxiliar function --------------------------- #


