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
#' function \emph{(PDF)}, different optimization algorithms can be used. Each
#' algorithm will return a different set of estimated parameter values. AIC and
#' BIC are not useful (in this case) to decide which parameter set of values is
#' the best. The goodness-of-fit tests \emph{(GOF)} can help in this case.
#' Please, see below the examples on how to use this function.
#'
#' @details The test is intended for mostly continuous distributions. Basically,
#' given the set of parameter values \strong{\emph{pars}} from distribution
#' \strong{\emph{distr}}, \strong{\emph{num.sampl}} sets of random samples will
#' be generated, each one of them with \strong{\emph{sample.size}} element. The
#' selected statistic \strong{\emph{pars}} will be computed for each randomly
#' generated set (\emph{b_stats}) and for sample \strong{\emph{varobj}}
#' (\emph{s_stat}). Next, the bootstrap \emph{p-value} will be computed as:
#' \eqn{mean(c(s_stat, b_stats) >= s_stat)}.
#'
#' If both variables, \strong{\emph{varobj}} and \strong{\emph{distr}}, are
#' numerical vectors, then \code{\link{tableBoots}} function is applied. That
#' is, the problem is confronted as a \eqn{n x m} contingency independence test,
#' since there is no way to prove that two arbitrary sequences of integer
#' numbers would follow the same probability distribution without provide
#' additional information/knowledge.
#'
#' If sampling size is lesser the size of the sample, then the test becomes a
#' Monte Carlo test. The test is based on the use of measures of goodness of
#' fit, statistics. The following statistics are available (and some
#' limitations for their application to continuous variables are given):
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
#'         \item Hellinger Divergence statistic (hd). Limitation: the same
#'               as 'chisq'.
#'        }
#' If the argument \strong{\emph{distr}} must be defined in
#' environment-namespace from any package or the environment defined by the
#' user. if  \eqn{missing( sample.size )}, then
#' \eqn{sample.size <- length(varobj) - 1}.
#'
#' Notice that 'chisq', 'rmse', and 'hd' tests can be applied to testing two
#' discrete probability distributions as well. However, here,
#' \strong{\emph{mcgoftest}} function is limited to continuous probability
#' distributions.
#'
#' Additionally, the only supported n-dimensional probability distribution is
#' Dirichlet Distribution (\emph{Dir}). The GOF for Dir is based on the fact
#' that if a variable \eqn{x = (x_1, x_2, ...x_n)} follows Dirichlet
#' Distribution with parameters \eqn{\alpha = \alpha_1, ... , \alpha_n} (all
#' positive reals), in short, \eqn{x ~ Dir(\alpha)}, then \eqn{x_i ~
#' Beta(\alpha_i, \alpha_0 - \alpha_i)}, where Beta(.) stands for the Beta
#' distribution and \eqn{\alpha_0 = \sum \alpha_i} (see Detail section,
#' \code{\link{dirichlet}} function, and example 5).
#'
#' @param varobj A a vector containing observations, the variable for which the
#' CDF parameters was estimated or the discrete absolute frequencies of each
#' observation category.
#' @param distr The possible options are:
#' \enumerate{
#'    \item The name of the cumulative distribution function (CDF).
#'    \item A concrete CDF, defined by the user,
#'          from where to estimate the cumulative probabilities.
#'    \item A numerical vector carrying discrete absolute frequencies for the
#'          same categories provided in \strong{\emph{varobj}}.
#' }
#'
#' @param pars CDF model parameters. A list of parameters to evaluate the CDF.
#' @param stat One string denoting the statistic to used in the testing:
#' \enumerate{
#'    \item "ks": Kolmogorov–Smirnov.
#'    \item "ad": Anderson–Darling statistic.
#'    \item "chisq: Pearson's Chi-squared.
#'    \item "rmse": Root Mean Square of the error.
#'    \item "hd": Hellinger Divergence statistics.
#' }
#' @param breaks Default is NULL. Basically, the it is same as in function
#' \code{\link[graphics]{hist}}. If \emph{breaks} = NULL, then function
#' 'nclass.FD' (see \code{\link[grDevices]{nclass}} is applied to estimate
#' the breaks.
#' @param par.names (Optional) The names of the parameters from
#' \strong{\emph{distr}} function. Some distribution functions would require to
#' provide the names of their parameters.
#' @param num.sampl Number of resamplings.
#' @param sample.size Size of the samples used for each sampling.
#' @param seed An integer used to set a 'seed' for random number generation.
#' @param num.cores,tasks Parameters for parallel computation using package
#' \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to use,
#' i.e. at most how many child processes will be run simultaneously (see
#' \code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
#' for Linux OS).
#' @param verbose if verbose, comments and progress bar will be printed.
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @importFrom stats ks.test rmultinom
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
#'         \item A. Basu, A. Mandal, L. Pardo, Hypothesis testing for two
#'             discrete populations based on the Hellinger distance. Stat.
#'             Probab. Lett. 80, 206–214 (2010).
#'     }
#' @seealso Distribution fitting: \code{\link{fitMixDist}},
#'     \code{\link[MASS]{fitdistr}}, \code{\link{fitCDF}}, and
#'     \code{\link{bicopulaGOF}}.
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
#' ## Define a mixture distribution to be evaluated with functions 'mixtdistr'
#' ## (see ?mixtdistr). In the current case, it will be mixture of a Log-Normal
#' ## and a Weibull distributions:
#'
#' phi = c( 0.37, 0.63) # Mixture proportions
#' args <- list(lnorm = c(meanlog = 0.837, sdlog = 0.385),
#'              weibull = c(shape = 2.7, scale = 5.8))
#'
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
#' ## ========= Example 5 ======
#' ## GoF for Dirichlet Distribution (these examples can be run,
#' ## the 'not-run' is only to prevent time consuming in R checking)
#' \dontrun{
#'         set.seed(1)
#'         alpha = c(2.1, 3.2, 3.3)
#'         x <- rdirichlet(n = 100, alpha = alpha)
#'
#'         mcgoftest(varobj = x, distr = "dirichlet",
#'                   pars = alpha, num.sampl = 999,
#'                   sample.size =  100, stat = "chisq",
#'                   par.names = "alpha",
#'                   num.cores = 4, breaks = 50,
#'                   seed = 1)
#'
#'         ## Now, adding some noise to the sample
#'         set.seed(1)
#'         x <- x + replicate(3, runif(100, max = 0.1))
#'         x <- x/ rowSums(x)
#'
#'         mcgoftest(varobj = x, distr = "dirichlet",
#'                   pars = alpha, num.sampl = 999,
#'                   sample.size =  100, stat = "chisq",
#'                   par.names = "alpha",
#'                   num.cores = 4, breaks = 50,
#'                   seed = 1)
#' }
#'
#' ## ========= Example 6 ======
#' ##  Testing whether two discrete probabilities vectors
#' ##  would come from different probability distributions.
#' set.seed(1)
#' ## Vector of absolute frequencies
#' n = 12
#' alpha = round(runif(n, min = 0.1, max = 0.8) * 10, 1)
#' x1 <- round(runif(n = n) * 100)
#'
#' ## Add some little noise
#' x2 <- x1 + round(runif(n, min = 0.1, max = 0.2) * 10)
#'
#' mcgoftest(varobj = x1, distr = x2,
#'           num.sampl = 999,
#'           stat = "chisq",
#'           seed = 1)
#'
#' ## Compare it with Chi-square test from R package 'stat'
#' chisq.test(rbind(x1,x2), simulate.p.value = TRUE, B = 2e3)$p.value
#'
#' ## Add bigger noise to 'x1'
#' x3 <- x1 + round(rnorm(n, mean = 5, sd = 1) * 10)
#' mcgoftest(varobj = x1, distr = x3,
#'           num.sampl = 999,
#'           stat = "chisq",
#'           seed = 1)
#'
#' ## Chi-square test from R package 'stat' is consistent with 'mcgoftest'
#' chisq.test(rbind(x1, x3), simulate.p.value = TRUE, B = 2e3)$p.value
#'
#' ## We will always fail in to detect differences between 'crude' probability
#' ## vectors without additional information.
#' chisq.test(x= x1, y = x3, simulate.p.value = TRUE, B = 2e3)$p.value
#'

mcgoftest <- function(
                    varobj,
                    distr,
                    pars = NULL,
                    num.sampl = 999,
                    sample.size,
                    stat = c("ks", "ad", "rmse", "chisq", "hd"),
                    breaks = NULL,
                    par.names = NULL,
                    seed = 1,
                    num.cores = 1,
                    tasks = 0,
                    verbose = TRUE) {
   ## A starting permutation script for permutation test used the idea
   ## published on \href{https://goo.gl/hfnNy8}{Alastair Sanderson}
   ## An idea presented in his post: \emph{Using R to analyse data statistical
   ## and numerical data analysis with R}. Herein, it is modified and extended.
   set.seed( seed )
   stat <- match.arg(stat)

   if (distr == "dirichlet" && is.element(stat, c("ks", "ad")))
      stop("\n*** 'ks' and 'ad' approaches are not available for ",
           "n-dimensional distributions")

   if (is.numeric(distr)) {
      if (is.element(stat, c("ks", "ad"))) {
         stat <- "chisq"
         warning("*** 'ks' and 'ad' approaches are not available for ",
                 "discrete distributions",
                 "\n'chisq' approach was applied")
      }
      if (!is.null(dim(distr)))
         stop(
            "*** 'distr' must be a character string naming a distribution ",
            "function e.g, 'norm' will permit accessing ",
            "functions 'pnorm' and 'rnorm', or a numerical vector"
         )
      if (!is.null(dim(varobj)))
         stop("\n*** if 'distr' is a numerical vector,",
               " then 'varobj' must be a numerical vector as well")
      if (length(distr) != length(varobj))
         stop("\n*** if 'distr' is a numeric vector,",
              " then length(distr) == length(varobj)")
      parametric <- FALSE
   }

   if (is.character(distr)) {
      is.discr <- is.element(distr, c("dirichlet", "multinom"))
      pdistr <- paste0("p", distr)
      pdistr <- try(get(pdistr, mode = "function",
                        envir = parent.frame()), silent = TRUE)
      if (inherits(pdistr, "try-error") && !is.discr)
         stop(
            "*** 'distr' must be a character string naming a distribution ",
            "function e.g, 'norm' will permit accessing ",
            "functions 'pnorm' and 'rnorm', or a numerical vector"
         )

      rdistr <- paste0("r", distr)
      rdistr <- try(get(rdistr, mode = "function",
                        envir = parent.frame()), silent = TRUE)

      if (inherits(rdistr, "try-error") && !is.discr)
         stop(
            "*** 'distr' must be a character string naming a ",
            "distribution function e.g, 'norm' will permit accessing ",
            "functions 'pnorm' and 'rnorm'"
         )
      parametric <- TRUE
   }

   if (inherits(varobj, c("matrix", "data.frame")))
      l <- nrow(varobj)
   else
      l <- length(varobj)

   if (missing( sample.size ))
       sample.size = l

   if (sample.size > l && !parametric)
      sample.size = l
   statis <- switch(stat,
          ks = "Kolmogorov-Smirnov",
          ad = "Anderson–Darling",
          rmse = "Root Mean Square",
          chisq = "Pearson's Chi-squared",
          hd = "Hellinger test"
   )

   if (verbose) {
       method <- ifelse(sample.size < length(varobj),
                        "Monte Carlo GoF testing",
                        "Permutation GoF testing")

       approach <- ifelse(parametric, "parametric approach",
                           "non-parametric approach")

       cat( "***", method,"based on" ,
            statis,"statistic",
            "(", approach, ")",
           " ...\n")
   }

   if (stat == "ks" && !parametric && sample.size < length(varobj)) {
      warning("*** ", method," based on ", statis," statistic ",
              "(", approach, ") ", "is not reliable. \n",
              "Use the parametric approach or permutation instead.")
   }

   if (parametric) {
      res <-    GoFtest(x = varobj,
                        R = num.sampl,
                        distr = distr,
                        szise = sample.size,
                        pars = pars,
                        stat = stat,
                        breaks = breaks,
                        parametric = parametric,
                        par.names = par.names,
                        num.cores = num.cores,
                        tasks = tasks,
                        verbose = verbose)
   }
   else {
      res <- tableBoots(rbind(varobj, distr),
                        stat = stat,
                        num.permut = num.sampl )
   }
   return(res)
}

# ======================= Anderson–Darling statistic ========================= #
ad_stat <- function(x, distr, pars = NULL, par.names = NULL) {
   x <- sort(x[complete.cases(x)])
   if (!missing(distr))  x <- distfn(x = x, dfn = distr,
                                     type = "p", arg = pars,
                                     par.names = par.names)

   n <- length(x)
   if (n < 8)
      stop("AD statistic the sample size must be greater than 7")

   h <- x * (1 - rev(x))
   h <- (2 * seq(x) - 1) * log(h)
   return(-n - mean(h))
}


# ============================================================================ #
# ======================= Auxiliary functions ================================ #
# ============================================================================ #

# ================ Auxiliary function to compute the frequencies ============= #

freqs <- function(x, distr, pars = NULL, breaks = NULL,
                    par.names = NULL) {
   n <- length(x)

   if (is.null(breaks))
      breaks <- nclass.FD(x)
   DENS <- hist(x, breaks = breaks, plot = FALSE)
   freq0 <- DENS$counts
   bounds <- data.frame(lower = DENS$breaks[ -length(DENS$breaks)],
                        upper = DENS$breaks[ -1 ])
   up_val <- distfn( x = bounds$upper, dfn = distr,
                     type = "p", arg = pars,
                     par.names = par.names)
   low_val <- distfn(x = bounds$lower, dfn = distr,
                     type = "p", arg = pars,
                     par.names = par.names)
   freq <- round(n * (up_val - low_val))

   return(data.frame(obsf = freq0, expf = freq))
}

# ------------- n-dimensional frequency based on marginals
# So far, only for Dirichlet distribution

freqnd <- function(x, distr, breaks, pars,
                   par.names = NULL) {


   fq <- switch(distr,
            beta = {
                    alfa <- sum(pars)
                    lapply(seq_len(ncol(x)),
                           function(k)
                              freqs(
                               x = x[, k],
                               distr = distr,
                               pars = list( shape1 = pars[k],
                                            shape2 = alfa - pars[k]),
                               breaks = breaks,
                               par.names = par.names))
                 },
            binom = {
                        lapply(seq_len(ncol(x)),
                            function(k)
                               freqs(
                                    x = x[, k],
                                    distr = distr,
                                    pars = list(size = pars[[1]],
                                                prob = pars[[2]][k]),
                                    breaks = breaks,
                                    par.names = par.names))
                    }
   )

   fq <- do.call(rbind, fq)
   rownames(fq) <- NULL
   return(fq)
}

FREQs <- function(x, distr, pars = NULL, breaks = NULL,
                  par.names = NULL) {
   if (is.element(distr, c("dirichlet", "multinom"))) {
      distr <- switch( distr,
                       dirichlet = "beta",
                       multinom = "binom")
      fq <- freqnd(x = x, distr = distr,
                   breaks = breaks, pars = pars,
                   par.names = NULL)
   }
   else
      fq <- freqs(x = x, distr, breaks = breaks,
                  pars = pars, par.names = NULL)
   return(fq)
}

# ====================== Root-Mean-Square statistic ========================== #

rmse <- function(x, distr, pars, breaks = NULL, par.names) {
   freq <- FREQs(x = x, distr = distr, pars = pars,
                    breaks = breaks, par.names = par.names)
   return(sqrt(sum((freq$obsf - freq$expf)^2, na.rm = TRUE))/length(freq$obsf))
}

# =================== Pearson's Chi-squared  statistic ======================= #

chisq <- function(x, distr, pars, breaks = NULL, par.names) {
   freq <- FREQs(x = x, distr = distr, pars = pars,
                breaks = breaks, par.names = par.names)
   if (any(freq$expf == 0)) {
      freq$expf <- freq$expf + 0.5
      freq$obsf <- freq$obsf + 0.5
   }
   return(sum((freq$obsf - freq$expf)^2/freq$expf))
}

# ====================== Kolmogorov-Smirnov  statistic ======================= #

ks_stat <- function(x, distr, pars, par.names) {
   cdf <- function(x, pars) distfn(x = x, dfn = distr, type = "p",
                                   arg = pars, par.names = par.names)
   res <- suppressWarnings(ks.test(x = x, y = cdf, pars))
   return(list(stat = unname(res$statistic), p.value = res$p.value))
}

# --------------------- Auxiliary function to get distribution --------------- #

distfn <- function(x, dfn, type = "r", arg, par.names = NULL) {
    if (!is.null(par.names))  {
       args <- list(x, arg)
       if (type == "r")
          names(args) <- c("n", par.names)
       else
          names(args) <- c("q", par.names)
    } else
       args <-  c(list(x), arg)
    do.call(paste0(type, dfn), args)
}

# ==================== Hellinger divergence  statistic ======================= #

hdiv <- function(x, distr, pars, breaks = NULL, par.names) {

   p <- FREQs(x = x, distr = distr, pars = pars,
            breaks = breaks, par.names = par.names)
   n1 <- sum(p$obsf, na.rm = TRUE)
   n2 <- sum(p$expf, na.rm = TRUE)
   n <- c(n1, n2)

   if (any(p$expf == 0)) {
      p$obsf <- (p$obsf + 1) / (n1 + length(p$obsf))
      p$expf <- (p$expf + 1) / (n2 + length(p$expf))
   }
   else {
      p$expf <- p$expf / n1
      p$obsf <- p$obsf / n2
   }

   w <- (2 * n[1] * n[2]) / (n[1] + n[2])
   sum_hd <- sapply(seq_along(x),
                    function(k) (sqrt(p[k, 1]) - sqrt(p[k, 2]))^2)

   return(w * sum(sum_hd, na.rm = TRUE))
}


# ======================= To compute the statistic ========================== #

# ------- To compute the statistic

stat_fun <- function(a, distr, pars, stat, breaks, parametric,
                     par.names) {
   switch(stat,
          ks = ks_stat(x = a, distr = distr, pars = pars,
                       par.names = par.names)$stat,
          ad = ad_stat(x = a, distr = distr, pars = pars,
                       par.names = par.names),
          rmse = rmse(x = a, distr = distr, pars = pars,
                      par.names = par.names,
                      breaks = breaks),
          chisq = chisq(x = a, distr = distr, pars = pars,
                        par.names = par.names, breaks = breaks),
          hd = hdiv(x = a, distr = distr, pars = pars,
                    par.names = par.names,
                    breaks = breaks)
   )
}


# ------- To iterate the statistic computation

DoIt <- function(r, x, distr, pars, stat, breaks,
                 parametric, par.names, szise) {
   if (parametric) a <- distfn(x = szise, dfn = distr,
                               type = "r", arg = pars,
                               par.names = par.names)

   # to test empirical versus theoretical values
   stat_fun(a = a, distr = distr, pars = pars, stat = stat,
            breaks = breaks, parametric = parametric,
            par.names = par.names)
}

# DoIt(1, distr= distr, pars=pars, stat=stat, breaks=breaks,
#      parametric=parametric, par.names= par.names)

# ======================= GoF test function caller =========================== #

GoFtest <- function(
                    x,
                    R,
                    distr,
                    szise,
                    pars,
                    stat,
                    breaks,
                    parametric,
                    par.names,
                    num.cores,
                    tasks,
                    verbose) {

   if (verbose) progressbar = TRUE else progressbar = FALSE
   if (num.cores > 1) {
      if (Sys.info()['sysname'] == "Linux") {
         bpparam <- MulticoreParam(workers = num.cores, tasks = tasks,
                                   progressbar = progressbar)
      } else bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                                  progressbar = progressbar)

      bstats <- unlist(bplapply(1:R, DoIt, x, distr = distr,
                                pars = pars, stat = stat,
                                breaks = breaks,
                                szise = szise,
                                parametric = parametric,
                                par.names = par.names,
                                BPPARAM = bpparam))
   }

   else {
      bstats <- c()
      for (k in seq_len(R)) {
         bstats <- c(bstats,
                     DoIt(1, x, distr, pars, stat, breaks,
                        parametric, par.names, szise)
                     )
      }
   }


   if (is.element(stat, c("ks", "ad")) && is.numeric(distr))
      stop("\n*** 'ks' and 'ad' are only available for continuous",
           " distribution")

   res <- switch(stat,
                 ks = {
                    statis = ks_stat(x = x, distr = distr, pars = pars,
                                     par.names = par.names)
                    p.value = statis$p.value
                    statis = statis$stat
                    c(KS.stat.D = statis,
                      mc_p.value = mean(c(statis, bstats) >= statis,
                                        na.rm = TRUE),
                      KS.stat.p.value = p.value,
                      sample.size = szise,
                      num.sampl = R)
                 },
                 ad = {
                    statis = ad_stat(x = x, distr = distr, pars = pars,
                                     par.names = par.names)
                    c(AD.stat = statis,
                      mc_p.value=mean(c(statis, bstats) >= statis,
                                      na.rm = TRUE),
                      sample.size = szise,
                      num.sampl = R)
                 },
                 rmse = {
                    statis = rmse(x = x, distr = distr,
                                  pars = pars ,
                                  par.names = par.names,
                                  breaks = breaks)
                    c(rmse = statis,
                      mc_p.value=mean(c(statis, bstats) >= statis,
                                      na.rm = TRUE),
                      sample.size = szise,
                      num.sampl = R)
                 },
                 chisq = {
                    statis = chisq(x = x, distr = distr, pars=pars,
                                   par.names = par.names,
                                   breaks = breaks)
                    c(Chisq = statis,
                      mc_p.value = mean(c(statis, bstats) >= statis,
                                        na.rm = TRUE),
                      sample.size = szise,
                      num.sampl = R)
                 },
                 hd = {
                    statis = hdiv(x = x, distr = distr, pars = pars,
                                  par.names = par.names,
                                  breaks = breaks)
                    c(hd = statis,
                      mc_p.value = mean(c(statis, bstats) >= statis,
                                        na.rm = TRUE),
                      sample.size = szise, num.sampl = R)
                 }
   )
   return(res)
}

# --------------------------- End auxiliary function ------------------------- #
