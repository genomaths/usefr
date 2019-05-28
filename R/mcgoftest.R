#' @rdname mcgoftest
#' @title Permutation test for Goodness of fit (GoF)
#' @description To accomplis the nonlinear fit of a probability distribution
#'     function (*PDF*), dIfferent optimization algorithms can be used. Each
#'     algorithm will return a different set of estimated parameter values. AIC
#'     and BIC are not useful (in this case) to decide which parameter set of
#'     values is the best. The goodness-of-fit tests (GOF) can help in this
#'     case.
#' @details The test is intended for continuos distributions. If sampling size
#'     is lesser the size of the sample, then the test becomes a Monte Carlo
#'     test. The thes is based on the use of measures of goodness of fit,
#'     statistics. The following statistics are availible:
#'
#'     \itemize{
#'         \item Kolmogorov- Smirnov statistic (ks). Limitations: sensitive to
#'               ties [1].
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
#'               theoretical value.
#'
#'         \item Pearson's Chi-squared statistic (chisq). Limitation: the sample
#'               must be discretized (partitioned into bins), which is could be
#'               a source of bias that leads to the rejection of the null
#'               hypothesis. Here, the discretization is done using function
#'               the resources from function \code{\link[graphics]{hist}}.
#'
#'         \item Root Mean Square statistic (rmst). Limitation: the same
#'               as 'chisq'.
#'     }
#' @param varobj A a vector containing observations, the variable for which the
#'     CDF parameters was estimated.
#' @param cdf The name of the cummulative distribution function (CDF) or a
#'     concrete CDF from where estimate the cummulative probabilities. The cdf
#'     must be defined in environment-namespace from any package or environment
#'     defined by user.
#' @param pars CDF model parameters. A list of parameters to evaluate the CDF.
#' @param stat One string denoting the statistic to used in the testing: "ks":
#'     Kolmogorov–Smirnov, "ad": Anderson–Darling statistic, "chisq: Pearson's
#'     Chi-squared, and "rmst": Root Mean Square statistic.
#' @param breaks Default is NULL. Basically, the it is same as in function
#'     \code{\link[graphics]{hist}}. If \emph{breaks} = NULL, then function
#'     \code{\link[grDevices]{nclass.FD}} is applied to estimate the breaks.
#' @param num.sampl Number of samplings/permutations. If sample.size < length(x)
#'     , then the test becomes a Monte Carlo test.
#' @param sample.size Size of the samples used for each sampling.
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @importFrom stats ks.test
#' @return A numeric vector with the following data:
#'     \enumerate{
#'
#'         \item Statistic value.
#'
#'         \item mc_p.value: the probability of finding the observed, or more
#'               extreme, results when the null hypothesis \eqn{H_0} of a study
#'               question is true obtained Monte Carlo/permutation approach.
#'     }
#' @export
#' @author Original permutation script idea from
#'     \href{https://goo.gl/hfnNy8}{Alastair Sanderson}. Modified and adapted by
#'     Robersy Sanchez (https://genomaths.com).
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
#' @examples
#' # Example 1
#' # Let us generate a random sample a from a specified Weibull distribution:
#' # Set a seed
#' set.seed( 1 )
#' # Random sample from Weibull( x | shape = 0.5, scale = 1.2 )
#' x = rweibull(10000, shape = 0.5, scale = 1.2)
#'
#' # MC KS test accept the null hypothesis that variable x comes
#' # from Weibull(x | shape = 0.5, scale = 1.2), while the standard
#' # Kolmogorov-Smirnov test reject the Null Hypothesis.
#' mcgoftest(x, cdf = pweibull, pars = c( 0.5, 1.2 ), num.sampl = 500,
#'         sample.size = 1000, num.cores = 4)
#'
#' # Example 2
#' # Let us generate a random sample a random sample from a specified Normal
#' # distribution:
#' # Set a seed
#' set.seed( 1 )
#' x = rnorm(10000, mean = 1.5, sd = 2)
#'
#' # MC KS test accept the null hypothesis that variable x comes
#' # from N(x | mean = 0.5, sd = 1.2), while the standard
#' # Kolmogorov-Smirnov test reject the Null Hypothesis.
#' mcgoftest(x, cdf = pnorm, pars = c(1.5, 2), num.sampl = 500,
#'           sample.size = 1000, num.cores = 1)

mcgoftest <- function(varobj, cdf, pars, num.sampl = 999, sample.size,
                   stat = c("ks", "ad", "rmst", "chisq"), breaks = NULL,
                   num.cores = 1, tasks = 0) {
   # The current version the permutation test for Kolmogorov-Smirnov
   # statistics, is based on \href{https://goo.gl/hfnNy8}{Alastair Sanderson}
   # idea presented in his post: \emph{Using R to analyse data statistical and
   # numerical data analysis with R}.

   stat <- match.arg(stat)

   if(is.character(cdf))
       cdf <- get(cdf, mode = "function", envir = parent.frame())
   if(!is.function(cdf))
       stop("*** 'cdf' must be numeric or a function ",
            "or a string naming a valid function")

   if (missing( sample.size ) || sample.size > round(length(varobj)))
       sample.size = round(length(varobj)/3)

   GoFtest <- function(x, R, szise, cdf, pars, stat, breaks, mc.cores) {
       myfun <- function(a, cdf, pars, stat, breaks) {
           switch(stat,
               ks = suppressWarnings(unname(ks.test(a, cdf, pars)$statistic)),
               ad = ad_stat(x = a, cdf = cdf, pars = pars),
               rmst = rmst(x, cdf, pars, breaks = breaks),
               chisq = chisq(x, cdf, pars, breaks = breaks)
          )
       }

       DoIt <- function(r, cdf, pars, stat, breaks) {
           i <- sample( length( x ), szise )
           # to test empirical versus theoretical values
           myfun(a = x[ i ], cdf = cdf, pars = pars, stat = stat, breaks = breaks)
       }

       if (Sys.info()['sysname'] == "Linux") {
           bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
       } else bpparam <- SnowParam(workers=num.cores, type = "SOCK")

       pstats <- unlist(bplapply(1:R, DoIt, cdf = cdf, pars = pars, stat = stat,
                               breaks = breaks, BPPARAM=bpparam))

       res <- switch(stat,
                   ks = {
                           stat = suppressWarnings(ks.test(x, cdf,pars))
                           p.value = unname(stat$p.value)
                           stat = unname(stat$statistic)
                           c(KS.stat.D = stat,
                               mc_p.value = mean(c(stat, pstats) >= stat,
                                               na.rm = TRUE),
                               KS.stat.p.value = p.value,
                               sample.size = sample.size, num.sampl = num.sampl)
                        },
                   ad = {
                           stat = ad_stat(x = x, cdf = cdf, pars = pars)
                           c(AD.stat = stat,
                               mc_p.value=mean(c(stat, pstats) >= stat,
                                               na.rm = TRUE),
                               sample.size = sample.size, num.sampl = num.sampl)
                        },
                   rmst = {
                           stat = rmst(x = x, cdf = cdf, pars = pars ,
                                       breaks = breaks)
                           c(rmst = stat,
                             mc_p.value=mean(c(stat, pstats) >= stat,
                                             na.rm = TRUE),
                             sample.size = sample.size, num.sampl = num.sampl)
                          },
                   chisq = {
                               stat = chisq(x = x, cdf = cdf, pars = pars,
                                            breaks = breaks)
                               c(rmst = stat,
                                   mc_p.value=mean(c(stat, pstats) >= stat,
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
          rmst = "Root Mean Square",
          chisq = "Pearson's Chi-squared"
   )

   method <- ifelse(sample.size < length(varobj), "Monte Carlo GoF testing",
                    "Permutation GoF testing")

   cat( "***", method,"based on" , statis,"statistic ...\n" )
   GoFtest(x = varobj, R = num.sampl, szise = sample.size, cdf = cdf, pars=pars,
           stat = stat, breaks = breaks, mc.cores = num.cores)
}

# ======================= Anderson–Darling statistic ========================= #
ad_stat <- function(x, cdf, pars) {
   x <- sort(x[complete.cases(x)])
   n <- length(x)
   if (n < 8)
      stop("AD statistic the sample size must be greater than 7")

   logp1 <- log(cdf(x, pars))
   logp2 <- log(1 - cdf(rev(x), pars))
   h <-  (2 * seq(1:n) - 1) * (logp1 + logp2)
   return(-n - mean(h))
}

# ================ Auxiliary funciton to compute the frequencies ============= #
freqs <- function(x, cdf, pars, breaks = NULL) {
   n <- length(x)
   if (is.null(breaks)) breaks <- nclass.FD(x)
   DENS <- hist(x, breaks = breaks, plot = FALSE)
   freq0 <- DENS$counts
   bounds <- data.frame(lower = DENS$breaks[ -length(DENS$breaks)],
                        upper = DENS$breaks[ -1 ])
   freq <- round(n * (cdf(bounds$upper, pars) - cdf(bounds$lower, pars)))
   return(data.frame(obsf = freq0, expf = freq))
}


# ===================== Root-Mean-Square statistic ========================== #
rmst <- function(x, cdf, pars, breaks = NULL) {
   freq <- freqs(x = x, cdf = cdf, pars = pars, breaks = breaks)
   return(sqrt(mean((freq$obsf - freq$expf) ^ 2)))
}

# ================== Pearson's Chi-squared  statistic ======================= #
chisq <- function(x, cdf, pars, breaks = NULL) {
   freq <- freqs(x = x, cdf = cdf, pars = pars, breaks = breaks)
   return(suppressWarnings(chisq.test(x = freq$obsf, y = freq$expf)$statistic))
}



