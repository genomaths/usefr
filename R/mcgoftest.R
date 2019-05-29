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
#'               ties [1]. Only the parametric Monte Carlo resampling or
#'               permutation (provided that there is not ties in the data) can
#'               be used.
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
#'         \item Root Mean Square statistic (rmst). Limitation: the same
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
#'     Chi-squared, and "rmst": Root Mean Square statistic.
#' @param breaks Default is NULL. Basically, the it is same as in function
#'     \code{\link[graphics]{hist}}. If \emph{breaks} = NULL, then function
#'     \code{\link[grDevices]{nclass.FD}} is applied to estimate the breaks.
#' @param parametric Logical object. If TRUE, then samples are drawn from the
#'     theoretical population described by \emph{distr}. Default: TRUE.
#' @param num.sampl Number of samplings/permutations. If sample.size < length(x)
#'     , then the test becomes a Monte Carlo test.
#' @param sample.size Size of the samples used for each sampling.
#' @param seed An integer used to set a 'seed' for random number generation.
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
#' mcgoftest(x, distr = pweibull, pars = c( 0.5, 1.2 ), num.sampl = 500,
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
#' mcgoftest(x, distr = pnorm, pars = c(1.5, 2), num.sampl = 500,
#'           sample.size = 1000, num.cores = 1)

mcgoftest <- function(varobj, distr, pars, num.sampl = 999, sample.size,
                   stat = c("ks", "ad", "rmst", "chisq"), breaks = NULL,
                   parametric = TRUE, seed = 1, num.cores = 1, tasks = 0) {
   # The current version the permutation test for Kolmogorov-Smirnov
   # statistics, is based on \href{https://goo.gl/hfnNy8}{Alastair Sanderson}
   # idea presented in his post: \emph{Using R to analyse data statistical and
   # numerical data analysis with R}.

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
               rmst = rmst(x = a, distr = distr, pars = pars, breaks = breaks),
               chisq = chisq(x = a, distr = distr, pars = pars, breaks=breaks)
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

       DoIt(1, distr=distr, pars=pars, stat=stat, breaks=breaks,
            parametric=parametric)

       if (Sys.info()['sysname'] == "Linux") {
           bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
       } else bpparam <- SnowParam(workers=num.cores, type = "SOCK")

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
                   rmst = {
                           statis = rmst(x = x, distr = distr, pars = pars ,
                                       breaks = breaks)
                           c(rmst = statis,
                             mc_p.value=mean(c(statis, pstats) >= statis,
                                             na.rm = TRUE),
                             sample.size = sample.size, num.sampl = num.sampl)
                          },
                   chisq = {
                               statis = chisq(x = x, distr = distr, pars=pars,
                                            breaks = breaks)
                               c(Chisq = unname(statis),
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
          rmst = "Root Mean Square",
          chisq = "Pearson's Chi-squared"
   )

   method <- ifelse(sample.size < length(varobj), "Monte Carlo GoF testing",
                    "Permutation GoF testing")

   if (sample.size == length(varobj)) parametric <- FALSE
   approach <- ifelse(parametric, "parametric approach",
                       "non-parametric approach")

   cat( "***", method,"based on" , statis,"statistic", "(", approach, ")",
       " ...\n")

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
ad_stat <- function(x, distr, pars) {
   x <- sort(x[complete.cases(x)])
   n <- length(x)
   if (n < 8)
      stop("AD statistic the sample size must be greater than 7")

   logp1 <- log(distfn(x = x, dfn = distr, type = "p", arg = pars))
   logp2 <- log(1 - distfn(x = rev(x), dfn = distr, type = "p", arg = pars))
   h <-  (2 * seq(1:n) - 1) * (logp1 + logp2)
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
rmst <- function(x, distr, pars, breaks = NULL) {
   freq <- freqs(x = x, distr = distr, pars = pars, breaks = breaks)
   return(sqrt(mean((freq$obsf - freq$expf) ^ 2)))
}

# ================== Pearson's Chi-squared  statistic ======================= #
chisq <- function(x, distr, pars, breaks = NULL) {
   freq <- freqs(x = x, distr = distr, pars = pars, breaks = breaks)
   return(suppressWarnings(chisq.test(x = freq$obsf, y = freq$expf)$statistic))
}

# ================== Kolmogorov-Smirnov  statistic ======================= #

ks_stat <- function(x, distr, pars) {
   cdf <- function(x, pars) distfn(x = x, dfn = distr, type = "p", arg = pars)
   res <- suppressWarnings(ks.test(x = x, y = cdf, pars))
   return(list(stat = unname(res$statistic), p.value = res$p.value))
}


# --------------------- Auxiliary function to get distribution --------------- #
distfn <- function(x, dfn, type = "r", arg, log = FALSE,
                   lower.tail = TRUE, log.p = FALSE) {
   switch(type,
          p = do.call(paste0(type, dfn),
                      c(list(x), arg, lower.tail = lower.tail, log.p = log.p)),
          r = do.call(paste0(type, dfn), c(list(x), arg))
   )
}
# -------------------------- End auxiliar function --------------------------- #


