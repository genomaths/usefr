#' @rdname mcgoftest
#' @title Permutation test for Goodness of fit
#' @description To accomplis the nonlinear fit of a probability distribution
#'     function (*PDF*), dIfferent optimization algorithms can be used. Each
#'     algorithm will return a different set of estimated parameter values. AIC
#'     and BIC are not useful (in this case) to decide which parameter set of
#'     values is the best. The goodness-of-fit tests (GOF) can help in this
#'     case. In addition, when the sample size goes beyond 150, standard GOFs
#'     can fail. In this case, permutation test can help to confront the issue.
#' @details For the current version the permutation test for Kolmogorov-Smirnov
#'     statistics, is based on \href{https://goo.gl/hfnNy8}{Alastair Sanderson}
#'     idea presented in his post: \emph{Using R to analyse data statistical and
#'     numerical data analysis with R}.
#' @param varobj A a vector containing observations, the variable for which the
#'     CDF parameters was estimated.
#' @param cdf The name of the cummulative distribution function (CDF) or a
#'     concrete CDF from where estimate the cummulative probabilities.
#' @param pars CDF model parameters. A list of parameters to evaluate the CDF
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
#'         \item mc_p.value: the probability of finding the observed, or more
#'               extreme, results when the null hypothesis \eqn{H_0} of a study
#'               question is true obtained Monte Carlo/permutation approach.
#'         \item KS: Kolmogorov-Smirnov statistic obtained with the application
#'               of the function \code{\link[stats]{ks.test}} on the
#'               variable *varobj*.
#'         \item KS_p.value: the probability of finding the observed, or more
#'               extreme, results when the null hypothesis \eqn{H_0} of a study
#'               question is true obtained with Kolmogorov-Smirnov test.
#'     }
#' @export
#' @author Original script from \href{https://goo.gl/hfnNy8}{Alastair Sanderson}
#'     . Modified and adapted by Robersy Sanchez (https://genomaths.com).
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
                   num.cores = 1, tasks=0) {

   if (missing( sample.size ) || sample.size >= round(length(varobj)))
       sample.size = round(length(varobj)/3)

   KStest <- function(x, R, szise, cdf, pars, mc.cores) {
       myfun <- function(a, cdf, pars) unname(ks.test(a, cdf, pars)$statistic)

       DoIt <- function(r, cdf, pars) {
           i <- sample( length( x ), szise )
           # to test empirical versus theoretical values
           myfun(x[ i ], cdf=cdf, pars=pars )
       }

       if (Sys.info()['sysname'] == "Linux") {
           bpparam <- MulticoreParam(workers=num.cores, tasks=tasks)
       } else bpparam <- SnowParam(workers=num.cores, type = "SOCK")

       pstats <- unlist(bplapply(1:R, DoIt, cdf=cdf, pars=pars,
                               BPPARAM=bpparam))
       stat <- ks.test(varobj, cdf, pars)
       p.val <- stat$p.value
       stat <- stat$statistic
       return(c(mc_p.value=mean(c(stat, pstats) >= stat, na.rm = TRUE),
           KS.stat = stat, KS_p.value = p.val))
   }

  cat( "*** Monte Carlo KS test,..\n" )
  KStest(varobj, R = num.sampl, szise = sample.size, cdf=cdf, pars=pars,
      mc.cores = num.cores)
}
