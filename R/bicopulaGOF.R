## Copyright (C) 2019 Robersy Sanchez <https://genomaths.com/>
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

#' @rdname bicopulaGOF
#' @title  Goodness of fit for Bidimensional Copula with Known Margins
#' @description Goodness-of-fit (GOF) tests for a two-dimensional copula based,
#'     by default, on the knowledge of the marginal probability distributions.
#'     Several functionalities/tools from \code{\link[copula]{copula-package}}
#'     are integrated to perform the GOF of copulas that includes specific
#'     margin parameter settings. In terms of
#'     \code{\link[copula]{copula-package}} vocabulary, these are GOF for copula
#'     objects from class \code{\link[copula]{Mvdc}} (also called non-free
#'     copulas).
#' @details Notice that \code{\link[copula]{copula-package}} already have
#'     function \code{\link[copula]{gofCopula}} to perform GOF. However,
#'     its use can be computational expensive for big datasets.
#' @param x Numerical vector with the observations from the first margin
#'     distribution.
#' @param y Numerical vector with the observations from the second margin
#'     distribution.
#' @param copula A copula object from class \code{\link[copula]{Mvdc}} or
#'     string specifying all the name for a copula from package
#'     \code{\link[copula]{copula-package}}.
#' @param margins A character vector specifying all the parametric marginal
#'     distributions. See details below.
#' @param paramMargins a list whose each component is a list (or numeric
#'     vectors) of named components, giving the parameter values of the marginal
#'     distributions. See details below.
#' @param sample.size The size of the samples used for each sampling. It is not
#'     required for the approaches: "Sn", "SnB", and "SnC"; see below.
#' @param nboots The number of booststrap resampling to perform.
#' @param approach a character string specifying the goodness-of-fit test
#'     statistic to be used, which has to be one (or a unique abbreviation) of
#'     following: "adchisq", "adgamma", "Sn", "SnB", "SnC", "chisq". With the
#'     exception of \emph{chisq}, all the other statistics are the same as in
#'     functions \code{\link[copula]{gofTstat}} and
#'     \code{\link[copula]{gofCopula}}. The test using \emph{chisq} implement
#'     the approach described in reference [1].
#' @param Rosenblatt The  Anderson–Darling statistic approach using Rosenblatt
#'     transformation is normally used for the GOF in function
#'     \code{\link[copula]{gofCopula}} from \code{\link[copula]{copula-package}}
#'     package. since, the current function applies a parametric bootstrap
#'     approach generating random variates from the analytical expression for
#'     the margin CDFs, the test does not depend on the theoretical distribution
#'     of the Anderson–Darling statistic. Simulations suggest, so far, that the
#'     application of Rosenblatt transformation may not be needed in this case.
#'     SO, the desicion on whether to apply the Rosenblatt transformation
#'     (computational expensive for big datasets) is left to the users.
#' @param breaks A single number giving the number of bins for the computation
#'     of the Pearson's Chi-squared statistic as suggested in reference [1].
#'     Bascally, it is used to split the unit square [0, 1]^2 into bins/regions.
#' @param method a character string specifying the estimation method to be used
#'     to estimate the dependence parameter(s); see
#'     \code{\link[copula]{fitCopula}}.
#' @param num.cores,tasks Paramaters for parallele computation using package
#'     \code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
#'     use, i.e. at most how many child processes will be run simultaneously
#'     (see \code{\link[BiocParallel]{bplapply}} and the number of tasks per job
#'     (only for Linux OS).
#' @param seed An integer used to set a 'seed' for random number generation.
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @importFrom copula pobs fitCopula mvdc cCopula htrafo describeCop gofTstat
#' @importFrom stats pgamma
#' @return The statistic value estimated for the observations, and the estimated
#'     bootstrap p.value.
#' @export
#' @examples
#' require(stats)
#'
#' set.seed(12)
#' margins = c("norm", "norm")
#' ## Random variates from normal distributions
#' X <- rnorm(2*1e3, mean = 0, sd = 10)
#' Y <- rnorm(2*1e3, mean = 0, sd = 10)
#'
#' parMargins = list( list(mean = 0, sd = 10),
#'                    list(mean = 0, sd = 10))
#'
#' bicopulaGOF(x = X, y = Y, copula = "normalCopula", sample.size = 1e2,
#'             margins = margins, paramMargins = parMargins, nboots = 999,
#'             Rosenblatt = TRUE, approach = "adgamma", breaks = 10,
#'             num.cores = 1L)
#'
#' bicopulaGOF(x = X, y = Y, copula = "normalCopula", sample.size = 1e2,
#'             margins = margins, paramMargins = parMargins, nboots = 999,
#'             Rosenblatt = FALSE, approach = "adgamma", breaks = 10,
#'             num.cores = 1L)
#'
#' ## --- Non-parallel expensive computation ---- -
#' # require(copula)
#' #
#' # U <- pobs(cbind(X, Y)) #' # Compute the pseudo-observations
#' # fit <- fitCopula(normalCopula(), U, method = 'ml')
#' # U <- cCopula(u = U, copula = fit@copula) #' # Rosenblatt transformation
#' #
#' # set.seed(123)
#' # system.time(
#' #   gofCopula(copula = fit@copula, x = U, N = 99, method = "Sn",
#' #             simulation = "pb")
#' # )
#' #
#' ## About
#' ##    user  system elapsed
#' ## 103.370   0.613 105.022
#' #
#' ## --- Parallel computation with 2 cores ---- -
#' ## Same algorithm as in 'gofCopula' adapted for parallel computation
#' # system.time(
#' #   bicopulaGOF(x = X, y = Y, copula = "normalCopula",
#' #               margins = margins, paramMargins = parMargins, nboots = 99,
#' #               Rosenblatt = TRUE, approach = "Sn", breaks = 10, seed = 123,
#' #               num.cores = 2L)
#' # )
#' ## About
#' ##  user  system elapsed
#' ## 2.491   0.100  51.185
bicopulaGOF <- function(x, y, copula = NULL, margins = NULL,
                       paramMargins = NULL, sample.size = NULL, nboots = 10,
                       approach = c("adchisq", "adgamma", "chisq",
                                    "Sn", "SnB", "SnC"),
                       Rosenblatt = FALSE, breaks = 12, method = 'ml',
                       num.cores = 1L, tasks = 0, seed = 123, ...) {
  if (is.null(copula))
       stop("*** A copula or a character string naming a copula must be given")
  if (is.character(copula)) {
       if (is.null(margins))
           stop("*** Provide names of probability distribution margins")
       if (is.null(paramMargins))
           stop("*** Provide parameters for the probability distribution margins")

       # Compute the pseudo-observations for the given data matrix through
       # the margin distributions
       u1 <- distfn(x, dfn = margins[1], type = "p", arg = paramMargins[[1]])
       u2 <- distfn(y, dfn = margins[2], type = "p", arg = paramMargins[[1]])
       U <- cbind(u1, u2)
       U <- pobs(U)
       copula = eval(parse(text=paste0("normalCopula", "()")))

       fit <- fitCopula(copula, U, method = method)
       copula = mvdc(fit@copula, margins = margins, paramMargins = paramMargins)

   } else {
       if (class(copula) != "mvdc")
           stop("*** 'copula' argument must be an object from 'mvdc' class")
   }

   d <- copula@copula@dimension
   approach <- match.arg(approach)

   # ------------------ If not Chi-squared approach ----------- -
   l <- grepl("ad", approach)
   t <- ifelse(is.element(approach, c("Sn", "SnB", "SnC")), TRUE, FALSE)
   if ((l || t) && !is.character(copula)) {
       # Compute the pseudo-observations for the given data matrix through
       # the margin distributions
       u1 <- distfn(x, dfn = copula@margins[1], type = "p",
                   arg = copula@paramMargins[[1]])
       u2 <- distfn(y, dfn = copula@margins[2], type = "p",
                   arg = copula@paramMargins[[1]])
       U <- cbind(u1, u2)
       U <- pobs(U)
   }

   if (l && Rosenblatt)  {
       U <- cCopula(u = U, copula = copula@copula)
   }
   # ------------------------------------------------------------ -

  if (Sys.info()['sysname'] == "Linux") {
    bpparam <- MulticoreParam(workers=num.cores, progressbar = TRUE,
                               tasks=tasks)
  } else bpparam <- SnowParam(workers=num.cores, progressbar = TRUE,
                               type = "SOCK")
   set.seed(seed)
   # ---- Parametric bootstrap with Anderson–Darling statistic ---- -
   if (is.null(sample.size) && !t) stop("*** Please provide a sample size")
   if (l) {
       pstats <- unlist(bplapply(1:nboots, statFun, n = sample.size,
                               copula = copula, d = d, approach = approach,
                               Rosenblatt = Rosenblatt, BPPARAM = bpparam))
   }

   if (approach == "chisq") {
       pstats <- unlist(bplapply(1:nboots, statFun, n = sample.size,
                               copula = copula, d = d, x = x, y = y,
                               approach = approach, breaks = breaks,
                               BPPARAM = bpparam))
   }

   # --- Parametric bootstrap with methods "Sn", "SnB", "SnC", and "Rn" --- -
   if (t) {
       # if (missing(estim.method)) estim.method <- "mpl"
       res <- .gofPB(copula = copula@copula, U, N = nboots,
                       method = approach, num.cores = num.cores,
                       sample.size = sample.size, ...)
   }

   if (!t) {
      res <- switch(approach,
                adchisq = {
                  stat <- ad_stat(x = pchisq(rowSums(qnorm(U)^2), df = d))
                  p.value <- mean(c(stat, pstats) >= stat, na.rm = TRUE)
                  c(AD.stat = stat, mc_p.value = p.value,
                    sample.size = sample.size, num.sampl = nboots)
                  },
                adgamma = {
                  stat <- ad_stat(x = pgamma(rowSums(-log(U)), shape = d))
                  p.value <- mean(c(stat, pstats) >= stat, na.rm = TRUE)
                  c(AD.stat = stat, mc_p.value = p.value,
                    sample.size = sample.size, num.sampl = nboots)
                  },
                chisq = {
                  fq <- freqs.(x = x, y = y, copula = copula, breaks = breaks)
                  stat <- sum((fq$obsf - fq$expf)^2/fq$expf)
                  p.value <- mean(c(stat, pstats) >= stat, na.rm = TRUE)
                  c(Chisq.stat = stat, mc_p.value = p.value,
                    sample.size = sample.size, num.sampl = nboots)
                }
         )
  }
  return(res)
}


# =================== Auxiliary function to compute frequencies ============== #
freqs. <- function(x, y, copula = NULL, breaks = NULL, unifnumb = 2 * 1e4) {
   n <- length(x)
   # r <- cbind(runif(unifnumb), runif(unifnumb))
   if (!is.null(copula)) {
       if (class(copula) != "mvdc")
           stop("*** object 'copula' must be a 'mvdc' class")
   }
   pdistr <- function(u, v, copula) pCopula(cbind(u,v), copula@copula)
   u1 <- distfn(x = x, dfn = copula@margins[1], type = "p",
               arg = copula@paramMargins[[1]])
   u2 <- distfn(x = y, dfn = copula@margins[2], type = "p",
               arg = copula@paramMargins[[1]])
   U <- cbind(u1, u2)
   U <- pobs(U) # Compute the pseudo-observations for the given data matrix

   # To split the interval [0, 1]^2 into segments
   labs <- levels(cut(seq(0,1, 1e-6), breaks))
   bounds <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
                   upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
   bounds[1, 1] <- 0

   # -------- Count events in inside each C-Volumen -------------- -
   # The expected probabilities are computed in the C-Volumen:
   # Vc([a,b] x [c,d]) = C(b,d) - C(a,d) - C(b.c) +  C(a,c) >= 0
   # Nelsen, R. Properties and applications of copulas: A brief survey. in
   # Proc. 1st Braz. Conf. on Stat. Modeling in Insurance and Finance,
   # 1–18 (2003).

   l <- nrow(bounds)
   i = 1
   obsf <- vector("integer", l * l)
   expf <- vector("integer", l * l)
   p <- vector("integer", l * l)
   for (k in 1:l) {
       for(j in 1:l) {
           # ----- Expected counts/p ----- -
           bd <- pdistr(u = bounds[k, 2], v = bounds[j, 2], copula = copula)
           ad <- pdistr(u = bounds[k, 1], v = bounds[j, 2], copula = copula)
           bc <- pdistr(u = bounds[k, 2], v = bounds[j, 1], copula = copula)
           ac <- pdistr(u = bounds[k, 1], v = bounds[j, 1], copula = copula)
           p[i] <- (bd - ad - bc + ac)
           expf[i] <- n * p[i]
           # ----- Observed counts ------ -
           a2b <- (U[, 1] > bounds[k, 1]) & (U[, 1] <= bounds[k, 2])
           c2d <- (U[, 2] > bounds[j, 1]) & (U[, 2] <= bounds[j, 2])
           obsf[i] <- sum(a2b & c2d)
           i = i + 1
       }
   }
   return(data.frame(obsf = obsf, expf = expf, prob = p))
}

# ======== Auxiliary function to compute the statistics for bootstrap ======== #
statFun <- function(r, n, x = NULL, y = NULL, copula, d, approach,
                    breaks = 12, Rosenblatt = FALSE, ...) {
  if (approach != "chisq") {
    U <- rMvdc(n, copula) # Random sampling from copula distribution
    U <- pobs(U)  # Compute the pseudo-observations
    if (Rosenblatt) U <- cCopula(u = U, copula = copula@copula)
  }
  switch(approach,
         adchisq = ad_stat(x = pchisq(rowSums(qnorm(U)^2), df = d)),

         adgamma = ad_stat(x = pgamma(rowSums(-log(U)), shape = d)),

         chisq = {
           idx <- sample.int(n = length(x), size = n)
           fq <- freqs.(x = x[idx], y = y[idx], copula = copula, breaks = breaks)
           sum((fq$obsf - fq$expf)^2/fq$expf)
         }
  )
}
# -------------------------- End auxiliary function -------------------------- #

# ======================= Anderson–Darling statistic ========================= #
ad_stat <- function(x, distr, pars = NULL) {
  x <- sort(x[complete.cases(x)])
  if (!missing(distr))  x <- distfn(x = x, dfn = distr, type = "p", arg = pars)

  n <- length(x)
  if (n < 8)
    stop("To compute AD statistic the sample size must be greater than 7")

  h <- x * (1 - rev(x))
  h <- (2 * seq(x) - 1) * log(h)
  return(-n - mean(h))
}
# -------------------------- End auxiliary function -------------------------- #

# =================== Auxiliary function to get distribution ================= #
distfn <- function(x, dfn, type = "r", arg, log = FALSE,
                   lower.tail = TRUE, log.p = FALSE) {
  switch(type,
         p = do.call(paste0(type, dfn),
                     c(list(x), arg, lower.tail = lower.tail, log.p = log.p)),
         r = do.call(paste0(type, dfn), c(list(x), arg))
  )
}
# -------------------------- End auxiliary function -------------------------- #

# ============================================================================ #
#
# ========== Parallel version of function .gofPB from 'copula' package ========
# <https://github.com/cran/copula/blob/master/R/gofCopula.R>
# mmaechler version 0.999-19       |    22b51c5 on Dec 21, 2018d

## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

.gofPB <- function(copula, x, N, method = c("Sn", "SnB", "SnC"),
                   estim.method = c("mpl", "ml", "itau", "irho", "itau.mpl"),
                   trafo.method = ifelse(method == "Sn", "none",
                                         c("cCopula", "htrafo")),
                   trafoArgs = list(), test.method = c("family", "single"),
                   verbose = interactive(), useR = FALSE,
                   ties = NA, ties.method = c("max", "average", "first",
                                              "last", "random", "min"),
                   fit.ties.meth = eval(formals(rank)$ties.method),
                   sample.size = NULL, num.cores = 1L, tasks = 0,
                   ...) {
  ## Checks -- NB: let the *generic* fitCopula() check 'copula'
  stopifnot(N >= 1)
  if(!is.matrix(x)) {
    warning("coercing 'x' to a matrix.")
    stopifnot(is.matrix(x <- as.matrix(x)))
  }

  stopifnot((d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == d)
  method <- match.arg(method)
  estim.method <- match.arg(estim.method)
  test.method <- match.arg(test.method)
  if(method != "Sn")
    trafo.method <- match.arg(trafo.method, c("cCopula", "htrafo"))
  if(trafo.method == "htrafo") {
    if(!is(copula, "outer_nacopula"))
      stop("'trafo.method' = \"htrafo\" only implemented for copula objects",
           " of type 'outer_nacopula'")
    if(length(copula@childCops))
      stop("currently, only Archimedean copulas are supported")
  }

  ## Input checks
  if (method != "Sn" && trafo.method == "none")
       stop(sprintf("'trafo.method' must be \"cCopula\" or \"htrafo\" with 'method'=\"%s\"", method))
  if (method == "Sn" && trafo.method != "none")
    stop(sprintf("'trafo.method' must be \"none\" with 'method'=\"%s\"", method))

  ## Ties: by default, if at least one column has at least one duplicated entry
  if (is.na(ties <- as.logical(ties))) {
    ties <- any(apply(x, 2, anyDuplicated))
    if (ties)
      warning("argument 'ties' set to TRUE")
  }

  ## Progress bar
  if(verbose) {
    pb <- txtProgressBar(max = N+1, style=if(isatty(stdout())) 3 else 1) # setup progress bar
    on.exit(close(pb)) # on exit, close progress bar
  }

  ## 1) Compute the pseudo-observations
  uhat <- pobs(x, ties.method = ties.method[1])
  uhat.fit <- if (ties == FALSE || ties.method == fit.ties.meth) uhat
  else pobs(x, ties.method = fit.ties.meth)

  ## 2) Fit the copula
  ##    (if test.method = "family", otherwise take the provided copula; this
  ##     is useful for testing random number generators)
  C.th.n <- if(test.method == "family") {
    fitter <- function(..., test.method)
      fitCopula(copula, uhat.fit, method = estim.method,
                estimate.variance = FALSE, ...)@copula
    fitter(...) # avoids passing on 'test.method' to optim()
    # [can be omitted if 'test.method' is a formal arg of gofCopula()];
    # see https://stackoverflow.com/questions/7028385/can-i-remove-an-element-in-dot-dot-dot-and-pass-it-on
  } else copula

  ## 3) Compute the realized test statistic
  # (only) transform if method != "Sn" and trafo.method given
  doTrafo <- (method != "Sn" && trafo.method != "none")
  u <- if(doTrafo) {
    stopifnot(is.list(trafoArgs))
    if(length(names(trafoArgs)) != length(trafoArgs))
      stop("'trafoArgs' must be a fully named list")
    switch(trafo.method,
           "cCopula"= do.call(cCopula, c(list(uhat, copula = C.th.n), trafoArgs)),
           "htrafo" = do.call(htrafo,  c(list(uhat, copula = C.th.n), trafoArgs)),
           stop("wrong transformation method"))
  } else uhat
  T. <- if(method == "Sn") gofTstat(u, method = method, copula = C.th.n,
                                   useR = useR)
  else gofTstat(u, method = method)
  if(verbose) setTxtProgressBar(pb, 1) # update progress bar

  ## 4) Simulate the test statistic under H_0

  ## If ties, get tie structure from x  (FIXME?: what if 'x' has no ties, but 'ties = TRUE' ?? )
  ## IK: If 'x' has no ties, but 'ties = TRUE', extra computations for nothing
  if (ties)
    ir <- apply(x, 2, function(y) rank(sort(y)))
  else ir <- NA

  if (Sys.info()['sysname'] == "Linux") {
    bpparam <- MulticoreParam(workers = num.cores, progressbar = TRUE,
                              tasks = tasks)
  } else bpparam <- SnowParam(workers = num.cores, progressbar = TRUE,
                              type = "SOCK")

  T0 <- unlist(bplapply(1:N, STAT, copula = copula, method = method,
                        ties.method = ties.method, C.th.n = C.th.n,
                        fit.ties.meth = fit.ties.meth, d = d, ir = ir,
                        estim.method = estim.method, ties = ties, useR = useR,
                        test.method = test.method, doTrafo = doTrafo, n = n,
                        trafo.method = trafo.method, trafoArgs = trafoArgs,
                        BPPARAM = bpparam))

  ## 5) Return result object
  tr.string <- if (trafo.method == "none") ""
  else sprintf(", 'trafo.method'=\"%s\"", trafo.method)
  text <- sprintf(", with 'method'=\"%s\", 'estim.method'=\"%s\"%s:",
                  method, estim.method, tr.string)
  structure(class = "htest",
            list(method = paste0(.gofTstr("Parametric", copula, test.method),
                                 text),
                 parameter = c(parameter = getTheta(C.th.n)),
                 statistic = c(statistic = T.),
                 # p.value = (sum(T0 >= T.) + 0.5) / (N + 1), # No!
                 # T. must be count, since it is an observed experimental
                 # output as well!
                 p.value = mean(c(T0, T.) >= T., na.rm = TRUE),
                 data.name = deparse(substitute(x))))
}

## ========== Auxiliary function for informative output ==========
## From copula package
.gofTstr <- function(type, copula, test) {
  paste(type,
        "bootstrap-based goodness-of-fit test of",
        if(test == "single") "a single", # single: *the* exception;
        ## strongly recommended default "family":
        # not mentioned (if only for back-compatib.)
        describeCop(copula, kind = "short"))
}

STAT <- function(k, copula, method, ties.method, fit.ties.meth, estim.method,
                 test.method, doTrafo, trafo.method, trafoArgs, C.th.n,
                 ties, d, ir, n, useR, ...) {

  ## Sample the fitted (if test.method = "family") copula
  U <- rCopula(n, C.th.n)
  if(ties) { ## Sample x may have ties -- Reproduce tie structure of x
    for (i in 1:d) {
      U <- U[order(U[,i]),]
      U[,i] <- U[ir[,i], i]
    }
  }
  Uhat <- pobs(U, ties.method = ties.method[1])
  Uhat.fit <- if (ties == FALSE || ties.method == fit.ties.meth) Uhat
  else pobs(U, ties.method = fit.ties.meth)

  ## Fit the copula (if test.method = "family"; see Step 2))
  C.th.n. <- if(test.method == "family") {
    fitter <- function(..., test.method)
      fitCopula(copula, Uhat.fit, method = estim.method,
                estimate.variance = FALSE, ...)@copula
    fitter(...) # see Step 2)
  } else copula

  ## Compute the test statistic
  u. <- if(doTrafo) { # (no checks needed; all done above)
    switch(trafo.method,
           "cCopula"= do.call(cCopula, c(list(Uhat, copula = C.th.n.),
                                       trafoArgs)),
           "htrafo" = do.call(htrafo,  c(list(Uhat, copula = C.th.n.),
                                       trafoArgs)))
  } else Uhat
  T0. <- if(method == "Sn") gofTstat(u., method = method, copula = C.th.n.,
                                   useR = useR)
  else gofTstat(u., method = method)

  T0. # return
}

