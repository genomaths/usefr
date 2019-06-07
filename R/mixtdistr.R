#' @rdname mixtdistr
#' @name mixtdistr
#' @aliases rmixtdistr
#' @aliases pmixtdistr
#' @aliases dmixtdistr
#' @aliases qmixtdistr
#' @title Mixture of Distribution Functions
#' @description Density, distribution function, quantile function and random
#'     generation for mixture of distributions
#' @param x, q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be
#'     the number required.
#' @param phi Numerical vector with mixture proportions, where (sum(phi) = 1).
#' @param arg A list of named vectors with the corresponding named distribution
#'     parameters values. The names of the vector of parameters and the
#'     parameter names must correspond to defined functions. For example, if
#'     one of the involved distributions is the gamma density
#'     (\code{\link[stats]{dgamma}}), then the corresponding vector of
#'     parameters must be gamma = c(shape = 'some value', scale = 'some value').
#'     See examples for more details.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ≤ x],
#'     otherwise, P[X > x].
#' @examples
#' set.seed(123) # set a seed for random generation
#' # A mixture of three distributions
#' phi = c(5/10, 3/10, 2/10) # Mixture proportions
#'
#' # Named vector of the corresponding distribution function parameters
#' # must be provided
#' args <- list(gamma = c(shape = 20, scale = 1/10),
#'             weibull = c(shape =  4, scale = 0.8),
#'             lnorm = c(meanlog = 1.2, sdlog = 0.08))
#'
#' #  Sampling from the specified mixture distribution
#' x <- rmixtdistr(n = 1e5, phi = phi , arg = args)
#'
#' # The graphics for the simulated dataset and the corresponding theoretical
#' # mixture distribution
#' hist(x, 100, freq = FALSE)
#' x1 <- seq(0, 10, by = 0.001)
#' lines(x1, dmixtdistr(x1, phi = phi, arg = args), col = "red")

#' @name dmixtdistr
#' @rdname mixtdistr
#' @title Mixture of Distribution Functions
#' @description NULL
#' @details NULL
#' @export
#'
dmixtdistr <- function(x, phi, arg,  log = FALSE,
                       lower.tail = TRUE, log.p = FALSE) {
   k <- length(phi)
   n <- numeric(length(x))
   dfn = names(arg)
   if (length(x) > 1) {
      d <- rowSums(vapply(1:k, function(i)
                                   phi[i] * distfn(x, dfn = dfn[i], type = "d",
                                               arg = arg[[i]], log = log), n))
   } else {
      d <- sum(sapply(1:k, function(i)
                                   phi[i] * distfn(x, dfn = dfn[i], type = "d",
                                               arg = arg[[i]], log = log)))
   }
   return(d)
}

#' @name pmixtdistr
#' @rdname mixtdistr
#' @title Mixture of Distribution Functions
#' @description NULL
#' @details NULL
#' @export
pmixtdistr <- function(q, phi, arg,  lower.tail = TRUE, log.p = FALSE) {
   k <- length(phi)
   n <- numeric(length(q))
   dfn = names(arg)
   if (length(q) > 1) {
      d <- rowSums(vapply(1:k, function(i)
                                   phi[i] * distfn(q, dfn = dfn[i], type = "p",
                                                   arg = arg[[i]],
                                                   lower.tail = lower.tail,
                                                   log.p = log.p), n))
   } else {
      d <- sum(sapply(1:k, function(i)
                                   phi[i] * distfn(q, dfn = dfn[i], type = "p",
                                                   arg = arg[[i]],
                                                   lower.tail = lower.tail,
                                                   log.p = log.p)))
   }
   return(d)
}

#' @name qmixtdistr
#' @rdname mixtdistr
#' @title Mixture of Distribution Functions
#' @description NULL
#' @details NULL
#' @export
qmixtdistr <- function(p, interval = c(0, 1000),
                       phi, arg, lower.tail = TRUE, log.p = FALSE) {
   k <- length(phi)
   n <- numeric(length(p))
   dfn = names(arg)
   qmixtfn <- function(p) {
       uniroot(function(q) {
           ifelse(p <= 0, 0, pmixtdistr(q, phi = phi, arg = arg,
                                       lower.tail = lower.tail,
                                       log.p = log.p) - p)
       }, interval, tol = 1e-10)$root
   }
   qmixtfn <- Vectorize(qmixtfn)
   return(qmixtfn(p))
}

#' @name rmixtdistr
#' @rdname mixtdistr
#' @title Mixture of Distribution Functions
#' @description NULL
#' @details NULL
#' @export
rmixtdistr <- function(n, phi, arg) {
   j <- sample.int(length(phi), n, replace = TRUE, prob = phi)
   k <- length(phi)
   dfn = names(arg)
   freqs <- sapply(1:k, function(i) sum(j == i))
   return(unlist(lapply(1:k, function(i)
                   distfn(freqs[i], dfn = dfn[i], type = "r",
                           arg = arg[[i]]))))
}


# ---------------------------Auxiliary function ------------------------------- #
distfn <- function(x, dfn, type = "d", arg, log = FALSE,
                   lower.tail = TRUE, log.p = FALSE) {
  switch(type,
         d = do.call(paste0(type, dfn), c(list(x), arg, log = log)),
         p = do.call(paste0(type, dfn),
                     c(list(x), arg, lower.tail = lower.tail, log.p = log.p)),
         q = do.call(paste0(type, dfn),
                     c(list(x), arg, lower.tail = lower.tail, log.p = log.p)),
         r = do.call(paste0(type, dfn), c(list(x), arg))
  )
}
# -------------------------- End auxiliar function --------------------------- #





