#' @name hnorm
#' @aliases hnorm
#' @aliases rhnorm
#' @aliases phnorm
#' @aliases dhnorm
#' @aliases qhnorm
#'
#' @title Half-Normal distribution
#' @description Probability density function (PDF), cummulative density function
#'     (CDF), quantile function and random generation for the Half-normal
#'     (hnorm) distribution.
#' @details An alternative parametrization to avoid issues when sigma is near
#'     zero is applied by using a scaled precision (inverse of the variance)
#'     obtained by setting \eqn{\theta=sqrt(\pi)/\sigma*sqrt(2)}. Details about
#'     these functions can be found in \href{https://goo.gl/yxMF6T}{Wikipedia}
#'     and in \href{http://mathworld.wolfram.com/Half-NormalDistribution.html}{
#'     MathWorld}. Notice that \eqn{\theta = 1} means
#'     \eqn{\sigma = \sqrt \pi/\sqrt 2}.
#' @param q numeric vector
#' @param n number of observations
#' @param theta numerical parameter, strictly positive (default 1).
#' @param sigma Standart deviation of the normal distribution. Here,
#'     \eqn{\sigma = \sqrt \pi/(\theta\sqrt 2)}
#' @param lower.tail logical; if TRUE (default), probabilities are P[X<=x],
#'     otherwise, P[X > x]
#' @param log.p logical; if TRUE, probabilities/densities p are returned as
#'     log(p).
#' @return Half-normal PDF values (theta parameter) for dhnorm,
#'     Half-normal probability for phnorm, quantiles or Half-normal random
#'     generated values for rhnorm.
#' @examples
#' set.seed(123) # set a seed
#' sigma = 1.2
#' theta = sigma2theta(sigma)
#' x <- rhnorm(n = 1e5, theta = theta)
#' hist(x, 100, freq = FALSE)
#' curve(dhnorm(x, theta = theta), col = "red", add = TRUE)
#'
#'#' # Checking the function outputs for the logarithms of probabilities
#' x <- rhalfnorm(n = 10, theta = sigma2theta(2))
#' x1 <- phnorm(x, theta = sigma2theta(2), log = TRUE)
#' x2 <- phnorm(x, theta = sigma2theta(2), log = FALSE)
#' all(round(x1, 8) == round(log(x2), 8))
#'
#' x3 <- dhnorm(x, theta = sigma2theta(2), log = TRUE)
#' x4 <- dhnorm(x, theta = sigma2theta(2), log = FALSE)
#' all(round(x3, 8) == round(log(x4), 8))

#' @name dhnorm
#' @rdname hnorm
#' @title Half-Normal distribution
#' @description NULL
#' @details NULL
#' @export
dhnorm <- function(x, theta = 1, log = FALSE) {
   x <- x * theta * sqrt(2/pi)
   if (log) {
       const <- log(2) + log(2)/2 - log(pi)/2
       d <- ifelse(x < 0, 0, const + log(theta) + dnorm(x, log = TRUE))
   } else d <- ifelse(x < 0, 0, 2 * theta * sqrt(2) * dnorm(x)/sqrt(pi))
   return(d)
}
#'
#' @name phnorm
#' @rdname hnorm
#' @title Half-Normal distribution
#' @description NULL
#' @details NULL
#' @export
phnorm <- function(q, theta = 1, lower.tail = TRUE, log = FALSE) {
   q <- q * theta * sqrt(2/pi)
   # 'p' is given in terms of the error function through 'pnorm'
   # erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE), see ?qnorm
   p <- ifelse(q < 0, 0, 2 * pnorm(q) - 1)
   if (!lower.tail) p = 1 - p
   if (log) p <- log(p)
   return(p)
}

#' @name qhnorm
#' @rdname hnorm
#' @title Half-Normal distribution
#' @description NULL
#' @details NULL
#' @export
qhnorm <- function(p, theta = 1, sigma = NULL, lower.tail=TRUE, log.p=FALSE) {
   if (log.p) p = exp(p)
   if (!lower.tail) p = 1 - p
   p <- (p + 1)/2
   if (is.null(sigma)) sigma <- theta2sigma(theta)
   q <- ifelse(p < 0, 0, qnorm(p, mean = 0, sd = sigma))
   return(q)
}

#' @name rhnorm
#' @rdname hnorm
#' @title Half-Normal distribution
#' @description NULL
#' @details NULL
#' @export
rhnorm <- function(n, theta = 1) {
   r <- abs(rnorm(n, sd = theta2sigma(theta)))
   return(r)
}

#' @name theta2sigma
#' @rdname hnorm
#' @title Half-Normal distribution
#' @description NULL
#' @details NULL
#' @export
theta2sigma <- function(theta) {
  return(sqrt(pi/2)/theta)
}

#' @name sigma2theta
#' @rdname hnorm
#' @title Half-Normal distribution
#' @description NULL
#' @details NULL
#' @export
sigma2theta <- function(sigma) {
  return(sqrt(pi/2)/sigma)
}





