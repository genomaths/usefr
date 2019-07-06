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

#' @rdname mutualinf
#' @name mutualinf
#' @aliases mutualinf
#' @title Mutual information Based on Multivariate Distributions Constructed
#'     from Copulas
#' @description Computes the mutual information for pairwise x and y marginal
#'     values based on their multivariate distribution constructed from a
#'     copula.
#' @details The mutual information of a pairwise x and y marginal values is
#'     defined as:
#'
#'     \deqn{I{x, y} = log(P(x,y)) - (log(P_1(x)) + log(P_2(y)))}
#'
#'     where P(x,y) is the multivariate distribution constructed from a
#'     copula, and P_1(x) and P_2(y) are the marginal CDFs.
#'
#'     The values \eqn{I{x, y}} expresses a measurement of the relative
#'     dependece/independece of x and y at the specified point value.
#'
#'     Notice that the above definition expresses the differences between two
#'     uncertainty variations. So, for values \eqn{I{x, y} > 0}, we shall say
#'     that at point (x, y) there is a gain of information for the association
#'     of the subjacent stochastic processes generating x and y in respect to
#'     the independent processes. Otherwise, for values \eqn{I{x, y} < 0} we
#'     shall say that at point (x, y) there is a loss of information for the
#'     association of the subjacent stochastic process generating x and y in
#'     respect to the independent processes. Or, equivallently, there is a gain
#'     of information for the independent processes in respect to
#'     their association.
#'
#' @param x,y marginal variates
#' @param copula A copula object from class \code{\link[copula]{Mvdc}} or
#'     string specifying all the name for a copula from package
#'     \code{\link[copula]{copula-package}}.
#' @param margins A character vector specifying all the parametric marginal
#'     distributions. See details below.
#' @param paramMargins A list whose each component is a list (or numeric
#'     vectors) of named components, giving the parameter values of the marginal
#'     distributions. See details below.
#' @param method A character string specifying the estimation method to be used
#'     to estimate the dependence parameter(s) (if the copula needs to be
#'     estimated) see \code{\link[copula]{fitCopula}}.
#' @return A list with a data frame carrying the estimated mutual information
#'     for each (x, y) pair, the joint and marginal probabilities, and the
#'     "mvdc" copula object.
#' @seealso \code{\link{ppCplot}}, \code{\link{bicopulaGOF}},
#'     \code{\link[copula]{gofCopula}}, \code{\link{fitCDF}},
#'     \code{\link[MASS]{fitdistr}}, and \code{\link{fitMixDist}}.
#' @export
#'
#' @examples
#' require(stats)
#' set.seed(12) # set a seed for random number generation
#' ## Random generation of a Normal distributed marginal variate
#' X <- rnorm(2000, mean = 1, sd = 0.2)
#'
#' ## Random generation of a Weibull-3P distributed marginal variate
#' Y <- X + rweibull3p(2000, shape = 2, scale = 0.85, mu = 1)
#'
#' ## Correlation test
#' cor.test(X, Y, method = "spearman")
#'
#' ## Non-linear model fit for 'Y' distribution values
#' fitY <- fitCDF(Y, distNames = 12) # 3P Weibull distribution model
#' coefs <- coef(fitY$bestfit) # model coefficients
#'
#' ## Goodness-of-fit test for the  Weibull-3P distribution model
#' mcgoftest(varobj = Y, distr = "weibull3p", pars = coefs, num.sampl = 99,
#'         sample.size =  1999, stat = "chisq", num.cores = 4, breaks = 200,
#'         seed = 123)
#'
#' ## Settngs to estimate the Mutual information
#' margins = c("norm", "weibull3p")
#' parMargins = list(list(mean = 1, sd = 0.2),
#'                 as.list(coefs))  # Notice "as.list" is used here, not "list"
#'
#' ## Finally esitmation of the mutual information
#'  mutual.Inf <- mutualinf(x = X, y = Y, copula = "normalCopula",
#'                         margins = margins,  paramMargins = parMargins )
#' head(mutual.Inf$stat)
#' ## The fitted copula is also returned, so, it can be used in downstream
#' ## analyses
#' mutual.Inf$copula@copula
#'
mutualinf <- function(x, y, copula = NULL, margins = NULL, paramMargins = NULL,
                       method = 'ml', ties.method = "max") {
   if (is.null(copula))
       stop("*** A copula or a character string naming a copula must be given")
   if (is.character(copula)) {
       if (is.null(margins))
           stop("*** Provide names of probability distribution margins")
       if (is.null(paramMargins))
           stop("*** Provide parameters for the probability distribution margins")

       # Compute the pseudo-observations for the given data matrix through
       # the margin distributions
       p1 <- do.call(paste0("p", margins[1]), c(list(x), paramMargins[[1]]))
       p2 <- do.call(paste0("p", margins[2]), c(list(y), paramMargins[[2]]))
       U <- cbind(p1, p2)
       V <- pobs(U)
       copula = eval(parse(text = paste0("copula::", copula, "()")))

       fit <- fitCopula(copula, V, method = method)
       copula = mvdc(fit@copula, margins = margins, paramMargins = paramMargins)
       jprob <- pCopula(u = pobs(U, ties.method = ties.method),
                        copula = copula@copula)
   } else {
       if (class(copula) != "mvdc")
           stop("*** 'copula' argument must be an object from 'mvdc' class")
   }

   if (!is.character(copula)) {
       # Compute the pseudo-observations for the given data matrix through
       # the margin distributions
       p1 <- do.call(paste0("p", copula@margins[1]),
                   c(list(x), copula@paramMargins[[1]]))
       p2 <- do.call(paste0("p", copula@margins[2]),
                   c(list(y), copula@paramMargins[[2]]))
       U <- cbind(p1, p2)
       # U <- pobs(U)
       jprob <- pCopula(u = pobs(U, ties.method = ties.method),
                        copula = copula@copula)
   }

   res <- list()
   res$stat <- data.frame(jprob = jprob, p1 = p1, p2 = p2, x = x, y = y,
                           mInf = logb(jprob) - (logb(p1) +  logb(p2)))
   res$copula <- copula
   return(res)
}

## ======================= Auxiliary function logb ============================
# Although log(0) is not defined for entroupy computation log(0) is made 0.

logb <- function(p, logbase = 2) {
   n <- length(p)
   if (n > 1) {
       logP <- integer(n)
       idx <- p > 0
       logP[idx] <- log(p[idx], base = logbase)
   } else {
       if (p > 0 ) logP <- log(p, base = logbase)
       logP <- 0
  }
  return(logP)
}



