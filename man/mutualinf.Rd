% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutualinf.R
\name{mutualinf}
\alias{mutualinf}
\title{Mutual information Based on Multivariate Distributions Constructed
    from Copulas}
\usage{
mutualinf(
  x,
  y,
  copula = NULL,
  margins = NULL,
  paramMargins = NULL,
  method = "ml",
  ties.method = "max"
)
}
\arguments{
\item{x, y}{marginal variates}

\item{copula}{A copula object from class \code{\link[copula]{Mvdc}} or
string specifying all the name for a copula from package
\code{\link[copula]{copula-package}}.}

\item{margins}{A character vector specifying all the parametric marginal
distributions. See details below.}

\item{paramMargins}{A list whose each component is a list (or numeric
vectors) of named components, giving the parameter values of the marginal
distributions. See details below.}

\item{method}{A character string specifying the estimation method to be used
to estimate the dependence parameter(s) (if the copula needs to be
estimated) see \code{\link[copula]{fitCopula}}.}
}
\value{
A list with a data frame carrying the estimated mutual information
    for each (x, y) pair, the joint and marginal probabilities, and the
    "mvdc" copula object.
}
\description{
Computes the mutual information for pairwise x and y marginal
    values based on their multivariate distribution constructed from a
    copula.
}
\details{
The mutual information of a pairwise x and y marginal values is
    defined as:

    \deqn{I{x, y} = log(P(x,y)) - (log(P_1(x)) + log(P_2(y)))}

    where P(x,y) is the multivariate distribution constructed from a
    copula, and P_1(x) and P_2(y) are the marginal CDFs.

    The values \eqn{I{x, y}} expresses a measurement of the relative
    dependece/independece of x and y at the specified point value.

    Notice that the above definition expresses the differences between two
    uncertainty variations. So, for values \eqn{I{x, y} > 0}, we shall say
    that at point (x, y) there is a gain of information for the association
    of the subjacent stochastic processes generating x and y in respect to
    the independent processes. Otherwise, for values \eqn{I{x, y} < 0} we
    shall say that at point (x, y) there is a loss of information for the
    association of the subjacent stochastic process generating x and y in
    respect to the independent processes. Or, equivallently, there is a gain
    of information for the independent processes in respect to
    their association.
}
\examples{
require(stats)
set.seed(12) # set a seed for random number generation
## Random generation of a Normal distributed marginal variate
X <- rnorm(2000, mean = 1, sd = 0.2)

## Random generation of a Weibull-3P distributed marginal variate
Y <- X + rweibull3p(2000, shape = 2, scale = 0.85, mu = 1)

## Correlation test
cor.test(X, Y, method = "spearman")

## Non-linear model fit for 'Y' distribution values
fitY <- fitCDF(Y, distNames = 12) # 3P Weibull distribution model
coefs <- coef(fitY$bestfit) # model coefficients

## Goodness-of-fit test for the  Weibull-3P distribution model
mcgoftest(
    varobj = Y, distr = "weibull3p", pars = coefs, num.sampl = 99,
    sample.size = 1999, stat = "chisq", num.cores = 4, breaks = 200,
    seed = 123
)

## Settngs to estimate the Mutual information
margins <- c("norm", "weibull3p")
parMargins <- list(
    list(mean = 1, sd = 0.2),
    as.list(coefs)
) # Notice "as.list" is used here, not "list"

## Finally estimation of the mutual information
mutual.Inf <- mutualinf(
    x = X, y = Y, copula = "normalCopula",
    margins = margins, paramMargins = parMargins
)
head(mutual.Inf$stat)
## The fitted copula is also returned, so, it can be used in downstream
## analyses
mutual.Inf$copula@copula

}
\seealso{
\code{\link{ppCplot}}, \code{\link{bicopulaGOF}},
    \code{\link[copula]{gofCopula}}, \code{\link{fitCDF}},
    \code{\link[MASS]{fitdistr}}, and \code{\link{fitMixDist}}.
}
