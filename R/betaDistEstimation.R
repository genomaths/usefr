## Copyright (C) 2019 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez
##
## This program is part 'usefr' R package.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
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

#' @rdname betaDistEstimation
#'
#' @title Nonlinear Parameter Estimation for Beta Distribution
#' @description This function perform a rough estimation of the shape
#'     parameters of beta distribution
#' @details In order to obtain the estimates for shape parameters of beta
#' distribution, the squared of the difference between the empirical cumulative
#' distribution function (ecdf) & the theoretical cdf is minimized using the
#' Non-Linear Minimization function \code{\link[stats]{nlm}} 'stats' package.
#'
#' If \code{\link[stats]{nlm}} function fails, then an estimation using
#' \code{\link[stats]{optim}} function is tried.
#'
#' @param q prior probabilities
#' @param init.pars  initial parameter values. Defaults to alpha = 1 &
#' beta = 1, which imply the parsimony pseudo-counts greater than zero.
#' @param force.optim Whether to force the use of \code{\link[stats]{optim}}
#' function for the parameter estimation. Default is FALSE.
#' @param hessian if TRUE, the hessian of f at the minimum is returned.
#' @param method,control,lower,upper,gr (Optional). In the case that
#' \code{\link[stats]{nlm}} function fails, the methods, list of control
#' parameters, and bounders to be used with \code{\link[stats]{optim}} to
#' accomplish the parameter estimation.
#' @param control (Optional). In the case that \code{\link[stats]{nlm}}
#' function fails, a list of control parameters to be used with
#' \code{\link[stats]{optim}} function (see function help: ?optim) accomplish
#' the parameter estimation.
#' @param ... Further parameter for \code{\link[stats]{nlm}} function.
#' @return A list with components, which would vary depending on whether the
#' estimation was performed with \code{\link[stats]{nlm}} or
#' \code{\link[stats]{optim}}. In all the cases the list element carrying
#' the estimated parameters values is named \strong{\emph{parameters}}.
#'
#' @importFrom stats ecdf optim nlm pbeta rbeta
#' @author Robersy Sanchez <https://genomaths.com>
#' @export
#' @seealso \code{\link{betaBinPost}} and \code{\link{estimateDirichDist}}
#' @examples
#' ### A random generation numerical values with Beta distribution
#' x1 = rbeta(n = 1000, shape1 = 2, shape2 = 3)
#'
#' ### Parameter estimation with "nlm" function
#' betaDistEstimation(q = x1, gradtol = 1e-12, hessian = T)
#'
#' ### Parameter estimation with "optim" function
#' betaDistEstimation(q = x1, force.optim = T, hessian = T)

betaDistEstimation <- function(
                               q,
                               init.pars = c(1, 1),
                               force.optim = FALSE,
                               hessian = FALSE,
                               method = "BFGS",
                               gr = NULL,
                               control = list(maxit = 500,
                                              abstol = (10^-8)),
                               lower = -Inf,
                               upper = Inf,
                               ...) {
    ## q: prior probabilities
    ## init.pars: initial parameter values. Defaults to alpha = 1 & beta = 1,
    ## which the parsimony pseudo-counts greater than zero.

    Q <- ecdf(q)  ## Empirical Cumulative Distribution Function
    dat <- data.frame(Q = Q(q), q = q)

    #### Beta fitting In order to obtain the estimates for
    #### shape paramaters the squared of the difference
    #### between the ecdf & the theoretical cdf is
    #### minimized
    min.RSS <- function(data, inits) {
        alpha <- inits[1]  ## This corresponds to the initial
        ## starting parameter for alpha
        beta <- inits[2]  ## This corresponds to the initial
        ## starting parameter for beta Because optim
        ## minimizes a function, the
        with(data, sum((Q - pbeta(q, alpha, beta))^2))
    }

    if (!force.optim) {
        pars <- try(suppressWarnings(
                                    nlm(f = min.RSS,
                                        p = init.pars,
                                        data = dat,
                                        hessian = hessian,
                                        ...)),
                    silent = TRUE)

        if (!inherits(pars, "try-error")) {
            nms <- names(pars)
            nms[nms == "estimate"] <- "parameters"
            names(pars) <- nms
            pars$opt.fun <- "nlm"
        }
    } else pars = init.pars


    if (inherits(pars, "try-error") || force.optim) {

        pars <- try(suppressWarnings(optim( par = init.pars,
                                            fn = min.RSS,
                                            gr = gr,
                                            data = dat,
                                            method = method,
                                            control = control,
                                            hessian = hessian,
                                            lower = lower,
                                            upper = upper)),
            silent = TRUE)
        if (!inherits(pars, "try-error")) {
            nms <- names(pars)
            nms[nms == "par"] <- "parameters"
            names(pars) <- nms
            pars$opt.fun <- "optim"
        }
    }
    return(pars)
}
