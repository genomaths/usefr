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

#' @rdname AICquasiPoisson
#'
#' @title AIC for Quasi-Poisson glm model
#' @description AIC for Quasi-Poisson glm model dpois
#'
#' @param fitObj a fitted model object
#'
#' @return numeric AIC for Quasi-Poisson glm model
#'
#' @importFrom stats predict coef dpois
#' @keywords internal
.AICquasiPoisson <- function(fitObj) {
   LogLike <- sum(dpois(fitObj$data$count, lambda=exp(predict(fitObj)),
                       log=TRUE))
   return(2 * (length(coef(fitObj)) - LogLike))
}
