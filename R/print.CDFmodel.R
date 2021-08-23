

## ########################################################################## #
##
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
##
## ########################################################################## #
##

#' @rdname print.CDFmodel
#' @aliases print.CDFmodel
#' @title Printing object from \emph{CDFmodel} class by simple print methods
#' @description An object from \emph{CDFmodel} classes is yielded by
#' function \emph{fitCDF}. This objects carries the information of a fitted
#' non-linear model.
#' @details The definition of these class makes less complex the downstream
#' analyses.
#' @param x Object from class \strong{\emph{CDFmodel}}.
#' @param digits Number of significant digits to be used.
#' @keywords internal
#' @export
print.CDFmodel <- function(x, digits = getOption("digits"), ...) {
    if (length(x$gof) > 2)
        gof <- cbind(
            Adj.R.Square = x$gof[1],
            rho = x$gof[2], R.Cross.val = x$gof[3], x$aic)
    else
        gof <- cbind(
            Adj.R.Square = x$gof[1],
            rho = x$gof[2], x$aic)
    rownames(gof) <- "gof"
    print(summary(x$bestfit), digits)
    cat("\nGoodness of fit:\n")
    print(gof, digits)
    invisible(x)
}


#' @rdname print.CDFmodel
#' @aliases print.CDFmodelList
#' @name print.CDFmodelList
#' @title Printing object from \emph{CDFmodelList} class by simple print methods
#' @description An object from \emph{CDFmodelList} classes is yielded by
#' function \emph{fitCDF}. This objects carries the information of one or more
#' fitted non-linear models.
#' @details The definition of these class makes less complex the downstream
#' analyses.
#' @param x Object from class \strong{\emph{CDFmodel}}.
#' @param digits Number of significant digits to be used.
#' @keywords internal
#' @export
print.CDFmodelList <- function(x, digits = getOption("digits"), ...) {
    cat("List of CDFmodel with", length(x), "elements\n")
    cat("------\n")
    print(x[[ 1L ]])
    cat("\n------\n")
    cat( length(x), "more 'CDFmodel' elements\n")
    invisible(x)
}



