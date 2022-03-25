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

#' @rdname coefs
#' @aliases coefs
#' @title Extract Model Coefficients
#' @param object an object for which the extraction of model coefficients is
#' meaningful.
#' @param column Only if \strong{object} is from \code{\link[base]{data.frame}}
#' class.
#' @description This extend \code{\link[stats]{coef}}
#' @details The extreme case is when \strong{object} is from
#' \code{\link[base]{data.frame}}. In this case, each row provides the
#' corresponding coefficient value and the row-names must be provided as well.
#' @param ... Additional parameter not in use yet.
setGeneric(
    "coefs",
    def = function(
        object,
        ...) stats::coef(object, ...))


setOldClass("DirchModel")

#' @rdname coefs
#' @aliases coefs
#' @export
setMethod("coefs", signature(object = "DirchModel"),
    function(
            object,
            output = c("alpha", "marginals")) {

        output <- match.arg(output)
        coefss <- switch(output,
                        alpha = object$alpha,
                        marginals = {
                            coefss <- object$marginals
                            coefss <- lapply(coefss, function(x) {
                                x$Estimate
                            })
                            if (length(coefss) > 1 && is.list(coefss))
                                coefss <- data.frame(do.call(rbind, coefss))
                            else {
                                if (is.list(object))
                                    coefss <- object[[1]]$Estimate
                                else
                                    coefss <- object$Estimate
                            }
                            colnames(coefss) <- c("shape1", "shape2")
                            coefss
                        })
        return(coefss)
    }
)

setOldClass(Classes = "BetaModel")

#' @rdname coefs
#' @aliases coefs
#' @export
setMethod("coefs", signature(object = "BetaModel"),
    function(object) {
        coefs <- object$Estimate
        names(coefs) <- c("shape1", "shape2")
        return(coefs)
    }
)


setOldClass(Classes = "nls.lm")

#' @rdname coefs
#' @aliases coefs
#' @export
setMethod("coefs", signature(object = "nls.lm"),
    function(object) {
        object$par
    }
)

setOldClass(Classes = c("CDFmodel", "CDFmodelList"))

#' @rdname coefs
#' @aliases coefs
#' @export
setMethod("coefs", signature(object = "CDFmodel"),
    function(object) {
        coefs(object$bestfit)
    }
)

#' @rdname coefs
#' @aliases coefs
#' @export
setMethod("coefs", signature(object = "CDFmodelList"),
    function(object) {
        object <- object[ seq_len(length(object) - 1) ]
        sapply(object, function(m) coefs(m))
    }
)

#' @rdname coefs
#' @aliases coefs
#' @export
setMethod("coefs", signature(object = "data.frame"),
    function(object,
             column = 1L) {
        coef <- object[, column]
        names(coef) <- attr(object,"row.names")
        return(coef)
    }
)
