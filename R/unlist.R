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

#' @rdname unlist
#' @title Flatten lists extended to any class
#' @description Given a list 'x' of R objects from the same class and same
#' format, 'unlist' simplifies it to produce a new R object which contains all
#' the initial components which in 'x' object.
#'
#' @param x Any list of R objects.
#' @details This function try to completely flattening a list. If it fail, then
#' the result for \code{\link[base]{unlist}} function (recursive = TRUE and
#' use.names = TRUE) is returned
#' @export

setGeneric("unlist", signature = "x")
unlist <- function(x) UseMethod("unlist", x)

#' @aliases unlist.default
#' @export
unlist.default <- function(x) {
    x0 <- try(suppressWarnings(do.call("c", unname(x))),
        silent = TRUE
    )
    if (!inherits(x0, "try-error")) {
        x <- x0
    } else {
        x <- base::unlist(x,
            recursive = TRUE,
            use.names = TRUE
        )
    }
    return(x)
}
