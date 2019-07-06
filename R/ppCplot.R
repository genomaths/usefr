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

#' @rdname ppCplot
#' @title P-P plot of Two-dimensional Copulas
#' @description The function build the P-P plot of Two-dimensional Copulas upon
#'     the knowledge of the margin distribution provided by the user. The
#'     empirical probabilities are computed using function
#'     \code{\link[copula]{empCopula}} from package
#'     \code{[copula-package]{copula}}.
#' @param X Numerical vector with the observations from the first margin
#'     distribution.
#' @details Empirical and theoretical probabilities are estimated using the
#'     quantiles generated with the margin quantile functions. Nonlinear fit of
#'     margin distributions can be previously accomplished using any of the
#'     functions \code{\link{fitCDF}}, \code{\link[MASS]{fitdistr}}, or function
#'     \code{\link{fitMixDist}} for the case where the margins are mixture of
#'     distributions. \emph{npoints} random uniform and iid numbers from the
#'     interval [0, 1] are generated and used to evaluate the quantile margin
#'     distribution functions. Next, the quantiles are used to compute the
#'     empirical and theoretical copulas, which will be used to estimate the
#'     corresponding probabilities.
#' @param Y Numerical vector with the observations from the second margin
#'     distribution.
#' @param copula A copula object from class \code{\link[copula]{Mvdc}} or
#'     string specifying all the name for a copula from package
#'     \code{\link[copula]{copula-package}}.
#' @param margins A character vector specifying all the parametric marginal
#'     distributions. See details below.
#' @param paramMargins A list whose each component is a list (or numeric
#'     vectors) of named components, giving the parameter values of the marginal
#'     distributions. See details below.
#' @param npoints Number of points used to build the P-P plot. The
#' @param method A character string specifying the estimation method to be used
#'     to estimate the dependence parameter(s) (if the copula needs to be
#'     estimated) see \code{\link[copula]{fitCopula}}.
#' @param smoothing character string specifying whether the empirical
#'     distribution function (for F.n()) or copula (for C.n()) is computed (if
#'     smoothing = "none"), or whether the empirical beta copula (smoothing =
#'     "beta") or the empirical checkerboard copula (smoothing = "checkerboard")
#'     is computed (see \code{\link[copula]{empCopula}}.
#' @param ties.method character string specifying how ranks should be computed
#'     if there are ties in any of the coordinate samples of x; passed to
#'     \code{\link[copula]{pobs}} (see \code{\link[copula]{empCopula}}.
#' @param xlab A label for the x axis, defaults to a description of x.
#' @param ylab A label for the y axis, defaults to a description of y.
#' @param glwd Grid line width.
#' @param bgcol Grid background color.
#' @param gcol Grid line color
#' @param dcol Diagonal line color.
#' @param dlwd Diagonal line color.
#' @param xlwd X-axis line width.
#' @param ylwd Y-axis line width.
#' @param xcol X-axis line color.
#' @param ycol Y-axis line color.
#' @param cex.xtitle Cex for x-axis title.
#' @param cex.ytitle Cex for y-axis title.
#' @param padj adjustment for each tick label perpendicular to the reading
#'     direction. For labels parallel to the axes, padj = 0 means right or top
#'     alignment, and padj = 1 means left or bottom alignment. This can be a
#'     vector given a value for each string, and will be recycled as necessary.
#' @param hadj adjustment (see par("adj")) for all labels parallel
#'     (‘horizontal’) to the reading direction. If this is not a finite value,
#'     the default is used (centring for strings parallel to the axis,
#'     justification of the end nearest the axis otherwise).
#' @param xcex,ycex A numerical value giving the amount by which axis labels
#'     should be magnified relative to the default.
#' @param tck The length of tick marks as a fraction of the smaller of the width
#'     or height of the plotting region. If tck >= 0.5 it is interpreted as a
#'     fraction of the relevant side, so if tck = 1 grid lines are drawn. The
#'     default setting (tck = NA) is to use tcl = -0.5.
#' @param tcl The length of tick marks as a fraction of the height of a line of
#'     text. The default value is -0.5; setting tcl = NA sets tck = -0.01 which
#'     is S' default.
#' @param xline,yline On which margin line of the plot the x & y labels must be
#'     placed, starting at 0 counting outwards
#'     (see \code{\link[graphics]{mtext}}).
#' @param xfont,yfont An integer which specifies which font to use for x & y
#'     axes titles (see \code{\link[graphics]{par}}).
#' @param family,lty,bty,col,xlim,ylim,pch,las,mar,font Graphical parameters
#'     (see \code{\link[graphics]{par}}).
#' @param cex A numerical value giving the amount by which plotting text and
#'     symbols should be magnified relative to the default. This starts as 1
#'     when a device is opened, and is reset when the layout is changed, e.g.
#'     by setting mfrow.
#' @param seed An integer used to set a 'seed' for random number generation.
#' @param ... Other graphical parameters to pass to functions:
#'     \code{\link[graphics]{abline}}, \code{\link[graphics]{mtext}} and
#'     \code{\link[graphics]{axis}}.
#' @importFrom copula pobs fitCopula mvdc pMvdc C.n
#' @return The P-P plot and invisible temporary object with the
#'     information to build the graphic which can be assigned to a variable
#'     to use in further plots or analyses.
#' @seealso \code{\link{fitCDF}}, \code{\link[MASS]{fitdistr}},
#'     \code{\link{fitMixDist}}, and \code{\link{bicopulaGOF}}.
#' @export
#' @author Robersy Sanchez (\url{https://genomaths.com}).
#'
#' @examples
#' set.seed(12)
#' margins = c("norm", "norm")
#' ## Random variates from normal distributions
#' X <- rlnorm(200, meanlog =- 0.5, sdlog = 3.1)
#' Y <- rnorm(200, mean = 0, sd = 6)
#' cor(X,Y) ## Correlation between X and Y
#'
#' parMargins = list( list(meanlog = 0.5, sdlog = 3.1),
#'                    list(mean = 0, sd = 10))
#'
#' copula = "normalCopula"
#' npoints = 100
#'
#' ## The information to build the graphic is stored in object 'g'.
#' g <- ppCplot(X = X, Y = Y, copula = "normalCopula", margins = margins,
#'              paramMargins = parMargins, npoints = 20)
#'
ppCplot <- function(X, Y, copula = NULL, margins = NULL, paramMargins = NULL,
               npoints = 100, method = "ml",
               smoothing = c("none", "beta", "checkerboard"),
               ties.method = "max", xlab = "Empirical probabilities",
               ylab = "Theoretical probabilities", glwd = 1.2, bgcol = "grey94",
               gcol = "white", dcol = "red", dlwd = 0.8, tck = NA, tcl = -0.3,
               xlwd = 0.8, ylwd = 0.8, xcol = "black", ycol = "black",
               cex.xtitle = 1.3, cex.ytitle = 1.3, padj = -1, hadj = 0.7,
               xcex = 1.3, ycex = 1.3, xline = 1.6, yline = 2.1, xfont = 3,
               yfont = 3, family = "serif", lty = 1, bty="n", col = "black",
               xlim = c(0, 1), ylim = c(0, 1), pch = 20, las = 1,
               mar = c(4, 4, 2, 1), font=3, cex = 1, seed = 132, ...) {
   if (is.null(copula))
       stop("*** A copula or a character string naming a copula must be given")
   n <- length(X)

   if (is.character(copula)) {
       if (is.null(margins))
           stop("*** Provide names of probability distribution margins")
       if (is.null(paramMargins))
           stop("*** Provide parameters for the margin CDFs")
       if (missing(X)) stop("*** Provide the numerical vector of X values")
       if (missing(Y)) stop("*** Provide the numerical vector of Y values")

       smoothing <- match.arg(smoothing)
       # Compute the pseudo-observations for the given data matrix through
       # the margin distributions
       u <- do.call(paste0("p", margins[1]), c(list(X), paramMargins[[1]]))
       v <- do.call(paste0("p", margins[2]), c(list(Y), paramMargins[[2]]))
       U <- cbind(u, v)
       copula = eval(parse(text=paste0("copula::",copula, "()")))

       V <- pobs(U, ties.method = ties.method)
       fit <- fitCopula(copula, V, method = method)
       copula = mvdc(fit@copula, margins = margins, paramMargins = paramMargins)
   } else {
       if (class(copula) != "mvdc")
           stop("*** 'copula' argument must be an object from 'mvdc' class")
       u <- do.call(paste0("p", copula@margins[1]),
                   c(list(X), copula@paramMargins[[1]]))
       v <- do.call(paste0("p", copula@margins[2]),
                   c(list(Y), copula@paramMargins[[2]]))
       U <- cbind(u, v)
       # U <- pobs(U, ties.method = ties.method)
   }

   set.seed(seed)
   if (missing(npoints) || is.null(npoints)) npoints <- 100
   if (is.numeric(npoints)) {
       u <- do.call(paste0("r", copula@margins[1]),
                    c(list(npoints), copula@paramMargins[[1]]))
       v <- do.call(paste0("r", copula@margins[2]),
                    c(list(npoints), copula@paramMargins[[2]]))
   }

   emprob <- C.n(u = pobs(cbind(u, v), ties.method = ties.method), X = U)
   thprob <- pCopula(u = pobs(cbind(u, v), ties.method = ties.method),
                       copula = copula@copula)
   # pMvdc(x = cbind(u, v), mvdc = copula)

   par(mar = mar, font = font, family = family)
   plot(x = emprob, y = thprob, pch = pch,
       panel.first = {points(0, 0, pch=16, cex=1e6, col = bgcol)
                      grid(col = gcol, lty = lty)},
       xaxt ="n", yaxt = "n", ann = FALSE, bty = bty, col = col, xlim = xlim,
       ylim = ylim, cex = cex, ...)
   abline(a = 0, b =  1, col = dcol, lwd = dlwd, ...)
   axis(1, padj = padj, tck = tck, tcl = tcl, lwd = xlwd, col = xcol,
       cex = xcex, ...)
   axis(2, hadj = hadj, las = las, tck = tck, tcl = tcl, lwd = ylwd,
       col = ycol, cex = ycex, ...)
   mtext(side = 1, text = xlab, line = xline, cex = cex.xtitle,
       font = xfont, ...)
   mtext(side = 2, text = ylab, line = yline, cex = cex.ytitle,
       font = yfont, ...)

   invisible(list(data = data.frame(emprob = emprob, thprob = thprob),
                   copula = copula))
}

