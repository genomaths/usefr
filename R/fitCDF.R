#' @rdname fitCDF
#' @title Nonlinear fit of a commulative distribution function
#' @description Usually the parameter estimation of a cummulative distribution
#'     function (*CDF*) are accomplished using the corresponding probability
#'     density function (*PDF*). DIfferent optimization algorithms can be used
#'     to accomplished this task and different algorithms can yield different
#'     esitmated parameters. Hence, why not try to fit the CDF directly?
#' @param varobj A a vector containing observations, the variable for which the
#'     CDF parameters will be estimated.
#' @param distNames a vector of distribution numbers to select from the listed
#'     below in details section, e.g. c(1:10, 15)
#' @param plot Logic. Default TRUE. Whether to produce the plots for the best
#'     fitted CDF.
#' @param plot.num The number of distributions to be plotted.
#' @param only.info Logic. Default TRUE. If true, only information about the
#'     parameter estimation is returned.
#' @param maxiter,maxfev,ptol Parameters to ontrol of various aspects of the
#'     Levenberg-Marquardt algorithm through function
#'     \code{\link[minpack.lm]{nls.lm.control}} from *minpack.lm* package.
#' @param verbose Logic. If TRUE, prints the function log to stdout
#' @details The nonlinear fit (NLF) problem for CDFs is addressed with
#'     Levenberg-Marquardt algorithm implemented in function
#'     \code{\link[minpack.lm]{nls.lm}} from package *minpack.lm*. This function
#'     is inspired in a script for the function
#'     \code{\link[propagate]{fitDistr}} from the package propagate [1]. Some
#'     parts or script ideas from function \code{\link[propagate]{fitDistr}}
#'     are used, but here we to estimate CDF and not the PDF as in the case of
#'     "\code{\link[propagate]{fitDistr}}. A more informative are also
#'     incorporated. The studentized residuals are provided as well.
#'     The list (so far) of possible CDFs is:
#'     \enumerate{
#'         \item Normal \href{https://goo.gl/xaEAdT}{(Wikipedia)}
#'         \item Log-normal \href{https://goo.gl/a7MtYq}{(Wikipedia)}
#'         \item Half-normal \href{https://goo.gl/yxMF6T}{(Wikipedia)}. An
#'             Alternatively using a scaled precision (inverse of the variance)
#'             parametrization (to avoid issues if \eqn{\sigma} is near zero),
#'             obtained by setting \eqn{\theta=sqrt(\pi)/\sigma*sqrt(2)}.
#'         \item Generalized Normal \href{https://goo.gl/EPk8mH}{(Wikipedia)}
#'         \item T-Generalized Normal [2].
#'         \item Laplace \href{https://goo.gl/fCykV9}{(Wikipedia)}
#'         \item Gamma \href{https://goo.gl/cYkvar}{(Wikipedia)}
#'         \item 3P Gamma [3].
#'         \item Generalized 4P Gamma [3]
#'                 \href{https://goo.gl/1n4kpW.}{(Wikipedia)}
#'         \item Generalized 3P Gamma [3].
#'         \item Weibull \href{https://goo.gl/WMXmQP}{(Wikipedia)}
#'         \item 3P Weibull \href{https://goo.gl/WMXmQP}{(Wikipedia)}
#'         \item Beta \href{https://goo.gl/893wzR}{(Wikipedia)}
#'         \item 3P Beta \href{https://goo.gl/893wzR}{(Wikipedia)}
#'         \item 4P Beta \href{https://goo.gl/893wzR}{(Wikipedia)}
#'         \item Beta-Weibull \href{https://goo.gl/dpaG8h}{ReliaWiki}
#'         \item Generalized Beta \href{https://goo.gl/UcVsct}{(Wikipedia)}
#'         \item Rayleigh \href{https://goo.gl/d9b3zv}{(Wikipedia)}
#'         \item Exponential \href{https://goo.gl/stVsi7}{(Wikipedia)}
#'         \item 2P Exponential \href{https://goo.gl/stVsi7}{(Wikipedia)}
#'     }
#' @return After return the plots, a list with following values is provided:
#'     \itemize{
#'         \item aic: Akaike information creterion
#'         \item fit: list of results of fitted distribution, with parameter
#'             values
#'         \item bestfit: the best fitted distribution according to AIC
#'         \item fitted: fitted values from the best fit
#'         \item rstudent: studentized residuals
#'         \item residuals: residuals
#'     }
#'     After x = fitCDF( varobj, ...), attributes( x$bestfit ) yields:
#'     $names
#'     [1] "par" "hessian" "fvec" "info" "message" "diag" "niter" "rsstrace"
#'         "deviance"
#'     $class
#'     [1] "nls.lm"
#'     And fitting details can be retrived with summary(x$bestfit)
#' @importFrom numDeriv grad
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @importFrom utils flush.console
#' @importFrom stats var sd quantile ecdf pgamma pnorm pbeta pexp pweibull
#'     plnorm na.omit splinefun qqnorm qqline
#' @importFrom graphics par grid lines mtext abline
#' @export
#' @author Robersy Sanchez (https://genomaths.com)
#' @references
#'   \enumerate{
#'     \item Andrej-Nikolai Spiess (2014). propagate: Propagation of
#'           Uncertainty. R package version 1.0-4.
#'           http://CRAN.R-project.org/package=propagate
#'     \item Abramowitz, M. and Stegun, I. A. (1972) Handbook of Mathematical
#'           Functions. New York: Dover. Chapter 6: Gamma and Related Functions.
#'     \item Hand-book on STATISTICAL DISTRIBUTIONS for experimentalists
#'           (pag 73) by Christian Walck. Particle Physics Group Fysikum.
#'           University of Stockholm (e-mail: walck@physto.se).
#'   }
#' @examples
#' set.seed(1230)
#' x1 <- rnorm(10000, mean = 0.5, sd = 1)
#' cdfp <- fitCDF(x1, distNames = "Normal", plot = FALSE)
#' summary(cdfp$bestfit)

fitCDF <- function (varobj, distNames, plot = TRUE, plot.num = 1,
                    only.info = FALSE, maxiter = 10^4, maxfev = 1e+5,
                    ptol = 1e-12, verbose = TRUE) {

   # "count" function from "propagate"
   counter <- function (i) {
       if (i%%10 == 0) cat(i) else cat(".")
       if (i%%50 == 0) cat("\n")
       flush.console()
   }

   nls.rstudent <- function( grad, residuals, num.par ) {
       # *** Studentized residuals for nls *** #
       # grad: derivative of model as obtained from function deriv
       # R.S
       n = length( residuals )
       h = diag( grad %*% solve( crossprod( grad ) ) %*% t( grad ) )
       s = sqrt( n*var( residuals, na.rm = T )/( n - num.par ) )
       residuals/( s * sqrt( 1 - h ) )
   }

   options(warn = -1)
   if (is.vector(varobj)) X <- sort( varobj )
   else stop( "varobj must be a numeric vector!" )
   MEAN <- mean( X, na.rm = TRUE)
   VAR <- var( X, na.rm = TRUE)
   SD <- sd(X, na.rm = TRUE)
   MIN <- min( X, na.rm = TRUE)
   MAX <- max( X, na.rm = TRUE)
   Q = unname( quantile( X, 0.632, na.rm = TRUE ) )
   Fy = ecdf( X )
   pX = Fy( X )

   # *** Generalized normal CDF *** #
   # https://en.(Wikipedia).org/wiki/Generalized_normal_distribution
   pgnorm <- function( q , mean = 0, sigma = 1, beta = 1 )  {
       y = ( abs( q - mean )/sigma )^beta
       # 1/2 + sign( q - mean ) * pgamma( y, 1/beta )/( 2*gamma( 1/beta ) )
       1/2 + sign( q - mean ) * pgamma( y, 1/beta )/2
   }
   # *** Half-normal CDF *** #
   phalfnorm <- function(q, theta = 1) {
       # Alternatively using a scaled precision (inverse of the variance)
       # parametrization (to avoid issues if sigma is near zero),
       # obtained by setting
       ifelse(q < 0, 0, 2 * pnorm(q * theta * sqrt(2)/sqrt(pi)) - 1)
   }

   # *** Thermodynamic based Generalized normal CDF *** #
   ptgnorm <- function( q , mean = 0, sigma = 1, beta = 1 ) {
       # pgamma is closely related to the incomplete gamma function.
       # As defined by Abramowitz and Stegun 6.5.1 (and by 'Numerical Recipes')
       # this is: P(a,x) = 1/Gamma(a) integral_0^x t^(a-1) exp(-t) dt
       # P(a, x) is pgamma(x, a).
       # Abramowitz, M. and Stegun, I. A. (1972) Handbook of Mathematical
       # Functions. New York: Dover. Chapter 6: Gamma and Related Functions.
       y = ( abs( q - mean )/( sigma ) )^beta
       1/2 + sign( q - mean ) * pgamma( y, 1/beta )/2
   }

   # *** Definition of Laplace CDF *** #
   plaplace <- function( q , mean = 0, sigma = 1 ) {
       1/sqrt( 2 ) + 1/sqrt( 2 ) * ( 1 - exp( -abs( q - mean )/sigma ) )
   }

   # *** Definition of 3P gamma CDF *** #
   pgamma3 <- function( q, shape, rate, mu ) pgamma( q - mu, shape, rate )

   # *** Definition of generalized gamma distribution *** #
   pggamma4 <- function (q, alpha = 1, scale = 1, mu = 0, psi = 1,
                       lower.tail = TRUE, log.p = FALSE ) {
       # Hand-book on STATISTICAL DISTRIBUTIONS for experimentalists ( pag 73 )
       # by Christian Walck. Particle Physics Group Fysikum
       # University of Stockholm (e-mail: walck@physto.se )
       y <- ( ( q - mu )/scale )^alpha
       # See section Note at ?pgamma
       p <- pgamma( y, psi, lower.tail = lower.tail, log.p = log.p )
       p[ alpha < 0 ] <- NaN
       return(p)
   }

   pggamma3 <- function(q, alpha = 1, scale = 1, psi = 1, lower.tail = TRUE,
                       log.p = FALSE )  {
       # Hand-book on  STATISTICAL DISTRIBUTIONS for experimentalists (pag 73)
       # by Christian Walck. Particle Physics Group Fysikum. University of
       # Stockholm (e-mail: walck@physto.se)
       y <- ( q / scale )^alpha
       # See section Note at ?pgamma
       p <- pgamma( y, psi, lower.tail = lower.tail, log.p = log.p )
       p[ alpha < 0 ] <- NaN
       return(p)
   }

   # *** Generalized beta CDF *** #
   pgbeta <- function(q, shape1, shape2, lambda = 1, lower.tail = TRUE,
                       log.p = FALSE ) {
       ans <- pbeta( q = 1/(1 + (1/q - 1)/lambda ), shape1 = shape1,
                   shape2 = shape2, lower.tail = lower.tail, log.p = log.p )
       ans[ lambda <= 0 ] <- NaN
       return(ans)
   }

   # *** 3P beta CDF *** #
   pbeta3 <- function(q, shape1 = 2, shape2 = 3, a = 0, lower.tail = TRUE,
                       log.p = FALSE ) {
       pbeta( q - a, shape1, shape2, lower.tail = lower.tail, log.p = log.p )
   }

   # *** 4P beta CDF *** #
   pbeta4 <- function(q, shape1 = 2, shape2 = 3, a = 0, b = 1,
                       lower.tail = TRUE, log.p = FALSE ) {
       pbeta( ( q - a )/( b - a ), shape1, shape2, lower.tail = lower.tail,
               log.p = log.p )
   }

   # *** Definition of 3P Weibull *** #
   pweibull3 <- function(q, shape, scale, mu) pweibull(q - mu, shape, scale)

   # *** Definition of Rayleigh distribution *** #
   prayleigh <- function(q, sigma )  1 - exp( -q^2/( 2 * sigma^2 ))

   # *** Definition of 2P Exponential *** #
   pexp2 <- function(q, rate, mu) pexp(q - mu , rate = 1, lower.tail = TRUE,
                                       log.p = FALSE)

   # *** Definition of Beta-Weibull distribution *** #
   # It is taken from R pacakge "Newdistns"
   pbweibull <- function( q, alpha, beta, shape, scale )
       pbeta(pweibull(q, shape = shape, scale = scale), shape1 = alpha,
           shape2 = beta)

   distNAMES <- c("Normal", "Log-normal", "Half-Normal", "Generalized Normal",
               "T-Generalized Normal", "Laplace", "Gamma", "3P Gamma",
               "Generalized 4P Gamma", "Generalized 3P Gamma","Weibull",
               "3P Weibull", "Beta", "3P Beta", "4P Beta", "Beta-Weibull",
               "Generalized Beta", "Rayleigh", "Exponential",
               "2P Exponential")

   funLIST <- list(pnorm, plnorm, phalfnorm, pgnorm, ptgnorm, plaplace, pgamma,
                   pgamma3, pggamma4, pggamma3, pweibull, pweibull3, pbeta,
                   pbeta3, pbeta4, pbweibull, pgbeta, prayleigh, pexp, pexp2 )

   parLIST <- list(norm = c( mean = MEAN, sd = SD ),
                   lnorm = c(meanlog = mean( log1p( X ), na.rm = TRUE ),
                           sdlog = sd( log1p( X ), na.rm = TRUE ) ),
                   halfnorm = c(theta = sqrt(pi/2)),
                   gnorm = c( mean = MEAN, sigma = SD, beta = 2 ),
                   tgnorm = c( mean = MEAN, sigma = SD, beta = 2 ),
                   laplace = c( mean = MEAN, sigma = sqrt( VAR ) ),
                   gamma = c( shape = MEAN^2/VAR, rate = MEAN/VAR ),
                   gamma3 = c( shape = MEAN^2/VAR, rate = MEAN/VAR, mu = 0 ),
                   ggamma4 = c(alpha = MEAN^2/VAR, scale = VAR/MEAN, mu = MIN,
                               psi = 1 ),
                   ggamma3 = c( alpha = MEAN^2/VAR, scale = VAR/MEAN, psi = 1 ),
                   weibull = c( shape = log( 2 ), scale = Q ),
                   weibull3 = c( mu = MIN, shape = log( 2 ), scale = Q ),
                   beta = c(shape1 = 1, shape2 = 2 ),
                   beta3 = c(shape1 = 1, shape2 = 2, a = MIN ),
                   beta4 = c(shape1 = 2, shape2 = 3, a=0.9 * MIN, b=1.1 * MAX),
                   bweibull = c(alpha=1, beta=2, shape = log( 2 ), scale = Q ),
                   gbeta = c( shape1 = 1, shape2 = 2, lambda = 1 ),
                   rayleigh = c( sigma = SD ),
                   exp = c( rate = 1 ),
                   exp2 = c( rate = 1, mu = 0 )
   )

   if( missing( distNames ) ) distNames = 1:length(distNAMES)
   if( length( distNames ) != length( distNAMES ) ) {
       if (is.character(distNames))
           distNames = as.numeric(na.omit(match(distNames, distNAMES)))
       if (sum(!is.element(distNames, 1:length(distNAMES))) > 0)
           stop("At least one CDF is not found between the possible selections")
       distNAMES = distNAMES[ distNames ]
       funLIST = funLIST[ distNames ]
       parLIST = parLIST[ distNames ]
   }

   optFun <- function(par, probfun, quantiles, prob, eval = FALSE) {
       START <- as.list(par)
       START$q <- quantiles
       EVAL <- try(do.call( probfun, START ), silent = TRUE)
       if (inherits( EVAL, "try-error" ) ) return( NA )
       EVAL[ is.nan( EVAL ) ] <- 0
       # nls.lm will minimize the sum of squares of vector RSS
       RSS <- abs(prob - EVAL)
       if( eval ) return( EVAL )
       else return( RSS )
   }

   fitAIC <- function( fitobj ) {
       RESID <- fitobj$fvec
       sse = sum( RESID^2, na.rm = TRUE )
       N <- length( RESID )
       N * (1 + log(2 * pi) + log(sse/N)) + 2 * (1L + length( fitobj$par))
   }

   fitLIST <- vector("list", length = length(distNAMES))
   AICS <- rep(NA, length(distNAMES))

   for (i in 1:length(distNAMES)) {
       if (verbose) cat("Fitting ", distNAMES[i], " distribution..", sep = "")
       GRID <- matrix(parLIST[[i]], nrow = 1)
       colnames(GRID) <- names(parLIST[[i]])
       rssVEC <- rep(NA, nrow(GRID))
       for (j in 1:nrow(GRID)) {
           if (verbose) counter(j)
           PARS <- GRID[j, ]
           FIT <- try(nls.lm(par = PARS, fn = optFun, probfun = funLIST[[i]],
                           quantiles = X, prob = pX,
                           control = nls.lm.control(maxiter=maxiter,
                                                    maxfev=maxfev, ptol=ptol)),
                   silent = TRUE)
           if (inherits(FIT, "try-error")) rssVEC[j] <- NA
           else rssVEC[j] <- FIT$deviance
       }
       if (!is.na(rssVEC[ j ]) ) {
           WHICH <- which.min( rssVEC[ j ] )
           bestPAR <- GRID[WHICH, ]
           FIT <- try(nls.lm( par = bestPAR, fn = optFun, probfun=funLIST[[i]],
                           quantiles = X, prob = pX,
                           control = nls.lm.control(maxiter=maxiter,
                                                    maxfev=maxfev, ptol=ptol)),
                   silent = TRUE)
       }
       if (inherits(FIT, "try-error")) {
           FIT <- NA
           if (verbose) cat("Error!\n")
       }
       else {
           fitLIST[[i]] <- FIT
           AICS[i] <- tryCatch( fitAIC( FIT ), error = function(e) NA )
           if (verbose) cat("Done.\n")
       }
   }

   ORDER <- order(AICS)
   aicDAT <- data.frame( Distribution = distNAMES, AIC = AICS)
   aicDAT <- aicDAT[ ORDER, ]

   distNAMES = distNAMES[ ORDER ]
   fitLIST = fitLIST[ ORDER ]
   bestFIT = fitLIST[[ 1 ]]
   bestFIT$info = distNAMES[ 1 ]

   if( only.info ){
       return(c(Distribution=bestFIT$info, Deviance=bestFIT$deviance,
                 bestFIT$par, aicDAT))
   } else {
       if( plot ) {
           for(k in 1:min(plot.num, length(distNames))) {
               cat(" * Estimating Studentized residuals for",
                   distNAMES[ k ], "distribution\n" )
               SEL <- ORDER[ k ]
               FITs <- fitLIST[[ k ]]
               evalLIST <- as.list( FITs$par )
               evalLIST$q <- X
               evalY <- do.call( funLIST[[ SEL ]], evalLIST )

               #   # Derivative of the best fit function
               evalLST <- as.list( FITs$par )
               evalLST$func <- funLIST[[ SEL ]]
               evalLST$x <- X
               evalLST$method <- "simple"
               gradient <- do.call(grad, evalLST )

               # NaN/missing values are replace by cubic spline interpolation
               grad.spline = splinefun( X, gradient )
               ind = which(is.na( gradient))
               gradient[ ind ] <- grad.spline( X[ ind ] )
               residuals = pX - evalY

               rstudent = nls.rstudent(gradient, residuals, length(bestFIT$par))
               if (k == 1) {
                   rbestFIT = residuals
                   rstbestFIT = rstudent
                   YbestFIT = evalY
               }

               cat( " * Plots for", distNAMES[ k ], "distribution...\n" )
               par(mfrow=c(2, 2), mar= c(2,2.5,2.2,1) + 0.2, mgp=c(1.2,0.4,0))
               plot(Fy, verticals=TRUE, col="blue", pch=20,
                   xlab=expression(italic( "x" )), ylab="CDF",
                   main = paste( distNAMES[ k ], "Distribution"), cex.main=0.9)
               grid( NULL,NULL, lty = 6, col = "cornsilk2" )
               lines( evalLIST$q, evalY, col = 2, lty = 2, lwd = 2)
               mtext(text=paste("AIC =", round(AICS[ SEL ], 3 ) ), cex = 0.6)

               ## PP-plot
               par(mar = c( 2, 2, 2.2, 1 ) + 0.2, mgp = c( 1.2, 0.4, 0 ) )
               plot(pX, evalY, col = "blue" , main = "P-P Plot", pch = 20,
                   cex = 0.4, xlab = "Empirical CDF",
                   ylab = "Theoretical CDF", cex.main = 0.9)
               grid( NULL,NULL, lty = 6, col = "cornsilk2" )
               abline( 0, 1, col= "red", lwd = 2 ) # Reference line y = x
               pars = FITs$par
               pars=t( cbind( names( pars ), format( round ( pars, 3 ), 3 )))
               mtext( text = paste( pars, collapse = " " ), cex = 0.6 )

               par(mar = c( 3.5, 3, 1, 1 ) + 0.07, mgp = c( 1.2, 0.4, 0 ) )
               plot( X, rstudent, pch = 20, xlab = expression( italic( "x" ) ),
                   ylab = "Studentized residuals", col = "blue", cex = 0.4,
                   cex.main = 0.9 )
               grid( NULL,NULL, lty = 6, col = "cornsilk2" )
               abline(h = 2, col= "red", lwd = 2, lty=2) # Reference line y = x
               abline(h = -2, col= "red", lwd = 2, lty=2) # Reference line y = x

               par( mar = c( 3.5, 2, 1, 1 ) + 0.07, mgp = c( 1.2, 0.4, 0 ) )
               qqnorm(rstudent,  col="blue" , pch = 20, cex=0.4, cex.main=0.9)
               grid( NULL,NULL, lty = 6, col = "cornsilk2" )
               qqline( rstudent, col= "red", lwd = 2 ) # Reference line y = x
               j = j + 1
           }
           names(fitLIST) <- distNAMES
           fitLIST = fitLIST[ as.character( aicDAT$Distribution ) ]
           cat( "** Done ***\n" )
           return(list(aic = aicDAT, fit = fitLIST, bestfit = bestFIT,
                   fitted=YbestFIT, rstudent=rstbestFIT, residuals=rbestFIT))
       } else {
           names(fitLIST) <- distNAMES
           fitLIST = fitLIST[ as.character( aicDAT$Distribution ) ]
           cat("** Done ***\n")
           return( list( aic = aicDAT, fit = fitLIST, bestfit = bestFIT ) )
       }
   }
}