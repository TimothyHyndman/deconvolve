#' Bandwidth Selectors for Deconvolution Kernel Density Estimation
#' 
#' Computes a bandwidth for use in deconvolution kernel density estimation of 
#' \eqn{X} from data \eqn{W = X + U}. If 'SIMEX' algorithm used, computes a 
#' bandwidth for use in deconvolution regression of data \eqn{(W, Y)}  where
#' \eqn{Y = g(X) + V} and \eqn{W = X + U}.
#' 
#' The function \code{bandwidth} chooses from one of four different methods 
#' depending on how the error distribution is defined and which algorithm is
#' selected.
#' 
#' \strong{PI for Homoscedastic Error:} If \code{algorithm = "PI"} and the errors 
#' are defined by either a single function \code{phiU}, or a single value 
#' \code{sigU} along with its \code{errortype}, then the method used is as
#' described in Delaigle and Gijbels 2002, and Delaigle and Gijbels 2004.
#' 
#' \strong{PI for Heteroscedastic Error:} If \code{algorithm = "PI"} and the
#' errors are defined by a either a vector of functions \code{phiU}, or a vector 
#' \code{sigU} along with its \code{errortype} then the method used is as 
#' described in Delaigle and Meister 2008.
#' 
#' \strong{CV:} If \code{algorithm = "CV"} then the method used is as described 
#' in Stefanski and Carroll 1990, and Delaigle and Gijbels 2004.
#' 
#' \strong{SIMEX:} If \code{algorithm = "SIMEX"} then the method used is as 
#' described in Delaigle and Hall 2008.
#' 
#' @inheritParams deconvolve
#' @param varX An estimate of the variance of \eqn{X}. Only required for 
#' heteroscedastic errors.
#' @param algorithm Either \code{"PI"} for plug-in estimator or \code{"CV"} for 
#' cross-validation estimator. If \code{"CV"} then the errors must be 
#' homoscedastic.
#' @param Y A vector of the univariate dependant data. Only required for 'SIMEX'
#' algorithm.
#' 
#' @return The bandwidth estimator. If using 'SIMEX' algorithm then returns a
#' list containing the bandwidth 'h' and ridge parameter 'rho'.
#' 
#' @section Warnings:
#' \itemize{
#' 	\item The arguments \code{phiK}, \code{muK2}, \code{RK}, and \code{tt} must
#' 	all be calculated from the same kernel. If you change one of these, you must
#' 	also change the rest to match.
#' }
#' 
#' @section References:
#' Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic 
#' error. \emph{Bernoulli}, 14, 2, 562-579.
#' 
#' Delaigle, A. and Gijbels, I. (2002). Estimation of integrated squared density 
#' derivatives from a contaminated sample. \emph{Journal of the Royal 
#' Statistical Society, B}, 64, 4, 869-886.
#'
#' Delaigle, A. and Gijbels, I. (2004). Practical bandwidth selection in 
#' deconvolution kernel density estimation. \emph{Computational Statistics and 
#' Data Analysis}, 45, 2, 249 - 267.
#' 
#' Stefanski, L. and Carroll, R.J. (1990). Deconvoluting kernel density 
#' estimators. \emph{Statistics}, 21, 2, 169-184.
#' 
#' Delaigle, A. and Hall, P. (2008). Using SIMEX for smoothing-parameter choice 
#' in errors-in-variables problems. \emph{Journal of the American Statistical 
#' Association}, 103, 481, 280-287 
#' 
#' @author Aurore Delaigle, Timothy Hyndman, Tianying Wang
#' 
#' @example man/examples/bandwidth_eg.R
#' 
#' @export

bandwidth <- function(W, errortype, sigU, phiU, Y = NULL, varX = NULL, 
					  algorithm = "PI", phiK = NULL, muK2 = 6, 
					  RK = 1024 / 3003 / pi, tt = seq(-1, 1, 2e-04)){
	
	n <- length(W)
	deltat <- tt[2] - tt[1]

	if(is.null(phiK)){
		phiK <- phiK2
	}

	# Check inputs 1 -----------------------------------------------------------

	if (missing(errortype) & missing(sigU) & missing(phiU)) {
		stop("You must define the error distribution.")
	}

	if (missing(phiU)) {
		if (missing(errortype)) {
			stop("You must supply the error type.")
		}

		if(missing(sigU)) {
				stop("You must provide the standard deviation of the errors.")
		}

		if ((errortype == "norm" | errortype == "Lap") == FALSE) {
			stop("errortype must be one of: 'norm', or 'Lap'.")
		}
	}

	if ((algorithm == "CV" | algorithm == "PI" | algorithm == "SIMEX") == FALSE) {
		stop("algorithm must be one of: 'PI', 'CV', or 'SIMEX'.")
	}

	# Determine Error Type Provided --------------------------------------------

	if (missing(phiU)) {
		if (length(sigU) > 1){
			errors <- "het"
			if ((length(sigU) == length(W)) == FALSE) {
				stop("sigU must be either length 1 for homoscedastic errors or 
					 have the same length as W for heteroscedastic errors.")
			}
		} else {
			errors <- "hom"
		}
	} else {
		if (length(phiU) > 1){
			errors <- "het"
			if ((length(phiU) == length(W)) == FALSE) {
				stop("phiU must be either length 1 for homoscedastic errors or 
					 have the same length as W for heteroscedastic errors.")
			}
		} else {
			errors <- "hom"
		}
	}

	# Check inputs 2 -----------------------------------------------------------

	if (algorithm == "CV") {
		if (errors == "het") {
			stop("Algorithm type 'CV' can only be used with homoscedastic 
				 errors.")
		}
	}

	if (is.null(varX) & errors == "het") {
		stop("You must supply an estimate for the variance of X when the errors are heteroscedastic.")
	}

	if (algorithm == "SIMEX") {
		if (is.null(Y)) {
			stop("You must supply Y to use SIMEX.")
		}
		if (errors == "het") {
			stop("Algorithm type 'SIMEX' can only be used with homoscedastic 
				 errors.")
		}
		if (missing(sigU)) {
			stop("Algorithm 'SIMEX' currently doesn't work with errors supplied
				 using phiU")
		}
	}

	# Convert errortype to phiU ------------------------------------------------

	if(missing(phiU)) {
		if(errortype == 'Lap' & errors == "hom") {
			phiU <- function(tt) {
				1 / (1 + sigU^2 * tt^2 / 2)
			}
		}

		if(errortype == 'norm' & errors == "hom") {
			phiU <- function(tt) {
				exp(-sigU^2 * tt^2 / 2)
			}
		}

		if(errortype == 'Lap' & errors == "het") {
			phiU <- c()
			for (sigUk in sigU){
				phiUk <- function(tt) {
					1 / (1 + sigUk^2 * tt^2 / 2)
				}
				phiU <- c(phiU, phiUk)
			}
		}

		if(errortype == 'norm' & errors == "het") {
			phiU <- c()
			for (sigUk in sigU){
				phiUk <- function(tt) {
					exp(-sigUk^2 * tt^2 / 2)
				}
				phiU <- c(phiU, phiUk)
			}
		}
	}

	# Perform appropriate bandwidth calculation --------------------------------

	if (algorithm == "CV"){
		output <- CVdeconv(n, W, errortype, sigU, phiU, phiK, muK2, RK, deltat, 
						   tt)
	}

	if (algorithm == "PI" & errors == "het") {
		output <- PI_deconvUknownth4het(n, W, varX, phiU, phiK, muK2, RK, 
										deltat, tt)
	}

	if (algorithm == "PI" & errors == "hom") {
		output <- PI_deconvUknownth4(n, W, phiU, phiK, muK2, RK, deltat, tt)
	}

	if (algorithm == "SIMEX") {
		output <- hSIMEXUknown(W, Y, errortype, sigU)
	}

	output
}