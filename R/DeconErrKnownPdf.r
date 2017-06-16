#' Deconvolution when the error distribution is known
#' 
#' Computes the deconvolution kernel density estimator (KDE) of \eqn{X} from 
#' data \eqn{W = X + U} when the distribution of \eqn{U} is known.
#' 
#' PUT DETAILS HERE
#' 
#' @param xx A vector of x values on which to compute the density.
#' @param W A vector of the univariate contaminated data.
#' @param h The bandwidth to use.
#' @param errortype The distribution type of \eqn{U}. Either "Lap" for Laplace 
#' errors or "norm" for normal errors. If you use this way of defining the error 
#' then you must also provide \code{sigU} but should not provide \code{phiU}.
#' @param sigU The standard deviation of \eqn{U}.
#' @param phiU A function giving the characteristic function of \eqn{U}. If you 
#' define the errors this way then you should not provide \code{errortype} or 
#' \code{sigU}.
#' @param rescale If \code{TRUE}, estimator is rescaled so that it 
#' integrates to 1. Rescaling requires xx to be a fine grid of equispaced x 
#' values that covers the whole range of x-values where the estimated density is 
#' significantly non zero.
#' @param phiK A function giving the fourier transform of the kernel. 
#' If supplied, \code{muK2}, \code{RK}, and \code{tt} must also be supplied. If 
#' not supplied it defaults to \eqn{(1 - t^2)^3} on the interval \eqn{[-1,1]}.
#' @param muK2 The second moment of the kernel, i.e. \eqn{\int x^2 K(x) dx}.
#' @param RK The integral of the square of the kernel, i.e. \eqn{\int K^2(x) dx}.
#' @param tt A vector of evenly spaced t values on which to approximate the 
#' integrals in the Fourier domain. If phiK is compactly supported, the first 
#' and last elements of \code{tt} must be the lower and upper bound of the 
#' support of phiK. If phiK is not compactly supported, the first and last 
#' elements of \code{tt} must be large enough for your discretisation of the 
#' integrals to be accurate.
#' 
#' @return A vector containing the deconvolution KDE evaluated at each point in 
#' xx.
#' 
#' @section Warnings:
#' \itemize{
#'	\item The arguments \code{phiK}, \code{muK2}, \code{RK}, and \code{tt} must
#' 	all be calculated from the same kernel. If you change one of these, you must
#' 	also change the rest to match.
#'	\item The kernel used here must match the kernel used to compute the 
#' 	bandwidth.
#'	\item The DKDE can also be computed using the Fast Fourier Transform, which 
#' 	is a bit more complex. See Delaigle, A. and Gijbels, I. (2007). Frequent 
#' 	problems in calculating integrals and optimizing objective functions: a case 
#' 	study in density deconvolution. \emph{Statistics and Computing}, 17, 
#' 	349-355. However if the grid of t-values is fine enough, the estimator can 
#' 	simply be computed like here without having problems with oscillations.
#' }
#' 
#' @section Author:
#' Aurore Delaigle
#' 
#' @section References:
#' 
#' @example man/examples/KnownError_eg.R
#' 
#' @export

DeconErrKnownPdf<-function(xx, W, h, errortype, sigU, phiU, rescale = FALSE, 
	phiK = NULL, muK2 = 6, RK = 1024 / 3003 / pi, tt = seq(-1, 1, 2e-04)){

	#--------------------------------------------------------------------------#
	# Check optional arguments
	#--------------------------------------------------------------------------#
	if(missing(errortype) & (missing(sigU)) & missing(phiU)){
		stop("You must define the error distribution")
	}

	if(missing(errortype) == F){
		if(missing(sigU)){
			stop("You must provide the standard deviation of the errors")
		}
		if(errortype == 'Lap'){
			phiU <- function(tt){
				1 / (1 + sigU^2 * tt^2 / 2)
			}
		} else if(errortype == 'norm'){
			phiU <- function(tt){
				exp(-sigU^2 * tt^2 / 2)
			} 
		} else {
			stop("Error type not a valid option")
		}
	}

	if(missing(h)){
		stop("You must provide the bandwidth")
	}

	#--------------------------------------------------------------------------#
	# Compute DKDE
	#--------------------------------------------------------------------------#
	if(is.null(phiK)){
		phiK <- phiK2
	}

	W <- as.vector(W)
	n <- length(W)
	deltat <- tt[2] - tt[1]

	# Make sure t is a vector in the right format
	dim(tt) <- c(length(tt), 1);

	OO <- outerop(tt/h, t(W), "*")
	phiUth <- phiU(tt/h)

	# Estimate real and imaginary parts of empirical characteristic function of 
	# W computed at tt/h, for each component of tt.
	rehatphiX <- rowSums(cos(OO)) / phiUth / n
	imhatphiX <- rowSums(sin(OO)) / phiUth / n

	xt <- outerop(tt / h, t(xx), "*")
	longx <- length(xx)

	# Compute the DKDE estimator
	fXdecUK <- cos(xt) * kronecker(matrix(1, 1, longx), rehatphiX) + sin(xt) *
			   kronecker(matrix(1, 1, longx), imhatphiX)
	fXdecUK <- apply(fXdecUK * kronecker(matrix(1, 1, longx), phiK(tt)), 2, sum) / 
			   (2 * pi) * deltat / h
	fXdecUK[which(fXdecUK < 0)] <- 0

	if (rescale == 1){
		dx <- xx[2] - xx[1]
		fXdecUK <- fXdecUK / sum(fXdecUK) / dx
	}


	return (fXdecUK)
}