#' Deconvolution when the errors are heteroscedastic with known distributions
#' 
#' Computes the deconvolution kernel density estimator (KDE) of \eqn{X} from 
#' data \eqn{W = X + U_j} when the distributions of the \eqn{U_j} are known.
#' 
#' PUT DETAILS HERE
#' 
#' @inheritParams DeconErrKnownPdf
#' @param errortype The distribution type of \eqn{U}. Either "Lap" for Laplace 
#' errors or "norm" for normal errors. If you use this way of defining the 
#' errors then you must also provide \code{sigUj}.
#' @param sigUj	A vector of length n which contains standard deviations of each 
#' of the errors.
#' @param phiUkvec A vector of n functions that give the characteristic 
#' functions of the errors. Produce this vector by c(func1,func2,...,funcn) 
#' where each funcj is a function of tt. If you define the errors this way then 
#' you should not provide \code{errortype} or \code{sigUj}.
#' 
#' @section Warnings:
#' \enumerate{
#'	\item Rescaling requires xx to be a fine grid of equispaced x-values that 
#'	covers the whole range of x-values where the estimated density is 
#'	significantly non zero.
#'	\item Changing the kernel: if you change one of the arguments among phiK, 
#'	muK2, RK, deltat and tt, you must change them all as they need to correspond 
#'	to the same kernel.
#'	\item	If phiK is compactly supported, the first and last elements of t 
#'	must be the lower and upper bound of the support of phiK.
#'	\item	If phiK is not compactly supported, the first and last elements of t 
#'	must be larger enough for your discretisation of the intergals to be 
#'	accurate
#'	\item The kernel K here must match the phiK used to compute the bandwidth 
#'	(PI, CV or other)
#'	\item The DKDE can also be computed using the Fast Fourier Transform, which 
#' 	is a bit more complex. See Delaigle, A. and Gijbels, I. (2007). Frequent 
#' 	problems in calculating integrals and optimizing objective functions: a case 
#' 	study in density deconvolution. Statistics and Computing, 17, 349-355
#'	\item However if the grid of t-values is fine enough, the estimator can 
#' 	simply be computed like here without having problems with oscillations.
#'	}
#' 
#' @return The outcome is the deconvolution kernel density estimator when the 
#' errors are heteroscedastic.
#' 
#' @section References:
#' Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic 
#' error. Bernoulli, 14, 562-579
#' 
#' @section Author:
#' Aurore Delaigle
#' 
#' @example man/examples/KnownHetError_eg.R
#' 
#' @keywords heteroscedastic deconvolution kernel density estimator
#' 
#' @export

DeconErrKnownHetPdf<-function(xx, W, h, errortype, sigUj, phiUkvec, rescale=0, 
	phiK = NULL, muK2 = 6, RK = 1024 / 3003 / pi, tt = seq(-1, 1, 2e-04)){

	#--------------------------------------------------------------------------#
	# Check optional arguments
	#--------------------------------------------------------------------------#
	if(missing(errortype) & missing(phiUkvec)){
		stop("You must define the error distributions")
	}

	if(missing(errortype) == F){
		if(missing(sigUj)){
			stop("You must provide the standard deviations of the errors")
		}
		if(errortype == 'Lap'){
			phiUk <- function(tt, k){
				1 / (1 + sigUj[k]^2 * tt^2 / 2)
			}
		} else if(errortype == 'norm'){
			phiUk <- function(tt, k){
				exp(-sigUj[k]^2 * tt^2 / 2)
			}
		} else {
			stop("Error type not a valid option")
		}
	}

	# Characteristic function of kth error 
	if(missing(errortype) & (missing(phiUkvec) == F)){
		phiUk <- function(tt, k){
			phiUkvec[[k]](tt, k)
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
	dim(tt) <- c(length(tt), 1)

	# Default values of phiU(t)=characteristic function of the errors
	# If you want to consider another error type, simply replace phiU by the 
	# characteristic function of your error type

	OO <- outerop(tt / h, t(W), "*")

	# Compute phiU_k(-t/h) for each k -- since phiU_k is symmetric, this is the 
	# same as phiU_k(t/h)
	matphiU <- OO;
	for (k in 1:n)
		matphiU[, k] <- phiUk(tt / h, k)

	# Sum by rows of the matrix. This produces a vector of size equal to that of 
	# tt
	phiUsqth <- apply(matphiU^2, 1, sum)

	# Estimate real and imaginary parts of empirical characteristic function of 
	# W computed at tt/h, for each component of tt.
	# Results = vectors of size length(tt)

	rehatphiX <- apply(cos(OO) * matphiU, 1, sum) / phiUsqth
	imhatphiX <- apply(sin(OO) * matphiU, 1, sum) / phiUsqth

	# Matrix of size length(tt) x length(xx)
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
