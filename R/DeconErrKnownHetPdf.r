# #' Deconvolution when the errors are heteroscedastic with known distributions
# #' 
# #' Computes the deconvolution kernel density estimator (KDE) of \eqn{X} from 
# #' data \eqn{W_j = X + U_j} when the distribution and standard deviations of the 
# #' \eqn{U_j} are known.
# #' 
# #' PUT DETAILS HERE
# #' 
# #' @inheritParams deconvolve
# #' @inheritParams DeconErrKnownPdf
# #' @param errortype The distribution type of the \eqn{U_j}. Either "Lap" for 
# #' Laplace errors or "norm" for normal errors. If you use this way of defining 
# #' the errors then you must also provide \code{sigUj} but should not provide
# #' \code{phiUkvec}.
# #' @param sigUj	A vector of length n which contains the standard deviations of 
# #' the errors.
# #' @param phiUkvec A vector of n functions that give the characteristic 
# #' functions of the errors. Produce this vector by c(func1,func2,...,funcn) 
# #' where each funcj is a function of tt. If you define the errors this way then 
# #' you should not provide \code{errortype} or \code{sigUj}.
# #' 
# #' @inheritSection deconvolve Warnings
# #' 
# #' @inherit deconvolve return
# #' 
# #' @section References:
# #' Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic 
# #' error. \emph{Bernoulli}, 14, 2, 562-579.
# #' 
# #' @section Author:
# #' Aurore Delaigle
# #' 
# #' @example man/examples/KnownHetError_eg.R
# #' 
# #' @export

DeconErrKnownHetPdf<-function(xx, W, h, errortype, sigUj, phiUkvec, 
	rescale = FALSE, phiK = NULL, muK2 = 6, RK = 1024 / 3003 / pi, 
	tt = seq(-1, 1, 2e-04)){

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
			phiUkvec[[k]](tt)
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
