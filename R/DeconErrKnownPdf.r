# #' Deconvolution when the error distribution is known
# #' 
# #' Computes the deconvolution kernel density estimator (KDE) of \eqn{X} from 
# #' data \eqn{W = X + U} when the distribution of \eqn{U} is known.
# #' 
# #' PUT DETAILS HERE
# #' 
# #' @inheritParams deconvolve
# #' @param h The bandwidth to use.
# #' @param errortype The distribution type of \eqn{U}. Either "Lap" for Laplace 
# #' errors or "norm" for normal errors. If you define the errors this way then 
# #' you must also provide \code{sigU} but should not provide \code{phiU}.
# #' @param sigU The standard deviation of \eqn{U}.
# #' @param phiU A function giving the characteristic function of \eqn{U}. If you 
# #' define the errors this way then you should not provide \code{errortype} or 
# #' \code{sigU}.
# #' 
# #' @inherit deconvolve return
# #' 
# #' @inheritSection deconvolve Warnings
# #' 
# #' @section Author:
# #' Aurore Delaigle
# #' 
# #' @section References:
# #' 
# #' @example man/examples/KnownError_eg.R
# #' 
# #' @export

DeconErrKnownPdf<-function(xx, W, h, phiU, rescale = FALSE, 
	phiK = NULL, muK2 = 6, RK = 1024 / 3003 / pi, tt = seq(-1, 1, 2e-04)){

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