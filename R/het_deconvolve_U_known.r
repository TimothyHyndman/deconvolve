#' Deconvolution Kernel Density Estimator when the errors are known and 
#' heteroscedastic
#' 
#' Computes the deconvolution kernel density estimator (KDE) of \eqn{X} from 
#' data \eqn{W = X + U} when the distributions of \eqn{U} are known and 
#' heteroscedastic.
#' 
#' The method used is as described in Delaigle and Meister 2008.
#' 
#' @param W A vector of the univariate contaminated data.
#' @param x A vector of x values on which to compute the density.
#' @param phiUkvec A vector of functions giving the characteristic functions of 
#' \eqn{U}. Must have the same length as \code{W}.
#' @param h The bandwidth to use.
#' @param rescale If \code{TRUE}, estimator is rescaled so that it 
#' integrates to 1. Rescaling requires \code{x} to be a fine grid of equispaced 
#' \eqn{x} values that covers the whole range of \eqn{x}-values where the 
#' estimated density is significantly non zero.
#' @param kernel_type The deconvolution kernel to use. The default kernel has
#' characteristic function \eqn{(1-t^2)^3}.
#'  
#' @return An object of class "\code{deconvolve}" containing the elements
#' \item{W}{The original contaminated data}
#' \item{x}{The values on which the deconvolution KDE is evaluated.}
#' \item{pdf}{A vector containing the deconvolution KDE evaluated at each point 
#' in \code{x}}
#' 
#' The function \code{plot} produces a plot of the deconvolution KDE.
#'
#' @section Warnings:
#' \itemize{
#'	\item You should ensure that the kernel used here matches the one you used 
#' 	to calculate your bandwidth.
#' \item The DKDE can also be computed using the Fast Fourier Transform, which 
#' 	is a bit more complex. See Delaigle and Gijbels 2007. However if the grid of 
#' 	t-values is fine enough, the estimator can simply be computed like here 
#' 	without having problems with oscillations.
#' }
#' 
#' @section References:
#' Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic 
#' error. \emph{Bernoulli}, 14, 2, 562-579.
#' 
#' Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating 
#' integrals and optimizing objective functions: a case study in density 
#' deconvolution. \emph{Statistics and Computing}, 17, 349-355.
#' 
#' @author Aurore Delaigle, Timothy Hyndman, Tianying Wang
#' 
#' @example man/examples/het_deconvolve_U_known_eg.R
#' 
#' @export

het_deconvolve_U_known <- function(W,
								   phiUkvec,
								   h,
								   x = seq(min(W), max(W), length.out = 100),
								   kernel_type = c("default", "normal", "sinc"),
								   rescale = FALSE) {

	kernel_type <- match.arg(kernel_type)

	if ((length(phiUkvec) == length(W)) == FALSE) {
		stop("phiUkvec must have the same length as W")
	}

	if (kernel_type == "normal") {
		warning("You should only use the 'normal' kernel when the errors are 
			Laplace or convolutions of Laplace.")
	}

	if (kernel_type == "sinc") {
		warning("You should ensure that you are not using a plug-in bandwidth 
			method for the bandwidth when using the sinc kernel.")
	}

	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]

	n <- length(W)

	phiUk <- function(tt,k) {
		phiUkvec[[k]](tt)
	}

	OO <- outer(tt / h, W)

	# Compute phiU_k(-t/h) for each k -- since phiU_k is symmetric, this is the 
	# same as phiU_k(t/h)
	matphiU <- OO
	for (k in 1:n) {
		matphiU[, k] <- phiUk(tt / h, k)
	}

	phiUsqth <- apply(matphiU^2, 1, sum)

	# Estimate real and imaginary parts of empirical characteristic function of 
	# W computed at tt/h, for each component of tt.
	rehatphiX <- apply(cos(OO) * matphiU, 1, sum) / phiUsqth
	imhatphiX <- apply(sin(OO) * matphiU, 1, sum) / phiUsqth

	xt <- outer(tt / h, x)
	longx <- length(x)

	# Compute the DKDE estimator
	fXdecUK <- cos(xt) * kronecker(matrix(1, 1, longx), rehatphiX) + sin(xt) * 
			   kronecker(matrix(1, 1, longx), imhatphiX)
	fXdecUK <- apply(fXdecUK * kronecker(matrix(1, 1, longx), phiK(tt)), 2, sum) / 
			   (2 * pi) * deltat / h
	fXdecUK[which(fXdecUK < 0)] <- 0

	if (rescale == 1){
		dx <- x[2] - x[1]
		fXdecUK <- fXdecUK / sum(fXdecUK) / dx
	}

	# Construct final output object --------------------------------------------
	output <- list("x" = x, "pdf" = fXdecUK, "W" = W)
	class(output) <- c("deconvolve", "list")
	output
}
