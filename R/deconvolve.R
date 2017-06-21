#' Deconvolution Kernel Density Estimator
#' 
#' Computes the deconvolution kernel density estimator (KDE) of \eqn{X} from 
#' data \eqn{W = X + U} when the distribution of \eqn{U} is known,
#' heteroscedastic, or symmetric.
#' 
#' PUT DETAILS HERE
#' 
#' @param W A vector of the univariate contaminated data.
#' @param xx A vector of x values on which to compute the density.
#' @param errortype The distribution type of \eqn{U}. Either "Lap" for Laplace 
#' errors or "norm" for normal errors. If you define the errors this way then 
#' you must also provide \code{sigU} but should not provide \code{phiU}.
#' @param sigU The standard deviations of \eqn{U}. 
#' @param phiU A function giving the characteristic function of \eqn{U}. If you 
#' define the errors this way then you should not provide \code{errortype} or 
#' \code{sigU}.
#' @param bw The bandwidth to use. If \code{NULL}, a bandwidth will be
#' calculated using an appropriate plug-in estimator.
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
#' \code{xx}.
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
#' @section References:
#' 
#' @example man/examples/deconvolve_eg.R
#' 
#' @export

deconvolve <- function(W, xx, errortype = NULL, sigU = NULL, phiU = NULL, 
					   bw = NULL, rescale = FALSE, phiK = NULL, muK2 = 6, 
					   RK = 1024 / 3003 / pi, tt = seq(-1, 1, 2e-04)){

	# Decide on type of deconvolution ------------------------------------------
	if (is.null(errortype) & is.null(phiU)) {
		decon_type <- "symmetric"
	} else if (length(sigU) > 1  | length(phiU) > 1){
		decon_type <- "heteroscedastic"
	} else {
		decon_type <- "known"
	}

	# Calculate Bandwidth if not supplied --------------------------------------
	if (is.null(bw) & (decon_type == "symmetric") == FALSE) {
		if (is.null(phiU)) {
			bw <- bandwidth(W, errortype, sigU, phiK = phiK, muK2 = muK2, 
							RK = RK, tt = tt)
		} else {
			bw <- bandwidth(W, phiU = phiU, phiK = phiK, muK2 = muK2, 
							RK = RK, tt = tt)
		}
	}
	
	# Use default PhiK if not supplied -----------------------------------------
	if(is.null(phiK)){
		phiK <- phiK2
	}

	# Perform appropriate deconvolution ----------------------------------------
	if (decon_type == "known"){
		output <- DeconErrKnownPdf(xx, W, bw, errortype, sigU, phiU, rescale, 
								   phiK, muK2, RK, tt)
	}

	if (decon_type == "heteroscedastic"){
		output <- DeconErrKnownHetPdf(xx, W, bw, errortype, sigU, phiU, rescale, 
									  phiK, muK2, RK, tt)
	}

	if (decon_type == "symmetric") {
		out <- DeconErrSymPmf(W)
		phi.W <- out$phi.W
		output <- DeconErrSymPmfToPdf(out, W, phi.W, xx, phiK, muK2, tt, 
									  rescale, bw)
	}

	# Output PDF ---------------------------------------------------------------
	output
}
