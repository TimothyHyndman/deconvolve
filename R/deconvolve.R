#' Deconvolution Kernel Density Estimator
#' 
#' Computes the deconvolution kernel density estimator (KDE) of \eqn{X} from 
#' data \eqn{W = X + U} when the distribution of \eqn{U} is known,
#' heteroscedastic, or symmetric.
#' 
#' The function \code{deconvolve} chooses from one of three different methods 
#' depending on how the error distribution is defined.
#' 
#' \strong{Symmetric Error:} If neither \code{errortype} and \code{sigU}, or 
#' \code{phiU} are supplied then the error is assumed symmetric and the 
#' deconvolution method is based on the method described in Delaigle and Hall 
#' 2016.
#' 
#' \strong{Homoscedastic Error:} If the errors are defined by either a single 
#' function \code{phiU}, or a single value \code{sigU} along with its 
#' \code{errortype} then the method used is as described in Stefanski and
#' Carroll 1990.
#' 
#' \strong{Heteroscedastic Errors:} If the errors are defined by a either a 
#' vector of functions \code{phiU}, or a vector \code{sigU} along with its 
#' \code{errortype} then the method used is as described in Delaigle and 
#' Meister 2008.
#' 
#' Errors can be defined by either a distribution type (\code{errortype}) along 
#' with the standard deviation(s) (\code{sigU}), or by the characteric 
#' function(s) of the errors (\code{phiU}). 
#' 
#' @param W A vector of the univariate contaminated data.
#' @param xx A vector of x values on which to compute the density.
#' @param errortype The distribution type of \eqn{U}. Either "Lap" for Laplace 
#' errors or "norm" for normal errors. If you define the errors this way then 
#' you must also provide \code{sigU} but should not provide \code{phiU}.
#' @param sigU The standard deviations of \eqn{U}. A single value for
#' homoscedastic errors and a vector having the same length as \code{W} for 
#' heteroscedastic errors.
#' @param phiU A function giving the characteristic function of \eqn{U}. A 
#' single value for homoscedastic errors and a vector having the same length as 
#' \code{W} for heteroscedastic errors. If you define the errors this way then 
#' you should not provide \code{errortype} or \code{sigU}.
#' @param bw The bandwidth to use. If \code{NULL}, a bandwidth will be
#' calculated using an appropriate plug-in estimator.
#' @param varX An estimate of the variance of \eqn{X}. Only required when 
#' \code{bw = NULL} and the errors are heteroscedastic.
#' @param rescale If \code{TRUE}, estimator is rescaled so that it 
#' integrates to 1. Rescaling requires \code{xx} to be a fine grid of equispaced 
#' \eqn{x} values that covers the whole range of \eqn{x}-values where the 
#' estimated density is significantly non zero.
#' @param pmf If \code{TRUE}, returns a probability mass function instead of a 
#' density as the estimator. This is quicker than estimating a density. To use
#' this option, the errors must not be provided.
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
#' @return An object of class "\code{deconvolve}".
#' 
#' The function \code{plot} produces a plot of the deconvolution KDE.
#' 
#' An object of class "\code{deconvolve}" is a list containing at least some of
#' the elements:
#' \item{W}{The original contaminated data}
#' \item{x}{The values on which the deconvolution KDE is evaluated.}
#' \item{pdf}{A vector containing the deconvolution KDE evaluated at each point 
#' in \code{x}}
#' \item{support}{The support of the pmf found when the errors are assumed
#' symmetric}
#' \item{probweights}{The probability masses of the pmf found when the errors
#' are assumed symmetric}
#' 
#' @section Warnings:
#' \itemize{
#'	\item The arguments \code{phiK}, \code{muK2}, \code{RK}, and \code{tt} must
#' 	all be calculated from the same kernel. If you change one of these, you must
#' 	also change the rest to match.
#'	\item If you supply your own bandwidth, then you should ensure that the
#' 	kernel used here matches the one you used to calculate your bandwidth.
#'	\item The DKDE can also be computed using the Fast Fourier Transform, which 
#' 	is a bit more complex. See Delaigle and Gijbels 2007. However if the grid of 
#' 	t-values is fine enough, the estimator can simply be computed like here 
#' 	without having problems with oscillations.
#' }
#' 
#' @section References:
#' Stefanski, L.A. and Carroll, R.J. (1990). Deconvolving kernel density
#' estimators. \emph{Statistics}, 21, 2, 169-184.
#' 
#' Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic 
#' error. \emph{Bernoulli}, 14, 2, 562-579.
#' 
#' Delaigle, A. and Hall, P. (2016). Methodology for non-parametric 
#' deconvolution when the error distribution is unknown. \emph{Journal of the 
#' Royal Statistical Society: Series B (Statistical Methodology)}, 78, 1, 
#' 231-252.
#' 
#' Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating 
#' integrals and optimizing objective functions: a case study in density 
#' deconvolution. \emph{Statistics and Computing}, 17, 349-355.
#' 
#' @author Aurore Delaigle, Timothy Hyndman, Tianying Wang
#' 
#' @example man/examples/deconvolve_eg.R
#' 
#' @export

deconvolve <- function(W, xx, errortype = NULL, sigU = NULL, phiU = NULL, 
					   bw = NULL, varX = NULL, rescale = FALSE, pmf = FALSE, 
					   phiK = NULL, muK2 = 6, RK = 1024 / 3003 / pi, 
					   tt = seq(-1, 1, 2e-04)){

	if(is.null(phiK)){
		phiK <- phiK2
	}

	# Decide on type of deconvolution ------------------------------------------
	if (is.null(errortype) & is.null(phiU)) {
		decon_type <- "symmetric"
	} else if (length(sigU) > 1  | length(phiU) > 1){
		decon_type <- "heteroscedastic"
	} else {
		decon_type <- "known"
	}

	if (is.null(errortype) == FALSE) {
		if ((errortype == "norm" | errortype == "Lap") == FALSE) {
			stop("errortype must be one of: 'norm', or 'Lap'.")
		}
	}

	# Calculate Bandwidth if not supplied --------------------------------------
	if (is.null(bw) & (decon_type == "symmetric") == FALSE) {
		if (is.null(phiU)) {
			bw <- bandwidth(W, errortype, sigU, varX = varX, phiK = phiK, 
							muK2 = muK2, RK = RK, tt = tt)
		} else {			
			bw <- bandwidth(W, phiU = phiU, varX = varX, phiK = phiK, 
							muK2 = muK2, RK = RK, tt = tt)
		}
	}

	# Check inputs -------------------------------------------------------------
	if (pmf & !(decon_type == "symmetric")){
		stop("Option pmf cannot be used when the error is provided.")
	}

	# Convert errortype to phiU ------------------------------------------------

	if (!(decon_type == "symmetric")){
		if(is.null(phiU)) {
			if(errortype == 'Lap' & decon_type == "known") {
				phiU <- function(tt) {
					1 / (1 + sigU^2 * tt^2 / 2)
				}
			}

			if(errortype == 'norm' & decon_type == "known") {
				phiU <- function(tt) {
					exp(-sigU^2 * tt^2 / 2)
				}
			}

			if(errortype == 'Lap' & decon_type == "heteroscedastic") {
				phiU <- c()
				for (sigUk in sigU){
					phiUk <- function(tt) {
						1 / (1 + sigUk^2 * tt^2 / 2)
					}
					phiU <- c(phiU, phiUk)
				}
			}

			if(errortype == 'norm' & decon_type == "heteroscedastic") {
				phiU <- c()
				for (sigUk in sigU){
					phiUk <- function(tt) {
						exp(-sigUk^2 * tt^2 / 2)
					}
					phiU <- c(phiU, phiUk)
				}
			}
		}
	}

	# Perform appropriate deconvolution ----------------------------------------
	if (decon_type == "known"){
		pdf <- DeconErrKnownPdf(xx, W, bw, phiU, rescale, phiK, muK2, RK, tt)
		output <- list("x" = xx, "pdf" = pdf, "W" = W)
	}

	if (decon_type == "heteroscedastic"){
		pdf <- DeconErrKnownHetPdf(xx, W, bw, phiU, rescale, phiK, muK2, RK, tt)
		output <- list("x" = xx, "pdf" = pdf, "W" = W)
	}

	if (decon_type == "symmetric") {
		out <- DeconErrSymPmf(W)
		if (!pmf) {
			phi.W <- out$phi.W
			pdf <- DeconErrSymPmfToPdf(out, W, phi.W, xx, phiK, muK2, tt, 
										  rescale, bw)
			output <- list("x" = xx, "pdf" = pdf, "support" = out$support, 
						   "probweights" = out$probweights, "W" = W)
		} else {
			output <- list("support" = out$support,
						   "probweights" = out$probweights,
						   "W" = W)
		}
	}

	# Output object of class "deconvolve" --------------------------------------
	class(output) <- "deconvolve"
	output
}
