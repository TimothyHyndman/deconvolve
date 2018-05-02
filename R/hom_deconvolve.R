#' Deconvolution KDE when the error is unknown
#' 
#' Computes the deconvolution kernel density estimator (KDE) of \eqn{X} from 
#' data \eqn{W_1 = X + U_1} when the distribution of \eqn{U} is either unknown, 
#' or estimated from replicates, \eqn{W_2 = X + U_2}.
#' 
#' The function \code{hom_deconvolve} chooses from one of two different methods.
#' 
#' \strong{Error from Replicates:} If both \code{W1} and \code{W2} are supplied
#' then the error is calculated using replicates. This method was prototyped in 
#' Delaigle, Hall, and Meister 2008 and then further refined in Delaigle and 
#' Hall 2016, and Camirand, Carroll, and Delaigle 2018.
#' 
#' \strong{Symmetric Error:} If \code{W2} isn't supplied then the error is 
#' assumed symmetric and the deconvolution method is based on the method 
#' described in Delaigle and Hall 2016.
#' 
#' @param W1 A vector of the univariate contaminated data.
#' @param W2 A vector of replicate measurements. If supplied, then the error 
#' will be estimated using replicates.
#' @param x A vector of x values on which to compute the density.
#' @param h The bandwidth to use. If \code{NULL}, a bandwidth will be
#' calculated using an appropriate plug-in estimator.
#' @param rescale If \code{TRUE}, estimator is rescaled so that it 
#' integrates to 1. Rescaling requires \code{x} to be a fine grid of equispaced 
#' \eqn{x} values that covers the whole range of \eqn{x}-values where the 
#' estimated density is significantly non zero.
#' @param kernel_type The deconvolution kernel to use. The default kernel has
#' characteristic function \eqn{(1-t^2)^3}.
#' @param pmf Only used when \code{W2 = NULL}. If \code{TRUE}, returns a 
#' probability mass function instead of a density as the estimator. This is 
#' quicker than estimating a density.
#' @param m Only used when \code{W2 = NULL}. The number of point masses to use 
#' to estimate the distribution of \eqn{X}.
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
#' 	\item The method for deconvolution when the error is unknown and assumed
#' 	symmetric (as described in Delaigle and Hall (2016)) requires solving a 
#' 	non-linear objective function with both linear and non-linear constraints. 
#' 	We are yet to find a package in R that can perform this reliably. We instead
#' 	recommend using the MATLAB code found at <URL> as it is both faster, and 
#' 	more reliable.
#'	\item If you supply your own bandwidth, then you should ensure that the
#' 	kernel used here matches the one you used to calculate your bandwidth.
#'	\item The DKDE can also be computed using the Fast Fourier Transform, which 
#' 	is a bit more complex. See Delaigle and Gijbels 2007. However if the grid of 
#' 	t-values is fine enough, the estimator can simply be computed like here 
#' 	without having problems with oscillations.
#' }
#' 
#' @section References: 
#' Delaigle, A., Hall, P., and Meister, A. (2008). On Deconvolution with 
#' repeated measurements. \emph{Annals of Statistics}, 36, 665-685
#' 
#' Delaigle, A. and Hall, P. (2016). Methodology for non-parametric 
#' deconvolution when the error distribution is unknown. \emph{Journal of the 
#' Royal Statistical Society: Series B (Statistical Methodology)}, 78, 1, 
#' 231-252.
#' 
#' Camirand, F., Carroll, R.J., and Delaigle, A. (2018). Estimating the 
#' distribution of episodically consumed food measured with errors. 
#' \emph{Manuscript.}
#' 
#' Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating 
#' integrals and optimizing objective functions: a case study in density 
#' deconvolution. \emph{Statistics and Computing}, 17, 349-355.
#' 
#' @author Aurore Delaigle, Timothy Hyndman, Tianying Wang
#' 
#' @example man/examples/hom_deconvolve_eg.R
#' 
#' @export

hom_deconvolve <- function(W1, 
						   W2 = NULL,
						   x = seq(min(W1), max(W1), length.out = 100),
						   h = NULL,
						   kernel_type = c("default", "normal", "sinc"),
						   rescale = FALSE,
						   pmf = FALSE,
						   m = 20) {

	kernel_type <- match.arg(kernel_type)

	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	muK2 <- kernel_list$muk2
	RK <- kernel_list$rk
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]

	n <- length(W1)

	if (!is.null(W2)) {
		errors <- "rep"
	} else {
		errors <- "sym"
		warning("The method for deconvolution when the error is unknown and assumed symmetric is slow and unreliable in R. Consider instead using the MATLAB code found at <URL>.")
	}

	if (is.null(h) & !(errors == "sym")) {
			h <- bandwidth(W1, W2, kernel_type = kernel_type)
	}

	if (errors == "rep") {
		t_search <- tt/h
		phiU_splined <- function(t){
			replicates_phiU(t, W1, W2, t_search)
		}
		W <- c(W1, W2)
		output <- hom_deconvolve_U_known(W, phiU_splined, h, x, kernel_type, rescale)
		output <- list("x" = x, "pdf" = output$pdf, "W1" = W1, "W2" = W2)
	}

	if (errors == "sym") {
		W <- W1
		out <- DeconErrSymPmf(W, m, kernel_type)
		h <- NULL

		if (!pmf) {
			phi_W <- out$phi.W
			pdf <- DeconErrSymPmfToPdf(out, W, phi_W, x, kernel_type, rescale, h)
			output <- list("x" = x, "pdf" = pdf, "support" = out$support, 
						   "probweights" = out$probweights, "W" = W)
		} else {
			output <- list("support" = out$support,
						   "probweights" = out$probweights,
						   "W" = W)
		}
	}


	class(output) <- c("deconvolve", "list")
	output
}