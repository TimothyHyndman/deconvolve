#' Deconvolution Kernel Density Estimator
#' 
#' Computes the deconvolution kernel estimator (KDE) of the density of \eqn{X} from 
#' data \eqn{W_{i1} = X_i + U_{i1}, i=1,...,n} when the distribution of \eqn{U_i} is known 
#' or unknown and estimated from replicates, \eqn{W_{i1} = X_i + U_{i1}} and 
#' \eqn{W_{i2} = X_i + U_{i2}}, or without replicates if replicates are not available.
#' All error densities need to be symmetric. In the homoscedastic error case, the codes
#'  are suitable only if the charactersitic function of the errors is nonzero everywhere. 
#' In the heteroscedastic error case, the pooled characteristic function of the errors 
#' used by Delaigle and Meister (2006) must be nonzero everywhere
#' 
#' The function \code{deconvolve} chooses from one of five different methods 
#' depending on how the error distribution is defined/computed:
#' 
#' \strong{Known homoscedastic error distribution:} If the error distribution 
#' is defined by either a single function \code{phiU}, or a single value 
#' \code{sd_U} along with its \code{errortype} then the estimator of the density
#' of \eqn{X} is the one in Stefanski and Carroll (1990). 
#' 
#' \strong{Known heteroscedastic error distributions:} If the error distributions are 
#' defined by a either a vector of functions \code{phiU}, or a vector 
#' \code{sd_U} along with their \code{errortype} (one single error type for all individuals)
#' then the method used is the one from Delaigle and Meister (2008).
#' 
#' \strong{Unknown homoscedastic error distribution when replicates are available:} If both 
#' \code{W1} and \code{W2} are supplied and \code{het_replicates} is \code{FALSE}, then the 
#' error distribution is estimated using the replicates as 
#' in Delaigle, Hall and Meister (2008) and the estimator of the density of \eqn{X} is 
#' computed as in Delaigle, Hall and Meister (2008) except that we do the tail correction of 
#' the estimated characteristic function of the errors as in Delaigle and Hall (2016)
#' and Camirand, Carroll and Delaigle (2018). 
#-------------------------------------------------------------------------------
# NEED TO SAY WHICH VERSION WE USE EXACTLY. DO YOU REMEMBER WHAT WE TAKE FROM CAMIRAND ET AL?
#-------------------------------------------------------------------------------
#' 
#' \strong{Unknown heteroscedastic error distribution when replicates are available:} If both  
#' \code{W1} and \code{W2} are supplied and \code{het_replicates} is \code{TRUE}, then these 
#' distributions are estimated using replicates, as in Delaigle and Meister (2008) 
#' and the estimator of the density of \eqn{X} is computed as in Delaigle and Meister (2008)
#' except that we do the tail correction of the estimated pooled characteristic function of 
#' the errors as in Delaigle and Hall (2016) and Camirand, Carroll and Delaigle (2018).
#-------------------------------------------------------------------------------
# NEED TO SAY WHICH REFINEMENTS EXACTLY. DO YOU REMEMBER?
#-------------------------------------------------------------------------------
#' 
#' \strong{Unknown homoscedastic error distribution estimated without replicates:} 
#' If none of \code{errortype}, \code{phiU}, or \code{W2} are supplied then the density 
#' of \eqn{X} is assumed asymmetric and to satisfy the identifiability conditions of
#' Delaigle and Hall (2016). The estimator of the density of \eqn{X} is computed as in 
#' Delaigle and Hall (2016). 
#' 
#' The order in which we choose the methods is as follows:
#' \enumerate{
#' 	\item If provided, use \code{phiU} to define the errors, otherwise
#' 	\item If provided use \code{errortype} and \code{sd_u} to define the errors, otherwise
#' 	\item If provided, use the vector of replicates \code{W2} to estimate the error distribution, otherwise
#' 	\item We use the method for unknown homoscedastic error distribution estimated without replicates.
#' }
#' 
#' Note that in both 1 and 2, if a vector of replicates \code{W2} is provided we
#' augment the data in \code{W1} with that in \code{W2}.
#' 
#' 
#' @param W1 A vector of size n containing the univariate contaminated data.
#' @param W2 (optional) A vector of size n containing replicate measurements for the same 
#' n individuals (in the same order) as W1. If supplied, then the error distribution
#' will be estimated using the replicates only if \code{phiU}, and both of \code{errortype}
#' and \code{sd_U} are not provided.
#' @param xx A vector of x values on which to compute the estimator of the density of \eqn{X}.
#' @param errortype A single string giving the distribution of \eqn{U}, either "laplace" or "normal". 
#' If you define the error distribution this way then you must also provide 
#' \code{sd_U} but should not provide \code{phiU}. Argument is case-insensitive
#' and partially matched.
#' @param sd_U The standard deviation(s) of \eqn{U}: a single value for
#' homoscedastic errors and a vector of length \eqn{n} for heteroscedastic
#' errors. This does not need to be provided if you define your error distribution
#' using\code{phiU} and provide \code{bw}.
#' @param phiU Function(s) giving the characteristic function of the errors. A 
#' single function for homoscedastic errors and a vector of \eqn{n} functions 
#' for heteroscedastic errors. If you define the errors this way then you
#' should not provide \code{errortype}.
#' @param bw The bandwidth to use when computing the kernel estimator of the density
#' of \eqn{X}. If \code{NULL}, a bandwidth will be calculated using a plug-in estimator.
#' @param rescale If \code{TRUE}, the estimator of the density of \eqn{X} is rescaled so 
#' that it integrates to 1. Rescaling requires \code{xx} to be a fine grid of equispaced 
#' \eqn{x} values that cover the whole range of \eqn{x}-values where the 
#' estimated density is significantly non zero.
#' @param kernel_type The kernel K to use when computing the estimator of the 
#' density of \eqn{X}. The default kernel has characteristic function 
#' \eqn{(1-t^2)^3} for \eqn{t \in [-1,1]}. The normal kernel is the standard normal density.
#' The sinc kernel has characteristic function equal to 1 for \eqn{t \in [-1,1]}
#-------------------------------------------------------------------------------
# IS IT THE CASE THE ONLY ONE OF THOSE THREE KERNELS CAN BE USED, I.E. THE USER CANNOT PROVIDE THEIR OWN KERNEL?
# IS THERE AN ERROR MESSAGE IF THE USER PROVIDES A WRONG KERNEL?
#-------------------------------------------------------------------------------
#' @param m The number of point masses to use to estimate the distribution of 
#' \eqn{X} when the error distribution is not supplied and we use the method of
#' Delaigle and Hall (2016).
#-------------------------------------------------------------------------------
# NEED TO MENTION YOU USE A MODIFIED VERSION OF DH (2016)
# AND WHAT IS DIFFERENT FROM THERE.
#-------------------------------------------------------------------------------
#' @param show_diagnostics If \code{TRUE}, then diagnostic messages are printed 
#' displaying the results of the various optimizations performed when the error
#' distribution is not supplied and estimated by the method in Delaigle and
#' Hall (2016). Intended to be used for developement only.
#' @param het_replicates If \code{TRUE}, then the errors are not assumed to be 
#' homoscedastic and the code for heteroscedastic errors is used. Only applicable 
#' if \code{W2} is supplied.
#' 
#' @return An object of class "\code{deconvolve}".
#' 
#' The function \code{plot} produces a plot of the deconvolution KDE of the density of \eqn{X} on the grid \code{xx}.
#' 
#' An object of class "\code{deconvolve}" is a list containing at least the 
#' following elements:
#' \item{W1}{The original vector of contaminated data}
#' \item{x}{The values on which the deconvolution KDE is evaluated.}
#' \item{pdf}{A vector containing the deconvolution KDE of the density of \eqn{X}, 
#' evaluated at each point in \code{xx}}
#' It may also contain
#' \item{W2}{The original vector of replicates}
#' 
#' @section Warnings:
#' \itemize{
#'  \item If your errors have a Laplace distribution, \code{sd_U} is not the usual parameter
#'  	of a Laplace distribution but its standard deviation.
#'	\item The method for deconvolution when the error distribution is unknown and assumed
#'   symmetric, and estimated without replicates, as in Delaigle and Hall (2016), requires solving 
#' 	 multiple non-linear optimizations with both linear and non-linear 
#' 	 constraints. The current implementation can be slow and unreliable. An 
#' 	 alternative MATLAB implementation can be found at 
#' 	 <github.com/TimothyHyndman/deconvolve-supp> which may work better in some 
#' 	 circumstances.
#'	\item If you supply your own bandwidth, then you should ensure that the
#' 	kernel used here matches the one you used to calculate your bandwidth.
#'	\item The DKDE can also be computed using the Fast Fourier Transform, which 
#' 	is a bit more complex. See Delaigle and Gijbels (2007). However if the grid of 
#' 	t-values is fine enough, the estimator can simply be computed like here 
#' 	without having problems with oscillations.
#' }
#' 
#' @section References:
#' 
#' Camirand, F., Carroll, R.J. and Delaigle, A. (2018). Estimating the  
#' distribution of episodically consumed food measured with errors.  

#' Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating 
#' integrals and optimizing objective functions: a case study in density 
#' deconvolution. \emph{Statistics and Computing}, 17, 349-355.
#' 
#' Delaigle, A. and Hall, P. (2016). Methodology for non-parametric 
#' deconvolution when the error distribution is unknown. \emph{Journal of the 
#' Royal Statistical Society: Series B (Statistical Methodology)}, 78, 1, 
#' 231-252.
#' 
#' Delaigle, A., Hall, P. and Meister, A. (2008). On Deconvolution with  
#' repeated measurements. \emph{Annals of Statistics}, 36, 665-685 
#' 
#' Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic 
#' error. \emph{Bernoulli}, 14, 2, 562-579.
#' 
#' Stefanski, L.A. and Carroll, R.J. (1990). Deconvolving kernel density
#' estimators. \emph{Statistics}, 21, 2, 169-184.
#' \emph{Manuscript.} 
#' 
#' @author Aurore Delaigle, Timothy Hyndman, Tianying Wang
#' 
#' @example man/examples/deconvolve_eg.R
#' 
#' @export

deconvolve <- function(W1, W2 = NULL, xx = seq(min(W1), max(W1), length.out = 100), 
					   errortype = NULL, sd_U = NULL, phiU = NULL, bw = NULL, 
					   rescale = FALSE,
					   kernel_type = c("default", "normal", "sinc"), 
					   het_replicates = FALSE,
					   m = 20,
					   show_diagnostics = FALSE){

	# Partial matching ---------------------------------------------------------
	dist_types <- c("normal", "laplace")
	if (!is.null(errortype)) {
		errortype <- dist_types[pmatch(tolower(errortype), dist_types)]
		if (is.na(errortype)) {
			stop("Please provide a valid errortype.")
		}
	}

	kernel_type <- match.arg(kernel_type)

	# Determine error type provided --------------------------------------------
	if (!is.null(phiU)) {
		if (length(phiU) > 1) {
			errors <- "het"
		} else {
			errors <- "hom"
		}
	} else if (!is.null(errortype) & !is.null(sd_U)) {
		if (length(sd_U) > 1) {
			errors <- "het"
		} else {
			errors <- "hom"
		}
	} else if (!is.null(W2)) {
		if (het_replicates) {
			errors <- "het_rep"
		} else {
			errors <- "rep"
		}
	} else {
		errors <- "sym"
		warning("The method for deconvolution when the error is unknown is slow and unreliable in R. Consider instead using the MATLAB code found at <github.com/TimothyHyndman/deconvolve-supp>.")
	}

	# Augment W1 with W2 if provided along with phiU or sd_U and errortype
	if ((errors == "het" | errors == "hom") & !is.null(W2)) {
		W1 <- c(W1, W2)
		W2 <- NULL
		if (!is.null(phiU)) {
			warning("Both phiU and W2 have been provided. Continuing using errors defined by phiU and augmenting W1 with the data in W2.")
		} else {
			warning("Errortype and sd_U as well as W2 have been provided. Continuing using errors defined by errortype and sd_U and augmenting W1 with the data in W2.")
		}
		
	}

	# Check inputs -------------------------------------------------------------
	if (errors == "het") {
		if (is.null(phiU)) {
			if ((length(sd_U) == length(W1)) == FALSE) {
				stop("sd_U must be either length 1 for homoscedastic errors or have the same length as W1 for heteroscedastic errors.")
			}
		} else {
			if ((length(phiU) == length(W1)) == FALSE) {
				stop("phiU must be either length 1 for homoscedastic errors or have the same length as W1 for heteroscedastic errors.")
			}
		}
	}

	if (!is.null(errortype) & is.null(sd_U)) {
		stop("You must provide sd_U along with errortype.")
	}

	if (!is.null(phiU) & is.null(bw) & is.null(sd_U)){
		stop("You must provide sd_U along with phiU if you do not provide bw.")
	}

	if (errors == "rep"){
		if (!(length(W1) == length(W2))) {
			stop("W1 and W2 must be the same length.")
		}
	}

	if (kernel_type == "normal") {
		warning("You should only use the 'normal' kernel when the errors are ordinary 
			smooth such as Laplace or convolutions of Laplace.")
	}

	if (kernel_type == "sinc") {
		warning("You should ensure that you are not using a plug-in bandwidth 
			method for the bandwidth.")
	}

	# Calculate Bandwidth if not supplied --------------------------------------
	if (is.null(bw) & !(errors == "sym")) {
	#DON'T YOU NEED TO PROVIDE THE ALGORITHM TO USE IN BANDWIDTH?
			bw <- bandwidth(W1, W2, errortype, sd_U, phiU, 
							kernel_type = kernel_type)
	}

	# --------------------------------------------------------------------------
	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	muK2 <- kernel_list$muk2
	RK <- kernel_list$rk
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]
	
	# Convert errortype to phiU ------------------------------------------------
	if ((errors == "hom") | (errors == "het")){
		if(is.null(phiU)) {
			phiU <- create_phiU(errors, errortype, sd_U)
		}
	}

	# Perform appropriate deconvolution ----------------------------------------
	if (errors == "hom"){
		pdf <- DeconErrKnownPdf(xx, W1, bw, phiU, kernel_type, rescale)
		output <- list("x" = xx, "pdf" = pdf, "W1" = W1)
	}

	if (errors == "het"){
		pdf <- DeconErrKnownHetPdf(xx, W1, bw, phiU, rescale, phiK, muK2, RK, tt)
		output <- list("x" = xx, "pdf" = pdf, "W1" = W1)
	}

	if (errors == "rep") {
		phi_U <- create_replicates_phi_U(W1, W2, tt/bw)
		pdf <- DeconErrKnownPdf(xx, c(W1, W2), bw, phi_U, kernel_type, rescale)
		output <- list("x" = xx, "pdf" = pdf, "W1" = W1, "W2" = W2)
	}

	if (errors == "het_rep") {
		pdf <- decon_err_het_replicates(xx, W1, W2, kernel_type, bw, rescale)
		output <- list("x" = xx, "pdf" = pdf, "W1" = W1, "W2" = W2)
	}

	if (errors == "sym") {
		out <- DeconErrSymPmf(W1, m, kernel_type, show_diagnostics = show_diagnostics)
		phi.W <- out$phi_W
		pdf <- DeconErrSymPmfToPdf(out, W1, phi.W, xx, kernel_type, rescale, 
								   bw)
		output <- list("x" = xx, "pdf" = pdf, "support" = out$support, 
					   "probweights" = out$probweights, "W1" = W1)
	}

	# Output object of class "deconvolve" --------------------------------------
	class(output) <- c("deconvolve", "list")
	output
}
