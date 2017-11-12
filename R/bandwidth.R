#' Bandwidth Selectors for Deconvolution Kernel Density Estimation
#'
#' Computes a bandwidth for use in deconvolution kernel density estimation of
#' \eqn{X} from data \eqn{W = X + U}. If 'SIMEX' algorithm used, computes a
#' bandwidth for use in deconvolution regression of data \eqn{(W, Y)}  where
#' \eqn{Y = g(X) + V} and \eqn{W = X + U}.
#'
#' The function \code{bandwidth} chooses from one of five different methods
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
#' \strong{PI for Unknown Error:} If \code{algorithm = "PI"} and the errors are
#' not supplied, then the error is estimated using the method described in
#' Delaigle and Hall 2016 and then the bandwidth is calculated using the method
#' described in Delaigle and Gijbels 2002, and Delaigle and Gijbels 2004.
#'
#' \strong{CV:} If \code{algorithm = "CV"} then the method used is as described
#' in Stefanski and Carroll 1990, and Delaigle and Gijbels 2004.
#'
#' \strong{SIMEX:} If \code{algorithm = "SIMEX"} then the method used is as
#' described in Delaigle and Hall 2008.
#'
#' @inheritParams deconvolve
#' @param algorithm Either \code{"PI"} for plug-in estimator or \code{"CV"} for
#' cross-validation estimator. If \code{"CV"} then the errors must be
#' homoscedastic.
#' @param Y A vector of the univariate dependant data. Only required for 'SIMEX'
#' algorithm.
#' @param n_cores Number of cores to use when using SIMEX algorithm. If
#' \code{NULL}, the number of cores to use will be automatically detected.
#' @param sigU The standard deviations of \eqn{U}. A single value for
#' homoscedastic errors and a vector having the same length as \code{W} for
#' heteroscedastic errors.
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
#' Delaigle, A. and Hall, P. (2016). Methodology for non-parametric
#' deconvolution when the error distribution is unknown. \emph{Journal of the
#' Royal Statistical Society: Series B (Statistical Methodology)}, 78, 1,
#' 231-252.
#'
#' @author Aurore Delaigle, Timothy Hyndman, Tianying Wang
#'
#' @example man/examples/bandwidth_eg.R
#'
#' @export

bandwidth <- function(W, errortype = NULL, sigU = NULL, phiU = NULL, Y = NULL,
					  algorithm = "PI", n_cores = NULL, kernel_type = "default"){

	# Determine error type provided --------------------------------------------
	if (is.null(errortype) & is.null(phiU)) {
		errors <- "sym"
	} else if (length(sigU) > 1  | length(phiU) > 1){
		errors <- "het"
	} else {
		errors <- "hom"
	}

	# Check inputs -------------------------------------------------------------
	if (errors == "het") {
		if (is.null(phiU)) {
			if ((length(sigU) == length(W)) == FALSE) {
				stop("sigU must be either length 1 for homoscedastic errors or have the same length as W for heteroscedastic errors.")
			}
		} else {
			if ((length(phiU) == length(W)) == FALSE) {
				stop("phiU must be either length 1 for homoscedastic errors or have the same length as W for heteroscedastic errors.")
			}
		}
	}

	if ((errors == "het" | errors == "hom")  & is.null(sigU)) {
		stop("You must provide sigU along with the errors.")
	}

	if (is.null(errortype) == FALSE) {
		if ((errortype == "norm" | errortype == "Lap") == FALSE) {
			stop("errortype must be one of: 'norm', or 'Lap'.")
		}
	}

	if ((algorithm == "CV" | algorithm == "PI" | algorithm == "SIMEX") == FALSE) {
		stop("algorithm must be one of: 'PI', 'CV', or 'SIMEX'.")
	}

	if (algorithm == "CV") {
		if (!(errors == "hom")) {
			stop("Algorithm type 'CV' can only be used with homoscedastic
				 errors.")
		}
	}

	# if ((errors == "sym") & (!(algorithm == "PI"))) {
	# 	stop("You must use the PI algorithm if the errors are not supplied.")
	# }

	if (missing(sigU) & (errors == "hom" & algorithm == "PI")) {
		stop("You must supply sigU to use the PI algorithm when the errors are homoscedastic.")
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

	# if (errors == "sym" & is.null(sigU)) {
	# 	stop("You must provide an estimate for sigU when the errors are estimated.")
	# }

	# --------------------------------------------------------------------------
	n <- length(W)

	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	muK2 <- kernel_list$muk2
	RK <- kernel_list$rk
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]

	# Convert errortype to phiU ------------------------------------------------

	if (!(errors == 'est')) {
		if(is.null(phiU)) {
			phiU <- create_phiU(errors, errortype, sigU)
		}
	}

	# Calculate varX if in het case --------------------------------------------
	if (errors == "het"){
		n <- length(W)
		varX <- max(mean(W^2) - (mean(W))^2 - sum(sigU^2) / n, 1/n)	#max 1/n avoid negative varianace
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
		output <- PI_deconvUknownth4(n, W, sigU, phiU, phiK, muK2, RK, deltat, tt)
	}

	if (algorithm == "PI" & errors == "sym") {
		d <- DeconErrSymPmf(W, 10)
		theta <- d$support
		p <- d$probweights
		t <- tt
		tt <- d$phi.W$t.values

		# Estimate Var(U) ----------------------------------------------------------
		tt.BB.length <- 200		# Use a finer grid than tt
		tt.BB <- seq(tt[1], tt[length(tt)], length.out = tt.BB.length)
		hat.var.U <- estimate_var_u(W, tt.BB, theta, p)

		# Estimate PhiX and PhiU ---------------------------------------------------
		phi.X <- ComputePhiPmf(theta, p, tt)
		phi.U <- d$phi.W$norm / Mod(phi.X)

		# Actually find bandwidth
		output <- PI_DeconvUEstTh4(W, phi.U, hat.var.U, tt, phiK, muK2, t)
	}

	if (algorithm == "SIMEX") {
		output <- hSIMEXUknown(W, Y, errortype, sigU, phiK, muK2, RK, deltat,
							   tt, n_cores)
	}

	output
}
