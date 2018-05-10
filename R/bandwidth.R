#' Bandwidth Selectors for Deconvolution Kernel Density Estimation
#'
#' Computes a bandwidth for use in deconvolution kernel density estimation of
#' \eqn{X} from data \eqn{W = X + U} when the distribution of \eqn{U} is known,
#' unknown, or estimated from replicates, \eqn{W_2 = X + U_2}. If 'SIMEX' 
#' algorithm used, computes a bandwidth for use in deconvolution regression of 
#' data \eqn{(W, Y)}  where \eqn{Y = g(X) + V} and \eqn{W = X + U}.
#'
#' The function \code{bandwidth} chooses from one of six different methods
#' depending on how the error distribution is defined and which algorithm is
#' selected.
#'
#' \strong{PI for Homoscedastic Error:} If \code{algorithm = "PI"} and the errors
#' are defined by either a single function \code{phiU}, or a single value
#' \code{sd_U} along with its \code{errortype}, then the method used is as
#' described in Delaigle and Gijbels 2002, and Delaigle and Gijbels 2004.
#'
#' \strong{PI for Heteroscedastic Error:} If \code{algorithm = "PI"} and the
#' errors are defined by a either a vector of functions \code{phiU}, or a vector
#' \code{sd_U} along with its \code{errortype} then the method used is as
#' described in Delaigle and Meister 2008.
#'
#' \strong{PI for Replicates:} If \code{algorithm = "PI"} and a replicate 
#' vector \code{W2} is supplied, then the error is estimated using replicates.
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
#' \strong{SIMEX for Replicates:} If \code{algorithm = "SIMEX"} and a 
#' replicate vector \code{W2} is supplied, then \eqn{phi_U} is calculated using 
#' replicates and SIMEX is performed as according to 
#' \code{use_alt_SIMEX_rep_opt}.
#'
#' @inheritParams deconvolve
#' @param algorithm One of \code{"PI"} for plug-in estimator, \code{"CV"} for
#' cross-validation estimator or \code{"SIMEX"}. If \code{"CV"} then the errors 
#' must be homoscedastic.
#' @param Y A vector of the univariate dependent data. Only required for 'SIMEX'
#' algorithm.
#' @param n_cores Number of cores to use when using SIMEX algorithm. If
#' \code{NULL}, the number of cores to use will be automatically detected.
#' @param sd_U The standard deviations of \eqn{U}. A single value for
#' homoscedastic errors and a vector having the same length as \code{W} for
#' heteroscedastic errors.
#' @param seed Set seed for SIMEX. Allows for reproducible results using SIMEX.
#' @param use_alt_SIMEX_rep_opt Only used with SIMEX using replicates. If 
#' \code{TRUE}, performs SIMEX on \eqn{W = (W1 + W2)/2} and samples \eqn{U*} 
#' from (W1 - W2). The default performs SIMEX on \eqn{W = (W1, W2)} and 
#' and samples \eqn{U*} from \eqn{(W1 - W2)/\sqrt 2}.
#'
#' @return The bandwidth estimator. If using 'SIMEX' algorithm then returns a
#' list containing the bandwidth 'h' and ridge parameter 'rho'.
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
#' Delaigle, A., Hall, P., and Meister, A. (2008). On Deconvolution with  
#' repeated measurements. \emph{Annals of Statistics}, 36, 665-685 
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

bandwidth <- function(W, 
					  W2 = NULL, 
					  errortype = NULL, 
					  sd_U = NULL, 
					  phiU = NULL, 
					  Y = NULL,
					  algorithm = c("PI", "CV", "SIMEX"), 
					  n_cores = NULL,
					  kernel_type = c("default", "normal", "sinc"), 
					  seed = NULL,
					  use_alt_SIMEX_rep_opt = FALSE){

	# Determine error type provided --------------------------------------------
	if (!is.null(W2)) {
		errors <- "rep"
	} else if (is.null(errortype) & is.null(phiU)) {
		errors <- "sym"
	} else if (length(sd_U) > 1  | length(phiU) > 1){
		errors <- "het"
	} else {
		errors <- "hom"
	}

	# Partial matching ---------------------------------------------------------
	dist_types <- c("normal", "laplace")
	if (!is.null(errortype)) {
		errortype <- dist_types[pmatch(tolower(errortype), dist_types)]
		if (is.na(errortype)) {
			stop("Please provide a valid errortype.")
		}
	}

	algorithm <- match.arg(algorithm)
	kernel_type <- match.arg(kernel_type)

	# Check inputs -------------------------------------------------------------
	if (errors == "het") {
		if (is.null(phiU)) {
			if ((length(sd_U) == length(W)) == FALSE) {
				stop("sd_U must be either length 1 for homoscedastic errors or have the same length as W for heteroscedastic errors.")
			}
		} else {
			if ((length(phiU) == length(W)) == FALSE) {
				stop("phiU must be either length 1 for homoscedastic errors or have the same length as W for heteroscedastic errors.")
			}
		}
	}

	if ((errors == "het" | errors == "hom")  & is.null(sd_U)) {
		stop("You must provide sd_U along with the errors.")
	}

	if (errors == "rep"){
		if (!(length(W) == length(W2))) {
			stop("W and W2 must be the same length.")
		}
	}

	# if ((algorithm == "CV" | algorithm == "PI" | algorithm == "SIMEX") == FALSE) {
	# 	stop("algorithm must be one of: 'PI', 'CV', or 'SIMEX'.")
	# }

	if (algorithm == "CV") {
		if (!(errors == "hom")) {
			stop("Algorithm type 'CV' can only be used with homoscedastic
				 errors.")
		}
	}

	# if ((errors == "sym") & (!(algorithm == "PI"))) {
	# 	stop("You must use the PI algorithm if the errors are not supplied.")
	# }

	if (is.null(sd_U) & (errors == "hom" & algorithm == "PI")) {
		stop("You must supply sd_U to use the PI algorithm when the errors are homoscedastic.")
	}

	if (algorithm == "SIMEX") {
		if (is.null(Y)) {
			stop("You must supply Y to use SIMEX.")
		}
		if (errors == "het") {
			stop("Algorithm type 'SIMEX' can only be used with homoscedastic
				 errors.")
		}
		if ((is.null(sd_U) | is.null(errortype)) & is.null(W2)) {
			stop("Algorithm 'SIMEX' requires that the errors are provided using either W2 or errortype and sd_U.")
		}
	}

	if (kernel_type == "normal") {
		warning("You should only use the 'normal' kernel when the errors are 
			Laplace or convolutions of Laplace.")
	}

	if (kernel_type == "sinc" & algorithm == "PI") {
		stop("You cannot use the 'sinc' kernel with the plug-in algorithm.")
	}
	# if (errors == "sym" & is.null(sd_U)) {
	# 	stop("You must provide an estimate for sd_U when the errors are estimated.")
	# }

	# --------------------------------------------------------------------------
	n <- length(W)

	if (algorithm == "SIMEX") {
		kernel_list <- kernel(kernel_type, coarse = TRUE)
	} else {
		kernel_list <- kernel(kernel_type)
	}

	phiK <- kernel_list$phik
	muK2 <- kernel_list$muk2
	RK <- kernel_list$rk
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]

	# Convert errortype to phiU ------------------------------------------------
	if ((errors == 'hom') | (errors == 'het') ) {
		if(is.null(phiU)) {
			phiU <- create_phiU(errors, errortype, sd_U)
		}
	}

	# Perform appropriate bandwidth calculation --------------------------------

	if (algorithm == "CV"){
		output <- CVdeconv(n, W, phiU, phiK, muK2, RK, deltat, tt)
	}

	if (algorithm == "SIMEX" & errors == "hom") {
		generate_U_star <- create_generate_U_star(W, W2, errortype, sd_U)
		output <- hSIMEXUknown(W, Y, generate_U_star, sd_U, phiU, kernel_type, 
							   n_cores, seed)
	}

	if (algorithm == "SIMEX" & errors == "rep") {
		diff <- W - W2
		diff <- diff[(W != 0) & (W2 != 0)]
		sd_U <- sqrt(stats::var(diff)/2)
		n <- length(W)
		hnaive <- ((8 * sqrt(pi) * RK/3/muK2^2)^0.2) * 
			sqrt(stats::var(W)) * n^(-1/5)
		h_min <- hnaive / 3
		t_search <- tt/h_min
		phi_U <- create_replicates_phi_U(W, W2, t_search)
		generate_U_star <- create_generate_U_star(W, W2, errortype, sd_U, use_alt_SIMEX_rep_opt)

		if (use_alt_SIMEX_rep_opt) {
			W_bar <- (W + W2)/2
			output <- hSIMEXUknown(W_bar, Y, generate_U_star, sd_U, phi_U, kernel_type, 
								   n_cores, seed)
		} else {
			W_full <- c(W, W2)
			Y_full <- c(Y, Y)
			output <- hSIMEXUknown(W_full, Y_full, generate_U_star, sd_U, phi_U, kernel_type, 
								   n_cores, seed)
		}
	}

	if (algorithm == "PI" & errors == "het") {
		n <- length(W)
		varX <- max(mean(W^2) - (mean(W))^2 - sum(sd_U^2) / n, 1/n)
		output <- PI_deconvUknownth4het(n, W, varX, phiU, phiK, muK2, RK, 
										deltat, tt)
	}

	if (algorithm == "PI" & errors == "hom") {
		sd_X <- max( !is.na(sqrt(stats::var(W) - sd_U^2)), 1/n)
		output <- plugin_bandwidth(W, phiU, sd_X, kernel_type)
	}

	if (algorithm == "PI" & errors == "rep") {
		diff <- W - W2
		diff <- diff[(W != 0) & (W2 != 0)]
		sd_U <- sqrt(stats::var(diff)/2)
		n <- length(c(W, W2))
		sd_X <- max( !is.na(sqrt(stats::var(c(W, W2)) - sd_U^2)), 1/n)
		hnaive <- ((8 * sqrt(pi) * RK/3/muK2^2)^0.2) * 
			sqrt(stats::var(c(W, W2))) * n^(-1/5)
		h_min <- hnaive / 3
		t_search <- tt/h_min
		phi_U <- create_replicates_phi_U(W, W2, t_search)

		output <- plugin_bandwidth(c(W, W2), phi_U, sd_X, kernel_type)
	}

	if (algorithm == "PI" & errors == "sym") {
		warning("The plug-in bandwidth method when the error is unknown and assumed symmetric is slow and unreliable in R. Consider instead using the MATLAB code found at <URL>.")
		d <- DeconErrSymPmf(W, 10, kernel_type)
		theta <- d$support
		p <- d$probweights
		t <- tt
		tt <- d$phi_W$t.values

		# Estimate Var(U) ------------------------------------------------------
		tt.BB.length <- 200		# Use a finer grid than tt
		tt.BB <- seq(tt[1], tt[length(tt)], length.out = tt.BB.length)
		sd_U <- sqrt(estimate_var_u(W, tt.BB, theta, p))

		# Estimate PhiX and PhiU -----------------------------------------------
		phi.X <- ComputePhiPmf(theta, p, tt)
		phi.U <- d$phi_W$norm / Mod(phi.X)

		t_cutoff <- tt[length(tt)]
		phi_U_splined <- function(t) {
			phiU_spline(t, sd_U, t_cutoff, tt, phi.U)
		}
		
		# Actually find bandwidth
		sd_X <- max(!is.na(sqrt(stats::var(W) - sd_U^2)), 1 / n)
		output <- plugin_bandwidth(W, phi_U_splined, sd_X, kernel_type)
	}

	output
}

create_generate_U_star <- function(W, W2, errortype, sd_U, use_alt_SIMEX_rep_opt) {
	# Generates vector of length n with same distribution as U
	n <- length(W)

	if (!is.null(W2)) {
		if (use_alt_SIMEX_rep_opt) {
			generate_U_star <- function() {
				U_star <- sample((W - W2), replace = TRUE)
				U_star	
			}
		} else {
			generate_U_star <- function() {
				U_star <- sample((W - W2)/sqrt(2), size = 2*n, replace = TRUE)
				U_star	
			}
		}
	}

	if (!is.null(errortype)) {
		if (errortype == "laplace") {
			generate_U_star <- function() {
				U_star <- rlap(sd_U/sqrt(2), 1, n)	
			}
		}

		if (errortype == "normal") {
			generate_U_star <- function() {
				U_star <- stats::rnorm(n, 0, sd_U)	
			}
		}
	}

	generate_U_star
}