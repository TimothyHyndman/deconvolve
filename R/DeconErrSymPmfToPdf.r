DeconErrSymPmfToPdf <- function(X_pmf, W, phi_W, xx, kernel_type, rescale, h){

	theta <- X_pmf$support
	p <- X_pmf$probweights
	tt <- phi_W$t.values
	n <- length(W)

	# Estimate sd_U ------------------------------------------------------------
	tt_BB_length <- 200		# Use a finer grid than tt
	tt_BB <- seq(tt[1], tt[length(tt)], length.out = tt_BB_length)
	sd_U <- sqrt(estimate_var_u(W, tt_BB, theta, p))

	# Estimate PhiX and PhiU ---------------------------------------------------
	phi_X <- ComputePhiPmf(theta, p, tt)
	phi_U <- phi_W$norm / Mod(phi_X)

	t_cutoff <- tt[length(tt)]	# We have already found t_cutoff earlier when
								# calculating phiW
	phi_U_splined <- function(t) {
		phiU_spline(t, sd_U, t_cutoff, tt, phi_U)
	}
	
	# Find Plug-In Bandwidth ---------------------------------------------------
	if (is.null(h)) {
		sd_X <- max(sqrt(stats::var(W) - sd_U^2), 1 / n)
		h <- plugin_bandwidth(W, phi_U_splined, sd_X, kernel_type)
	}
	
	# --------------------------------------------------------------------------
	fX <- DeconErrKnownPdf(xx, W, h, phi_U_splined, kernel_type, rescale)
}

estimate_var_u <- function(W, tt_BB, theta, p){
	#---------------------------------------------------------------------------
	# Estimate var(U): approximate phi_U by poly of degree 2, and estimate varU 
	# by -2 * second order coefficient
	#---------------------------------------------------------------------------
	phi_W_BB <- ComputePhiEmp(W, tt_BB)
	phi_X_BB <- ComputePhiPmf(theta, p, tt_BB)
	phi_U_BB <- phi_W_BB$norm / Mod(phi_X_BB)

	t_vec <- 		  tt_BB[ phi_U_BB >= 0.95 ]
	phi_U_t_vec <- phi_U_BB[ phi_U_BB >= 0.95 ]

	pp <- stats::lm(phi_U_t_vec ~ stats::poly(t_vec, 2, raw = TRUE))
	
	-2 * pp$coefficients[[3]]
}