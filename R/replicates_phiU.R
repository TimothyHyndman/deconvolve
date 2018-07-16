create_replicates_phi_U  <- function(W1, W2, t_search) {
	# Estimate phiU from replicates
	diff <- W1 - W2
	sd_U <- sqrt(stats::var(diff)/2)

	n_diff <- length(diff)
	tout <- outer(t_search, diff)
	phi_U_rep <- rowSums(cos(tout))/n_diff
	phi_U_rep[phi_U_rep < 0] <- 0
	phi_U_rep <- sqrt(phi_U_rep)

	# Find range of t for which phi_U_rep is reliable
	t_cutoff <- find_t_cutoff(phi_U_rep, t_search)

	phi_U <- function(t) {
		phiU_spline(t, sd_U, t_cutoff, t_search, phi_U_rep)
	}

	phi_U
}

phiU_spline <- function(t, sd_U, t_cutoff, t_phiU, phiU) {
	ind1 <- (t >= -t_cutoff) & (t <= t_cutoff)
	ind2 <- !ind1

	phiU_lap <- function(t){
		1 / (1 + sd_U^2 / 2 * t^2)
	}

	y <- 0*t

	if (any(ind1)) {
		y[ind1] <- stats::spline(t_phiU, phiU, xout = t[ind1])$y 	
	}
	y[ind2] <- phiU_lap(t[ind2])

	y
}