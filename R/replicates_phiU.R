replicates_phiU <- function(t, W1, W2, t_search) {
	
	# Estimate phiU from replicates
	diff <- W1 - W2
	diff <- diff[(W1 != 0) & (W2 != 0)]
	sd_U <- sqrt(stats::var(diff)/2)

	n_diff <- length(diff)
	tout <- outer(t_search, diff)
	phiU <- rowSums(cos(tout))/n_diff
	phiU[phiU < 0] <- 0
	phiU <- sqrt(phiU)

	# Find range of t for which phiU is reliable
	tmp <- t_search[phiU < n_diff^(-1/4)]
	if (length(tmp) > 1) {
		t_cutoff <- min(tmp[tmp > 0])	
	} else {
		t_cutoff <- Inf
	}
	

	phiU_spline(t, sd_U, t_cutoff, t_search, phiU)

}

phiU_spline <- function(t, sd_U, t_cutoff, t_phiU, phiU) {
	ind1 <- (t >= -t_cutoff) & (t <= t_cutoff)
	ind2 <- !ind1

	phiU_lap <- function(t){
		1 / (1 + sd_U^2 / 2 * t^2)
	}

	y <- 0*t
	y[ind1] <- stats::spline(t_phiU, phiU, xout = t[ind1])$y 	
	y[ind2] <- phiU_lap(t[ind2])

	y
}