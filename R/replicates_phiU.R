replicates_phiU <- function(t, W1, W2, kernel_type) {
	
	kernel_list <- kernel(kernel_type)
	mu_K2 <- kernel_list$muk2
	RK <- kernel_list$rk
	tt <- kernel_list$tt

	W <- c(W1, W2)
	n <- length(W)
	hnaive <- ((8 * sqrt(pi) * RK/3/mu_K2^2)^0.2) * sqrt(stats::var(W)) * 
		n^(-1/5)

	h_min <- hnaive / 3

	# Estimate phiU from replicates
	diff <- W1 - W2
	diff <- diff[(W1 != 0) & (W2 != 0)]
	sd_U <- sqrt(stats::var(diff)/2)

	n_diff <- length(diff)
	tU <- tt/h_min
	tout <- outer(tU, diff)
	phiU <- rowSums(cos(tout))/n_diff
	phiU[phiU < 0] <- 0
	phiU <- sqrt(phiU)

	# Find range of t for which phiU is reliable
	tmp <- tU[phiU < n_diff^(-1/4)]
	t_cutoff <- min(tmp[tmp > 0])

	phiU_spline(t, sd_U, t_cutoff, tU, phiU)

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