create_replicates_phi_U  <- function(W1, W2, t_search) {
	# Estimate phiU from replicates as in Delaigle, Hall and Meister (2008) but with tail adjustments
	# as in Camirand, Carroll and Delaigle (2018).
	
	#Compute U1-U2
	diff <- W1 - W2
	#Estimate standard deviation of U
	sd_U <- sqrt(stats::var(diff)/2)
	#Number of replicates
	n_diff <- length(diff)
	#Estimate characteristic function of U at t_search from replicates
	tout <- outer(t_search, diff)
	phi_U_rep <- rowSums(cos(tout))/n_diff
	phi_U_rep[phi_U_rep < 0] <- 0
	phi_U_rep <- sqrt(phi_U_rep)

	# Find range of t for which phi_U_rep is reliable. Plays the role of the interval A 
	# at page 678 of Delaigle, Hall and Meister (2008). 
	# Cutoff is computed as in Camirand, Carroll and Delaigle (2018), which is a refined
	# version of the cutoff proposed by Delaigle and Hall (2016) 
	
	t_cutoff <- find_t_cutoff(phi_U_rep, t_search)
	
	# Outside the reliable t-range, pretend the error is Laplace. Elsewhere 
	#produce a smooth estimator of phiU via a spline approximation.
	phi_U <- function(t) {
		phiU_spline(t, sd_U, t_cutoff, t_search, phi_U_rep)
	}

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
	#Outside the reliable t-range, pretend the error is Laplace
	y[ind2] <- phiU_lap(t[ind2])

	y
}