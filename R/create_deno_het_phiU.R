create_deno_het_phi_U  <- function(W1, W2, t_search) {
	# Estimate phiU from replicates as in Delaigle, Hall and Meister (2008) but with tail adjustments
	# as in Camirand, Carroll and Delaigle (2018).
	
	#Compute means and differences
	diff2 <- (W1 - W2)/2
	#Estimate 1/2 var of U
	varU2 <- stats::var(diff2)
	#Number of replicates
	n_diff <- length(diff2)
	#Estimate characteristic function of U at t_search from replicates
	tout <- outer(t_search, diff2)
	het_deno <- abs(rowSums(cos(tout)))
	het_deno[het_deno < 0] <- 0

	# Find range of t for which het_deno is reliable. Plays the role of the interval A 
	# at page 678 of Delaigle, Hall and Meister (2008). 
	# Cutoff is computed as in Camirand, Carroll and Delaigle (2018), which is a refined
	# version of the cutoff proposed by Delaigle and Hall (2016) 
	
	#divide deno by n when computing the cutoff as otherwise it is not an empirical characteristic function
	t_cutoff <- find_t_cutoff(het_deno/n_diff, t_search,n_diff)
	
	# Outside the reliable t-range, pretend the error is Laplace. Elsewhere 
	#produce a smooth estimator of phiU via a spline approximation.
	phi_U <- function(t) {
		phiU_splineB(t, varU2, t_cutoff, t_search, het_deno,n_diff)
	}

}

phiU_splineB <- function(t, varU2, t_cutoff, t_phiU, het_deno,n_diff) {
	ind1 <- (t >= -t_cutoff) & (t <= t_cutoff)
	ind2 <- !ind1
	
	
	phiU_lap2 <- function(t){
		n_diff / (1 + varU2 * (t/2)^2)^2
	}

	y <- 0*t

	if (any(ind1)) {
		y[ind1] <- stats::spline(t_phiU, het_deno, xout = t[ind1])$y 	
	}
	#Outside the reliable t-range, pretend the error is Laplace
	y[ind2] <- phiU_lap2(t[ind2])

	y
}
