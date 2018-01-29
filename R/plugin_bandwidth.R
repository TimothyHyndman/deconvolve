plugin_bandwidth <- function(W, phi_U, sd_X, kernel_type) {
	 

	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	mu_K2 <- kernel_list$muk2
	RK <- kernel_list$rk
	tt <- kernel_list$tt


	# Grid of h on which to search for a solution
	maxh <- (max(W) - min(W)) / 10
	hnaive <- ((8 * sqrt(pi) * RK/3/mu_K2^2)^0.2) * sqrt(stats::var(W)) * 
	n^(-1/5)
	hgrid <- seq(hnaive / 3, maxh, length.out = 100)

	# Quantities that will be needed multiple times below
	toverh <- tt %*% t(1 / hgrid)
	phiKsq <- phiK(tt)^2
	phiUsq <- phiU(toverh)^2

	# theta 4
	th4 <- (sd_X^(-9) * factorial(8) ) / (2^9 * factorial(4) * sqrt(pi))

	# h3
	rr <- 3
	term1 <- -hgrid^2 * mu_K2 * th4
	term2 <- kronecker(matrix(1, 1, lh), tt^(2*rr) * phiKsq) / phiUsq
	term2 <- apply(term2, 2, sum) * deltat
	term2 <- term2 / (2 * pi * n * hgrid^(2 * rr + 1))

	ABias2 <- (term1 + term2)^2

	indh3 <- which.min(ABias2)
	if (indh3==1) {
		cat("\n minimum of Abias2 for rr=3 is the first element of the grid of 
		bandwidths. Consider enlarging the grid \n")
	}
	if(indh3==length(hgrid)){
		cat("\n minimum of Abias2 for rr=3 is the last element of the grid of 
		bandwidths. Consider enlarging the grid \n")
	}
	h3 <- hgrid[indh3]
}