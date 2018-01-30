plugin_bandwidth <- function(W, phi_U, sd_X, kernel_type) {
	 

	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	mu_K2 <- kernel_list$muk2
	RK <- kernel_list$rk
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]


	# Grid of h on which to search for a solution
	maxh <- (max(W) - min(W)) / 10
	hnaive <- ((8 * sqrt(pi) * RK/3/mu_K2^2)^0.2) * sqrt(stats::var(W)) * 
		n^(-1/5)
	lh <- 100
	hgrid <- seq(hnaive / 3, maxh, length.out = lh)

	# Quantities that will be needed multiple times below
	toverh <- tt %*% t(1 / hgrid)
	phiKsq <- phiK(tt)^2
	phiUsq <- phiU(toverh)^2

	calculate_indh <- function(rr, th) {
		term1 <- -hgrid^2 * mu_K2 * th
		term2 <- kronecker(matrix(1, 1, lh), tt^(2*rr) * phiKsq) / phiUsq
		term2 <- apply(term2, 2, sum) * deltat
		term2 <- term2 / (2 * pi * n * hgrid^(2 * rr + 1))

		ABias2 <- (term1 + term2)^2

		indh <- which.min(ABias2)
		if (indh == 1) {
			warning("Minimum of Abias2 for rr=3 is the first element of the grid of 
			bandwidths. Consider enlarging the grid.")
		}
		if(indh == lh){
			warning("Minimum of Abias2 for rr=3 is the last element of the grid of 
			bandwidths. Consider enlarging the grid.")
		}

		indh
	}

	# theta 4
	th4 <- (sd_X^(-9) * factorial(8) ) / (2^9 * factorial(4) * sqrt(pi))

	# h3
	indh3 <- calculate_indh(3, th4)
	h3 <- hgrid[indh3]

	# theta 3
	phi_W <- ComputePhiEmp(W, tt/h3)
	th3 <- sum(tt^(2 * rr) * phi_W$norm^2 * phiKsq / phiUsq[, indh3])
	th3 <- th3 * deltat / (2 * pi * h3^(2 * rr + 1))

	# h2
	indh2 <- calculate_indh(2, th3)
	h2 <- hgrid[indh2]

	# theta 2
	phi_W <- ComputePhiEmp(W, tt/h2)
	th2 <- sum(tt^(2 * rr) * phi_W$norm^2 * phiKsq / phiUsq[, indh2])
	th2 <- th2 * deltat / (2 * pi * h2^(2 * rr + 1))

	# h1

	term1 <- hgrid^4 * mu_K2^2 * th2 / 4
	term2 <- kronecker(matrix(1,1,lh), phiKsq) / phiUsq
	term2 <- apply(term2, 2, sum) * deltat / (2 * pi * n * hgrid)
	AMISE <- term1 + term2
	indh1 <- which.min(AMISE)
	
	hgrid[indh1]

}