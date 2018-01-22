plugin_bandwidth <- function(W, phi_U, stdevx, phi_K, mu_K2, RK) {
	  
	# Grid of h on which to search for a solution
	maxh <- (max(W) - min(W)) / 10
	hnaive <- ((8 * sqrt(pi) * RK/3/mu_K2^2)^0.2) * sqrt(stats::var(W)) * n^(-1/5)
	hgrid <- seq(hnaive / 3, maxh, length.out = 100)


	# theta 4
	th4 <- (st.dev.X^(-9) * factorial(8) ) / (2^9 * factorial(4) * sqrt(pi))

	# h3
	rr <- 3
	term1 <- -hgrid^2 * muK2 * th4
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