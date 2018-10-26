PI_deconvUestth4het <- function(W1,W2,hnaive,stdevx,Den, kernel_type){

#compute 2-stage plug-in bandwidth for heteroscedastic kerndel deconvolution estimator as in:
#Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic error. Bernoulli, 14, 562-579
#in the case where the error distributions are estimated through replicates
# !!! This code is only valid for a kernel of order 2 !!!!
#Den=Den of the estimator in Delaigle and Meister (2008)
#in the case where the error distributions are unknown and estimated from replicates.


	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	mu_K2 <- kernel_list$muk2
	RK <- kernel_list$rk
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]
	n <- length(W1)

	#Grid of h on which to search for a solution. If you did not find a minimum on that grid you can redefine it
	maxh=(max(W1)-min(W1))/10
	lh <- 101
	hgrid <- seq(hnaive / 3, maxh, length.out = lh)

	#Quantities that will be needed several times in the computations below
	toverh <- tt %*% t(1 / hgrid)
	phiKsq <- phiK(tt)^2
	deno_th <- Den(toverh)^2

	calculate_indh <- function(rr, th) {
		term1 <- -hgrid^2 * mu_K2 * th
		term2 <- kronecker(matrix(1, 1, lh), tt^(2*rr) * phiKsq) / deno_th
		term2 <- n*colSums(term2) * deltat  / (2 * pi * hgrid^(2 * rr + 1))
		ABias2 <- (term1 + term2)^2

		indh <- which.min(ABias2)
		if (indh == 1) {
			warning(paste("Minimum of Abias2 for rr =", as.character(rr), "is 
				the first element of the grid of bandwidths. Consider enlarging 
				the grid."))
		}
		if(indh == lh){
			warning(paste("Minimum of Abias2 for rr =", as.character(rr), "is 
				the last element of the grid of bandwidths. Consider enlarging 
				the grid."))
		}

		indh
	}

	calculate_th <- function(h, indh) {
		
		phi_W <- ComputePhiEmp((W1+W2)/2, tt/h)
		th <- sum(tt^(2 * rr) * (n*phi_W$norm)^2 * phiKsq / deno_th[, indh])

		th * deltat / (2 * pi * h^(2 * rr + 1))
	}


	#---------------------------------------------------------------------------
	# Start with theta4 and iterate down to get h1
	#---------------------------------------------------------------------------

	# Estimate theta4 by normal reference method     

	th4 = stdevx^(-9)*105/(32*sqrt(pi)) 

	# h3

	rr=3
	indh3 <- calculate_indh(rr, th4)
	h3 <- hgrid[indh3]

	# theta 3
	th3 <- calculate_th(h3, indh3)

	# h2
	rr <- 2
	indh2 <- calculate_indh(rr, th3)
	h2 <- hgrid[indh2]

	# theta 2
	th2 <- calculate_th(h2, indh2)



	# ------------------------------------------------------------------------------------------------------
	# Finally, compute the bandwidth that minimises the AMISE of the deconvolution kernel density estimator
	# ------------------------------------------------------------------------------------------------------

	term1=hgrid^4*mu_K2^2*th2/4
	term2=kronecker(matrix(1,1,lh),phiKsq)/deno_th
	term2=apply(term2,2,sum)*deltat/(2*pi*hgrid)
	AMISE=term1+n*term2

	indh=which.min(AMISE)
	hPI = hgrid[indh]
	return (hPI)
	
}

