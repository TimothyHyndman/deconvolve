PI_deconvUknownth4het <- function(W, varX, phiUkvec, kernel_type){
#compute 2-stage plug-in bandwidth for heteroscedastic kerndel deconvolution estimator as in:
#Delaigle, A. and Meister, A. (2008). Density estimation with heteroscedastic error. Bernoulli, 14, 562-579
# !!! This code is only valid for a kernel of order 2 !!!!

#phiUveck: vector of n functions that give the characteristic functiona of the n errors. Produce this vector by c(func1,func2,...,funcn) where each funcj is a function of tt



	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	mu_K2 <- kernel_list$muk2
	RK <- kernel_list$rk
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]
	n <- length(W)



	W=as.vector(W)

	#Grid of h on which to search for a solution. If you did not find a minimum on that grid you can redefine it
	maxh=(max(W)-min(W))/10
	#normal reference bandwidth of the naive KDE estimator (estimator that ignores the errors) using the same kernel as above
	hnaive=((8*sqrt(pi)*RK/3/mu_K2^2)^0.2)*sqrt(stats::var(W))*n^(-1/5)
	hgrid=seq(hnaive/3,maxh,(maxh-hnaive/3)/100)
	lh = length(hgrid)
	dim(hgrid)=c(1,lh)


	#Estimator of the standard deviation of X
	stdevx = max(sqrt(varX),1/n)

	#Quantities that will be needed several times in the computations below
	toverh=tt%*%(1/hgrid)
	phiKsq=(phiK(tt))^2
	phiKsq=as.vector(phiKsq)

	# Convert vector of functions to single function ---------------------------
	phiUk <- function(tt,k) {
		phiUkvec[[k]](tt)
	}

	#Compute sum over k of phiU_k(tt/h)^2
	phiUsq=matrix(0,length(tt),length(hgrid))
	for (k in seq_len(n))
		{
		matphiU=phiUk(toverh,k)
		phiUsq=phiUsq+matphiU^2
		}


	calculate_indh <- function(rr, th) {
		term1 <- -hgrid^2 * mu_K2 * th
		term2 <- kronecker(matrix(1, 1, lh), tt^(2*rr) * phiKsq) / phiUsq
		term2 <- colSums(term2) * deltat  / (2 * pi * hgrid^(2 * rr + 1))
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
		
		#Estimate empirical characteristic function of W at t/h
		OO <- outer(tt / h, W)
		
		# Compute phiU(-t/h) -- since phiU is symmetric, this is the same as phiU(t/h)
		matphiU=OO
		for (k in seq_len(n))
			{matphiU[,k]=phiUk(tt/h,k)}
		
		#Estimate empirical characteristic function of X at t/h
		rehatphiX=apply(cos(OO)*matphiU,1,sum)
		imhatphiX=apply(sin(OO)*matphiU,1,sum)
		
		rm(OO)
		
		#Compute th
		normhatphiX2=(rehatphiX^2+imhatphiX^2)/(phiUsq[,indh])^2
		th = sum(tt^(2*rr) * normhatphiX2 * phiKsq)
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
	term2=kronecker(matrix(1,1,lh),phiKsq)/phiUsq
	term2=apply(term2,2,sum)*deltat/(2*pi*hgrid)
	AMISE=term1+term2

	indh=which.min(AMISE)
	hPI = hgrid[indh]
	return (hPI)
	
}

