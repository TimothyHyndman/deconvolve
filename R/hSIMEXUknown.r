#' @importFrom foreach %dopar%
#'
# Author: Aurore Delaigle
# Computes bandwidth h and ridge parameter rho using a version of the SIMEX
# method of
# Delaigle, A. and Hall, P. (2008). Using SIMEX for smoothing-parameter choice
# in errors-in-variables problems. JASA, 103, 280-287
#
# WARNING: these are not the codes used in the original paper. This is a
# simplified version of those codes.
#
# Use the function NWDecUknown to compute the regression estimator with this rho
# and this h
# outputs
# h: bandwidth
# rho: ridge parameter.

hSIMEXUknown <- function(W, Y, generate_U_star, sd_U, phiU, kernel_type, n_cores, 
						 seed){

	if (is.null(n_cores)){
		n_cores <- parallel::detectCores()
		n_cores <- max(n_cores - 1, 1)
	}
	cl <- parallel::makeCluster(n_cores)
	doParallel::registerDoParallel(cl)

	# --------------------------------------------------------------------------
	# Preliminary calculations and initialisation of functions
	# --------------------------------------------------------------------------
	n <- length(W)
	# number of bins used to compute CV in each SIMEX world
	nbin <- min(100, n)

	# Define a grid where to search for the SIMEX bandwidth. By default we take
	# [h/2,2h], where h=PI bandwidth for density estimation.
	# Increase the grid if too small
	
	sd_X <- max( !is.na(sqrt(stats::var(as.vector(W)) - sd_U^2)), 1/n)
	hPIfX <- plugin_bandwidth(W, phiU, sd_X, kernel_type)
	gridh <- seq(hPIfX / 2, 2 * hPIfX, length.out = 21)


	# Define a defaul grid where to search for rho.
	# Recall that rho prevents the denominator of the NW estimator from being
	# too small. In the SIMEX world, the denominator estimates the contaminated
	# density f_W
	# This is what motivates the default grid for rho used here.

	# Estimator of fW(q_{0.05}) and fW(q_{0.95}) using standard (error-free) KDE
	# and normal reference bandwidth, where q_{alpha} denotes the alpha
	# empirical quantile of the W_i's.
	# W <- as.vector(W)
	hW <- 1.06 * sqrt(stats::var(W)) * n^(-1 / 5)
	ab <- stats::quantile(W, probs = c(0.05, 0.95))
	xout <- outer(ab, W, "-")
	fWEF <- c(0, 0)
	fWEF[1] <- mean(stats::dnorm(xout[1, ], 0, hW))
	fWEF[2] <- mean(stats::dnorm(xout[2, ], 0, hW))
	gridrho <- min(fWEF) * seq(0.025, 4, 0.025)


	lh <- length(gridh)
	lrho <- length(gridrho)

	#---------------------------------------------------------------------------
	# Step 1: find the ridge parameter using only the first level of SIMEX
	#---------------------------------------------------------------------------

	# Bin the W data to speed up the computations
	midbin <- unlist(BinData(W, nbin)[1])
	indbin <- matrix(unlist(BinData(W, nbin)[2]), nrow = n)

	BB <- 20
	bb <- NULL	# This is here to remove a note from check
	outcome_SIMEX1 <- foreach::foreach(bb = 1:BB, .packages = c("stats")) %dopar% {
	# for (bb in 1:BB) {
	    #set seed
	    if (!is.null(seed)) {
	        set.seed(seed + bb)
	    } else {
	    	set.seed(bb+1000)
		}
		CVrho <- matrix(0, lh, lrho)
		# Generate SIMEX data Wstar
		Wstar <- W + generate_U_star()

		# For each h in the grid of h-candidates, compute the CV criterion for
		# the data Wstar (this will automatically consider all rho candiates)
		for (kh in 1:lh){
			h <- gridh[kh]
			CVrho[kh, ] <- NWDecridgeL1OCUknown(n, Wstar, Y,
						   phiU, h, gridrho, midbin, indbin, nbin, kernel_type)
		}
		CVrho
	}
	# outcome_SIMEX1 <- CVrho

	CVrho <- matrix(0, lh, lrho)
	for (i in 1:BB)
	{
		CVrho <- CVrho + outcome_SIMEX1[[i]]
	}


	# find which pair of (h,rho) minimizes CV
	minCV <- which.min(CVrho)
	indh <-arrayInd(minCV, dim(CVrho))[1]
	indrho <-arrayInd(minCV, dim(CVrho))[2]
	# Rigdge parameter
	rho <- gridrho[indrho]

	# h from SIMEX level 1
	h1 <- gridh[indh]


	#---------------------------------------------------------------------------
	# Step 2: Keep rho fixed and find h SIMEX
	#---------------------------------------------------------------------------


	outcome_SIMEX2 <- foreach::foreach(bb = 1:BB, .packages = c("stats")) %dopar% {

	    if (!is.null(seed)){
	        set.seed(seed + bb + BB)
	    }else{
	        set.seed(bb + 1000 + BB)}

		CVhstar_tmp <- 0 * gridh
		# Generate SIMEX data Wstar2
		Wstar <- W + generate_U_star()
		Wstar2 <- Wstar + generate_U_star()

		# Bin the Wstar data to speed up the computations
		midbin <- unlist(BinData(Wstar, nbin)[1])
		indbin <- unlist(BinData(Wstar, nbin)[2])
		# Compute CV for each h in the grid, using the ridge parameter rho found
		# above
		for (kh in 1:lh){
			h <- gridh[kh]
			CVhstar_tmp[kh] <- NWDecridgeL1OCUknown(n, Wstar2, Y,
							   phiU, h, rho, midbin, indbin, nbin, kernel_type)
		}
		CVhstar_tmp
	}

	parallel::stopCluster(cl)

	CVhstar <- 0 * gridh
	for (i in 1:BB){
		CVhstar <- CVhstar + outcome_SIMEX2[[i]]
	}
	indh <- which.min(CVhstar)


	# h from SIMEX level 2
	h2 <- gridh[indh]

	# Finally deduce h SIMEX
	h <- h1^2 / h2
	outcome <- list(h, rho)
	names(outcome) <- c('h', 'rho')

	outcome
}

BinData <- function(W, nbin){
	# Author: Aurore Delaigle
	# This program bins the data W into nbins

	n <- length(W)

	# Compute extremities of the bins
	ab <- stats::quantile(W, c(0.05, 0.95))

	# Bin widths
	delta <- (ab[2] - ab[1]) / nbin

	# Bin centers
	midbin <- ab[1] + delta / 2 + delta * (seq(0, (nbin - 1)))

	# Bin left and right extremities
	Abin <- midbin - delta / 2
	Bbin <- midbin + delta / 2

	# Find in which bin each observation lies
	Wmat <- kronecker(matrix(1, 1, nbin), W)
	Amat <- matrix(rep(Abin, n), nrow = n, byrow = TRUE)
	Bmat <- matrix(rep(Bbin, n), nrow = n, byrow = TRUE)
	indice <- matrix(rep(seq(1, nbin), n), nrow = n, byrow = TRUE)
	indice <- apply(indice * ((Wmat > Amat) & (Wmat <= Bmat)), 1, sum)

	# Put those beyond the extremities at the extremities
	indice[which(W <= Abin[1])] <- 1
	indice[which(W > Bbin[nbin])] <- nbin

	list2 <- list(midbin, indice)

	list2
}
