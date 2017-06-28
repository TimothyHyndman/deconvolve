#' @export

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

# W: vector of contaminated data W_1,...,W_n
# Y: vector of data Y_1,...,Y_n
# h: bandwidth
# 
# errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other 
# error distributions, simply redefine phiU below 
# sigU: parameter of Laplace or normal errors used only to define phiU.
# rho: ridge parameter. 

hSIMEXUknown <- function(W, Y, errortype, sigU){
  
	W <- as.vector(W)
	n <- length(W)

	# --------------------------------------------------------
	# Preliminary calculations and initialisation of functions
	# --------------------------------------------------------

	# Default values of phiU(t)=characteristic function of the errors
	# If you want to consider another error type, simply replace phiU by the 
	# characteristic function of your error type
	if (errortype == "Lap") {
		phiU <- function(t) {
			1 / (1 + sigU^2 * t^2)
		}
		varU <- 2 * sigU^2
	}
	if (errortype == "norm") {
		phiU <- function(t) {
			exp(-sigU^2 * t^2 / 2)
		}
		varU <- sigU^2
	}

	# phiK: Fourier transform of the kernel K. You can change this if you wish 
	# to use another kernel but make sure 
	# you change the range of t-values, which should correspond to the support 
	# of phiK
	phiK <- function(t) {
		(1 - t^2)^3
	}


	# Range of t-values (must correspond to the domain of phiK)
	dt <- .0002
	t <- seq(-1, 1, dt)
	dim(t) <- c(length(t), 1)




	dim(W) <- c(1, n)
	dim(Y) <- c(1, n)

	# number of bins used to compute CV in each SIMEX world
	nbin <- min(100, n)	

	# Number of SIMEX samples
	BB <- 20
	# BB <- 1 	#For testing to make it faster

	# Define a grid where to search for the SIMEX bandwidth. By default we take 
	# [h/2,2h], where h=PI bandwidth for density estimation.
	# Increase the gird if too small
	hPIfX <- PI_deconvUknownth4(n, W, errortype, sigU)
	a <- hPIfX / 2
	b <- 2 * hPIfX
	gridh <- seq(a, b, (b - a) / 20)


	# Define a defaul grid where to search for rho. 
	# Recall that rho prevents the denominator of the NW estimator from being 
	# too small. In the SIMEX world, the denominator estimates the contaminated 
	# density f_W
	# This is what motivates the default grid for rho used here.

	# Estimator of fW(q_{0.05}) and fW(q_{0.95}) using standard (error-free) KDE 
	# and normal reference bandwidth, where q_{alpha} denotes the alpha 
	# empirical quantile of the W_i's.
	W <- as.vector(W)
	hW <- 1.06 * sqrt(stats::var(W)) * n^(-1 / 5)
	ab <- stats::quantile(W, probs = c(0.05, 0.95))
	xout <- outerop(ab, W, "-")
	fWEF <- c(0, 0)
	fWEF[1] <- mean(stats::dnorm(xout[1, ], 0, hW))
	fWEF[2] <- mean(stats::dnorm(xout[2, ], 0, hW))
	gridrho <- min(fWEF) * seq(0.025, 4, 0.025)


	lh <- length(gridh)
	lrho <- length(gridrho)
	CVrho <- matrix(0, lh, lrho)



	#---------------------------------------------------------------------
	# Step 1: find the ridge parameter using only the first level of SIMEX
	#---------------------------------------------------------------------

	# Bin the W data to speed up the computations
	midbin <- unlist(BinData(W, nbin)[1])
	indbin <- matrix(unlist(BinData(W, nbin)[2]), nrow = n)

	for (bb in 1:BB) {
		# Generate SIMEX data Wstar
		if (errortype == "Lap") {
			Wstar <- W + rlap(sigU, 1, n)
		}
		if (errortype == "norm") {
			Wstar <- W + stats::rnorm(n, 0, sigU)
		}

		# For each h in the grid of h-candidates, compute the CV criterion for 
		# the data Wstar (this will automatically consider all rho candiates)
		for (kh in 1:lh){
			h <- gridh[kh]
			CVrho[kh, ] <- CVrho[kh, ] + NWDecridgeL1OCUknown(n, Wstar, Y, 
				errortype, sigU, h, gridrho, midbin, indbin, nbin)
		}
	}


	# find which pair of (h,rho) minimizes CV
	minCV <- which.min(CVrho)
	indh <-arrayInd(minCV, dim(CVrho))[1]
	indrho <-arrayInd(minCV, dim(CVrho))[2]
	# Rigdge parameter
	rho <- gridrho[indrho]

	# h from SIMEX level 1
	h1 <- gridh[indh]


	#----------------------------------------
	# Step 2: Keep rho fixed and find h SIMEX 
	#----------------------------------------

	CVhstar <- 0 * gridh

	for (bb in 1:BB){    
	# Generate SIMEX data Wstar2
		if (errortype == "Lap"){
			Wstar <- W + rlap(sigU, 1, n)
			Wstar2 <- Wstar + rlap(sigU, 1, n)
		}
		if (errortype == "norm") {
			Wstar <- W + stats::rnorm(n, 0, sigU)
			Wstar2 <- Wstar + stats::rnorm(n, 0, sigU)
		}

		# Bin the Wstar data to speed up the computations
		midbin <- unlist(BinData(Wstar, nbin)[1])
		indbin <- unlist(BinData(Wstar, nbin)[2])
		# Compute CV for each h in the grid, using the ridge parameter rho found 
		# above
		for (kh in 1:lh){
			h <- gridh[kh]
			CVhstar[kh] <- CVhstar[kh] + NWDecridgeL1OCUknown(n, Wstar2, Y, 
				errortype, sigU, h, rho, midbin, indbin, nbin)
		}
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
	dim(W) = c(n, 1)

	# Compute extremities of the bins
	ab = stats::quantile(W, c(0.05, 0.95))

	# Bin widths
	delta = (ab[2] - ab[1]) / nbin

	# Bin centers
	midbin = ab[1] + delta / 2 + delta * (seq(0, (nbin - 1)))

	# Bin left and right extremities
	Abin = midbin - delta / 2
	Bbin = midbin + delta / 2

	# Find in which bin each observation lies
	Wmat = kronecker(matrix(1, 1, nbin), W)
	Amat = matrix(rep(Abin, n), nrow = n, byrow = T)
	Bmat = matrix(rep(Bbin, n), nrow = n, byrow = T)
	indice = matrix(rep(seq(1, nbin), n), nrow = n, byrow = T)
	indice = apply(indice * ((Wmat > Amat) & (Wmat <= Bmat)), 1, sum)

	# Put those beyond the extremities at the extremities
	indice[which(W <= Abin[1])] = 1
	indice[which(W > Bbin[nbin])] = nbin

	list2 <- list(midbin, indice)
	
	list2
}