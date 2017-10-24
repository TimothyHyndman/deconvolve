#' @export
# Author: Aurore Delaigle
# Compute the measurement error version of the Nadaraya-Watson regression 
# estimator
# Goal: estimate m where Y=m(X)+epsilon,  and we observe data on (W, Y),  where 
# W=X+U.
# See Fan,  J.,  and Truong,  Y. K. (1993),  Nonparametric Regression With 
# Errors in Variables,  The Annals of Statistics,  21,  19001925

# xx: vector of x-values where to compute the regression estimator
# W: vector of contaminated data W_1, ..., W_n
# Y: vector of data Y_1, ..., Y_n
# h: bandwidth
# 
# 
# errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other 
# error distributions,  simply redefine phiU below 
# sigU: parameter of Laplace or normal errors used only to define phiU.
# 
# rho: ridge parameter. See
# Delaigle,  A. and Hall,  P. (2008). Using SIMEX for smoothing-parameter choice 
# in errors-in-variables problems.  JASA,  103,  280-287 



#  -----------------------------------------------------------------------------
# 								WARNINGS:
#  -----------------------------------------------------------------------------
# 
# The range of t-values -1 and 1 correspond to the support of phiK. 
# If you change phiK and take a kernel for which phiK is not supported on 
# [-1, 1] you have to change -1 and 1 accordingly.
# 
# The phiK here must match the phiK used to compute the bandwidth (SIMEX or 
# other).
# 
# The estimator can also be computed using the Fast Fourier Transform,  which 
# is faster,  but more complex. 
# See Delaigle,  A. and Gijbels,  I. (2007). Frequent problems in calculating 
# integrals and optimizing objective functions: a case study in density 
# deconvolution.   Statistics and Computing,   17,   349 - 355
# However if the grid of t-values is fine enough,  the estimator can simply be 
# computed like here without having problems with oscillations.
# 
#  -----------------------------------------------------------------------------

NWDecUknown <- function(n, xx, W, Y, sigU, phiU, h, rho) {

	# --------------------------------------------------------
	# Preliminary calculations and initialisation of functions
	# --------------------------------------------------------


	# Default values of phiU(t)=characteristic function of the errors
	# If you want to consider another error type,  simply replace phiU by the 
	# characteristic function of your error type
	W <- as.vector(W)




	# phiK: Fourier transform of the kernel K. You can change this if you wish 
	# to use another kernel but make sure 
	# you change the range of t-values,  which should correspond to the support 
	# of phiK

	phiK <- function(t) {
		(1 - t^2)^3
	}


	# Range of t-values (must correspond to the domain of phiK)
	dt <- .0002
	t <- seq(-1, 1, dt)
	longt <- length(t)
	# dim(t) <- c(length(t), 1)



	# Compute the empirical characteristic function of W (times n) at t/h: 
	# \hat\phi_W(t/h)
	OO <- t(outer(t/h, t(W)))
	csO <- cos(OO)
	snO <- sin(OO)
	rm(OO)

	rehatphiW <- apply(csO, 2, sum)
	imhatphiW <- apply(snO, 2, sum)


	# Compute \sum_j Y_j e^{itW_j/h}
	dim(Y) <- c(1, n)
	renum <- Y %*% csO
	imnum <- Y %*% snO


	# Compute numerator and denominator of the estimator separately
	# Numerator: real part of 
	# (2*pi*h)^(-1) \int e^{-itx/h} n^(-1)\sum_j Y_j e^{itW_j/h} \phi_K(t)/\phi_U(t/h) dt
	# 
	# Denominator: real part of 
	# (2*pi*h)^(-1) \int e^{-itx/h} \hat\phi_W(t/h) \phi_K(t)/\phi_U(t/h) dt


	xt <- outer(t / h, t(xx))
	cxt <- cos(xt)
	sxt <- sin(xt)
	rm(xt)

	phiUth <- phiU(t / h)
	matphiKU <- phiK(t) / phiUth
	dim(matphiKU) <- c(1, longt)

	Den <- (rehatphiW * matphiKU) %*% cxt + (imhatphiW * matphiKU) %*% sxt
	Num <- (renum * matphiKU) %*% cxt + (imnum * matphiKU) %*% sxt
	Num <- Num * (dt / (n * h * 2 * pi))
	Den <- Den * (dt / (n * h * 2 * pi))


	# If denomintor is too small,  replace it by the ridge rho
	dd <- Den
	dd[which(dd < rho)] <- rho

	# Finally obtain the regression estimator
	y <- Num / dd

	as.vector(y)
}

