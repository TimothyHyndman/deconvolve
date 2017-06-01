#' Create Probability Density Function from Probability Mass Function
#' 
#' Takes the output of DeconErrSymPmf along with the contaminated data, W, and
#' uses a kernel estimator to find a probability density function for X.
#' 
#' Details here
#' 
#' @param X.pmf A list with parts "support" and "probweights"
#' @param W A vector of the contaminated data
#' @param phi.W A list with parts "complex", "re", "im", "norm" and "t.values" 
#' containing phi.W, Re(phi.W), Im(phi.W), Norm(phi.W) and the t values on which
#' they were calculated respectively.
#' 
#' @return A list with components:
#' \item{x}{The values on which \eqn{f_X} was calculated}
#' \item{y}{The density \eqn{f_X(x)} evaluated at each x}
#' 
#' @example man/examples/SymmetricError_eg.R
#' 
#' @export

DeconErrSymPmfToPdf <- function(X.pmf, W, phi.W){

	theta <- X.pmf$support
	p <- X.pmf$probweights

	# Estimate phi.X and phi.U
	tt <- phi.W$t.values
	# a <- phi.W
	

	#---------------------------------------------------------------------------
	# Estimate var(U): approximate phi.U by poly of degree 2, and estimate varU 
	# by -2 * second order coefficient
	#---------------------------------------------------------------------------

	tt.BB.length <- 200		# Use a finer grid than tt
	tt.BB <- seq(tt[1], tt[length(tt)], length.out = tt.BB.length)

	phi.W.BB <- ComputePhiEmp(W, tt.BB)
	phi.X.BB <- ComputePhiPmf(theta, p, tt.BB)
	phi.U.BB <- phi.W.BB$norm / Mod(phi.X.BB)

	t.vec <- 		  tt.BB[ phi.U.BB >= 0.95 ]
	phi.U.t.vec <- phi.U.BB[ phi.U.BB >= 0.95 ]

	pp <- stats::lm(phi.U.t.vec ~ poly(t.vec, 2, raw = T))
	hat.var.U <- -2 * pp$coefficients[[3]]

	# Find Plug-In Bandwidth
	phi.X <- ComputePhiPmf(theta, p, tt)
	phi.U <- phi.W$norm / Mod(phi.X)
	h.PIc <- PI_DeconvUEstTh4(W, phi.U, hat.var.U, tt)

	# Kernel stuff
	dt <- 0.0002
	t <- seq(-1,1, by = dt)
	PhiK <- function(t){
		(1 - t^2)^3
	}

	phi.U.PI <- PhiUSpline(t/h.PIc, hat.var.U, phi.U, tt)
	phi.W.PI <- ComputePhiEmp(W, t/h.PIc)

	phi.X.re <- phi.W.PI$re / phi.U.PI
	phi.X.im <- phi.W.PI$im / phi.U.PI

	xx.length <- 100
	xx <- seq(min(W), max(W), length.out = xx.length)
	dx <- xx[2] - xx[1]
	xt <- outerop(t/h.PIc, xx, "*")

	fX <- cos(xt) * matrix( rep(phi.X.re, length(xx)), ncol = length(xx)) + 
	sin(xt) * matrix( rep(phi.X.im, length(xx)), ncol = length(xx))

	fX <- colSums(fX * matrix( rep(PhiK(t), length(xx)), ncol = length(xx))) / 
	(2 * pi) * dt / h.PIc

	fX[fX < 0] <- 0
	fX <- fX / sum(fX) / dx

	return(list("x" = xx, "y" = fX))

}