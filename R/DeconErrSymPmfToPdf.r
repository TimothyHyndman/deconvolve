#' Create Probability Density Function from Probability Mass Function
#' 
#' Takes the output of DeconErrSymPmf along with the contaminated data, W, and
#' uses a kernel estimator to find an estimate for the probability density 
#' function of X.
#' 
#' PUT DETAILS HERE
#' 
#' @param X.pmf A list with parts "support" and "probweights"
#' @param W A vector of the contaminated data
#' @param phi.W A list with parts "complex", "re", "im", "norm" and "t.values" 
#' containing phi.W, Re(phi.W), Im(phi.W), Norm(phi.W) and the t values on which
#' they were calculated respectively.
#' @param xx A vector of x values on which to compute the density.
#' 
#' @return A vector containing the density \eqn{f_X(x)} evaluated at each x in 
#' \code{xx}
#' 
#' @example man/examples/SymmetricError_eg.R
#' 
#' @section Authors:
#' Aurore Delaigle
#' 
#' @section References:
#' Delaigle, A. and Hall, P. (2016). Methodology for non-parametric 
#' deconvolution when the error distribution is unknown. \emph{Journal of the 
#' Royal Statistical Society: Series B (Statistical Methodology)}, 78, 1, 
#' 231-252.
#' 
#' @export

DeconErrSymPmfToPdf <- function(X.pmf, W, phi.W, xx, PhiK, muK2, t){

	theta <- X.pmf$support
	p <- X.pmf$probweights
	dt <- t[2] - t[1]

	# Estimate Var(U) ----------------------------------------------------------
	tt.BB.length <- 200		# Use a finer grid than tt
	tt.BB <- seq(tt[1], tt[length(tt)], length.out = tt.BB.length)

	hat.var.U <- estimate_var_u(W, tt.BB, theta, p)

	# Estimate PhiX and PhiU ---------------------------------------------------
	tt <- phi.W$t.values
	phi.X <- ComputePhiPmf(theta, p, tt)
	phi.U <- phi.W$norm / Mod(phi.X)

	# Find Plug-In Bandwidth ---------------------------------------------------
	h.PIc <- PI_DeconvUEstTh4(W, phi.U, hat.var.U, tt, PhiK, muK2, t)

	# --------------------------------------------------------------------------
	phi.U.PI <- PhiUSpline(t/h.PIc, hat.var.U, phi.U, tt)
	phi.W.PI <- ComputePhiEmp(W, t/h.PIc)

	phi.X.re <- phi.W.PI$re / phi.U.PI
	phi.X.im <- phi.W.PI$im / phi.U.PI

	dx <- xx[2] - xx[1]
	xt <- outerop(t/h.PIc, xx, "*")

	fX <- cos(xt) * matrix( rep(phi.X.re, length(xx)), ncol = length(xx)) + 
	sin(xt) * matrix( rep(phi.X.im, length(xx)), ncol = length(xx))

	fX <- colSums(fX * matrix( rep(PhiK(t), length(xx)), ncol = length(xx))) / 
	(2 * pi) * dt / h.PIc

	fX[fX < 0] <- 0
	fX <- fX / sum(fX) / dx

	fX
}

estimate_var_u <- function(W, tt.BB, theta, p){
	#---------------------------------------------------------------------------
	# Estimate var(U): approximate phi.U by poly of degree 2, and estimate varU 
	# by -2 * second order coefficient
	#---------------------------------------------------------------------------
	phi.W.BB <- ComputePhiEmp(W, tt.BB)
	phi.X.BB <- ComputePhiPmf(theta, p, tt.BB)
	phi.U.BB <- phi.W.BB$norm / Mod(phi.X.BB)

	t.vec <- 		  tt.BB[ phi.U.BB >= 0.95 ]
	phi.U.t.vec <- phi.U.BB[ phi.U.BB >= 0.95 ]

	pp <- stats::lm(phi.U.t.vec ~ stats::poly(t.vec, 2, raw = T))
	hat.var.U <- -2 * pp$coefficients[[3]]
}