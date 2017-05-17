# X.pmf is a list with parts "support" and "probweights"
#' @export
# Just for testing
SymmetricErrorPmfToPdf <- function(X.pmf, W, phi.W){

	tt <- phi.W$t.values
	theta <- X.pmf$support
	p <- X.pmf$probweights

	# Estimate phi.X and phi.U
	phi.X <- ComputePhiPmf(theta, p, tt)
	phi.U <- phi.W$norm / Mod(phi.X)

	#---------------------------------------------------------------------------
	# Estimate var(U): approximate phi.U by poly of degree 2, and estimate varU 
	# by -2 * second order coefficient
	#---------------------------------------------------------------------------

	tt.BB.length <- 200		# Use a finer grid than tt
	tt.BB <- seq(tt[1], tt[length(tt)], length.out = tt.BB.length)

	phi.W.BB <- ComputePhiEmp(W, tt.BB)
	phi.X.BB <- ComputePhiPmf(theta, p, tt.BB)
	phi.U.BB <- phi.W.BB$norm / Mod(phi.X)

	t.vec <- 		  tt.BB[ phi.U.BB >= 0.95 ]
	phi.U.t.vec <- phi.U.BB[ phi.U.BB >= 0.95 ]

	pp <- stats::lm(phi.U.t.vec ~ poly(t.vec, 2, raw = T))
	hat.var.U <- -2 * pp$coefficients[[3]]

	#---------------------------------------------------------------------------
	# Compute density estimator as indicated in the paper
	#---------------------------------------------------------------------------

	t.limits <- c(min(tt), max(tt))
	# pp.phi.U <- spline(tt, hat.phi.U)	#No need to do this calculation here.
	# We'll do it later on in phiUspline or PI_deconvUestth4

	h.PIc <- PI_deconvUestth4(W, tlim, phi.U, hat.var.U, phi.U, tt)

	# fX.dec <- fXKernDec2()

	# fX.dec[fX.dec < 0] <- 0
	# fX.dec <- fX.dec / sum(fX.dec)

}