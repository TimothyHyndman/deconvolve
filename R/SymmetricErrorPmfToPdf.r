# X.pmf is a list with parts "support" and "probweights"
#' @export
# Just for testing
SymmetricErrorPmfToPdf <- function(theta, p, W, phi.W){

	# theta <- X.pmf$support
	# p <- X.pmf$probweights

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
	xt <- outerop(t/h.PIc, xx, "*")

	fX <- cos(xt) * matrix( rep(phi.X.re, length(xx)), ncol = length(xx)) + 
	sin(xt) * matrix( rep(phi.X.im, length(xx)), ncol = length(xx))

	fX <- colSums(fX * matrix( rep(PhiK(t), length(xx)), ncol = length(xx))) / 
	(2 * pi) * dt / h.PIc

	fX[fX < 0] <- 0
	fX <- fX / sum(fX)

	return(list("x" = xx, "y" = fX))

}