DeconErrSymPmfToPdf <- function(X.pmf, W, phi.W, xx, phiK = NULL, muK2 = 6, 
								t = seq(-1, 1, 2e-04), rescale = FALSE, 
								h = NULL){

	theta <- X.pmf$support
	p <- X.pmf$probweights
	dt <- t[2] - t[1]
	tt <- phi.W$t.values
	
	if(is.null(phiK)){
		phiK <- phiK2
	}

	# Estimate Var(U) ----------------------------------------------------------
	tt.BB.length <- 200		# Use a finer grid than tt
	tt.BB <- seq(tt[1], tt[length(tt)], length.out = tt.BB.length)

	hat.var.U <- estimate_var_u(W, tt.BB, theta, p)

	# Estimate PhiX and PhiU ---------------------------------------------------
	phi.X <- ComputePhiPmf(theta, p, tt)
	phi.U <- phi.W$norm / Mod(phi.X)

	# Find Plug-In Bandwidth ---------------------------------------------------
	if (is.null(h)) {
		h.PIc <- PI_DeconvUEstTh4(W, phi.U, hat.var.U, tt, phiK, muK2, t)	
	} else {
		h.PIc <- h
	}
	

	# --------------------------------------------------------------------------
	phi.U.PI <- PhiUSpline(t/h.PIc, hat.var.U, phi.U, tt)
	phi.W.PI <- ComputePhiEmp(W, t/h.PIc)

	phi.X.re <- phi.W.PI$re / phi.U.PI
	phi.X.im <- phi.W.PI$im / phi.U.PI

	xt <- outer(t/h.PIc, xx)

	fX <- cos(xt) * matrix( rep(phi.X.re, length(xx)), ncol = length(xx)) + 
	sin(xt) * matrix( rep(phi.X.im, length(xx)), ncol = length(xx))

	fX <- colSums(fX * matrix( rep(phiK(t), length(xx)), ncol = length(xx))) / 
	(2 * pi) * dt / h.PIc

	fX[fX < 0] <- 0

	if (rescale) {
		dx <- xx[2] - xx[1]
		fX <- fX / sum(fX) / dx	
	}
	
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

	t_vec <- 		  tt.BB[ phi.U.BB >= 0.95 ]
	phi.U.t_vec <- phi.U.BB[ phi.U.BB >= 0.95 ]

	pp <- stats::lm(phi.U.t_vec ~ stats::poly(t_vec, 2, raw = T))
	hat.var.U <- -2 * pp$coefficients[[3]]
}