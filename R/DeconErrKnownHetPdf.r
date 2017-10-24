DeconErrKnownHetPdf<-function(xx, W, h, phiUkvec, 
	rescale = FALSE, phiK = NULL, muK2 = 6, RK = 1024 / 3003 / pi, 
	tt = seq(-1, 1, 2e-04)){

	# Convert vector of functions to single function ---------------------------
	phiUk <- function(tt,k) {
		phiUkvec[[k]](tt)
	}

	#--------------------------------------------------------------------------#
	# Compute DKDE
	#--------------------------------------------------------------------------#

	if(is.null(phiK)){
		phiK <- phiK2
	}

	W <- as.vector(W)
	n <- length(W)
	deltat <- tt[2] - tt[1]

	# Default values of phiU(t)=characteristic function of the errors
	# If you want to consider another error type, simply replace phiU by the 
	# characteristic function of your error type

	OO <- outer(tt / h, W)

	# Compute phiU_k(-t/h) for each k -- since phiU_k is symmetric, this is the 
	# same as phiU_k(t/h)
	matphiU <- OO;
	for (k in 1:n)
		matphiU[, k] <- phiUk(tt / h, k)

	# Sum by rows of the matrix. This produces a vector of size equal to that of 
	# tt
	phiUsqth <- apply(matphiU^2, 1, sum)

	# Estimate real and imaginary parts of empirical characteristic function of 
	# W computed at tt/h, for each component of tt.
	# Results = vectors of size length(tt)

	rehatphiX <- apply(cos(OO) * matphiU, 1, sum) / phiUsqth
	imhatphiX <- apply(sin(OO) * matphiU, 1, sum) / phiUsqth

	# Matrix of size length(tt) x length(xx)
	xt <- outer(tt / h, xx)
	longx <- length(xx)

	# Compute the DKDE estimator
	fXdecUK <- cos(xt) * kronecker(matrix(1, 1, longx), rehatphiX) + sin(xt) * 
			   kronecker(matrix(1, 1, longx), imhatphiX)
	fXdecUK <- apply(fXdecUK * kronecker(matrix(1, 1, longx), phiK(tt)), 2, sum) / 
			   (2 * pi) * deltat / h
	fXdecUK[which(fXdecUK < 0)] <- 0

	if (rescale == 1){
		dx <- xx[2] - xx[1]
		fXdecUK <- fXdecUK / sum(fXdecUK) / dx
	}


	return (fXdecUK)
}
