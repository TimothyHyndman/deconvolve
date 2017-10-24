DeconErrKnownPdf<-function(xx, W, h, phiU, rescale = FALSE, 
	phiK = NULL, muK2 = 6, RK = 1024 / 3003 / pi, tt = seq(-1, 1, 2e-04)){

	#--------------------------------------------------------------------------#
	# Compute DKDE
	#--------------------------------------------------------------------------#
	if(is.null(phiK)){
		phiK <- phiK2
	}

	W <- as.vector(W)
	n <- length(W)
	deltat <- tt[2] - tt[1]

	OO <- outer(tt/h, W)
	phiUth <- phiU(tt/h)

	# Estimate real and imaginary parts of empirical characteristic function of 
	# W computed at tt/h, for each component of tt.
	rehatphiX <- rowSums(cos(OO)) / phiUth / n
	imhatphiX <- rowSums(sin(OO)) / phiUth / n

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