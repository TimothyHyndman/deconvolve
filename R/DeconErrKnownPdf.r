DeconErrKnownPdf <- function(xx, W, h, phiU, kernel_type, rescale){

	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]

	n <- length(W)

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


	fXdecUK
}
