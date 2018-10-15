decon_err_het_replicates <- function(xx, W1, W2,deno_U, kernel_type, h, rescale) {

	kernel_list <- kernel(kernel_type)
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]
	phiK <- kernel_list$phik

	OO_sum <- outer(tt / h, (W1 + W2)/2)
	num <- rowSums(exp(1i*OO_sum))
	den <- deno_U(tt / h)

	phi_X <- num/den
	rehatphiX <- Re(phi_X)
	imhatphiX <- Im(phi_X)

	# Matrix of size length(t) x length(xx)
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