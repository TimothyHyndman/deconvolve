decon_err_het_replicates <- function(xx, W1, W2, kernel_type, h, rescale) {

	kernel_list <- kernel(kernel_type)
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]
	phiK <- kernel_list$phik

	OO_diff <- outer(tt / h, (W1 - W2)/2)
	OO_sum <- outer(tt / h, (W1 + W2)/2)

	num <- rowSums(exp(1i*OO_sum))
	den <- rowSums(cos(OO_diff))

	# Use laplace where den is too small as in Delaigle, Hall and Meister (2008) but replacing
	# the interval A at page 678 of that paper by the cutoff used by Camirand, Carroll and Delaigle (2018),
	# which is a refined version of the cutoff proposed by Delaigle and Hall (2016) 
	t_cutoff <- find_t_cutoff(den, tt)
	lap_ind <- (tt < t_cutoff) | (tt > t_cutoff)
	sd_U <- sqrt(stats::var(W1 - W2)/2)
	phiU_lap <- function(t){
		1 / (1 + sd_U^2 / 2 * t^2)
	}
	den[lap_ind] <- phiU_lap(tt[lap_ind]/2/h)^2


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