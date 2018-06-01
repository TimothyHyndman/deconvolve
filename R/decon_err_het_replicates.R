decon_err_het_replicates <- function(xx, W1, W2, kernel_type, h, rescale) {

	kernel_list <- kernel(kernel_type)
	t <- kernel_list$tt
	deltat <- t[2] - t[1]
	phiK <- kernel_list$phik

	OO_diff <- outer(t / h, (W1 - W2)/2)
	OO_sum <- outer(t / h, (W1 + W2)/2)

	num <- rowSums(exp(1i*OO_sum))
	den <- rowSums(exp(1i*OO_diff))

	# Use laplace where den is too small
	t_cutoff <- find_t_cutoff(Re(den), t)		#Re() isn't meant to be there...
												#How to fix with t_cutoff...?
	lap_ind <- (t < t_cutoff) | (t > t_cutoff)
	sd_U <- sqrt(stats::var(W1 - W2)/2)
	phiU_lap <- function(t){
		1 / (1 + sd_U^2 / 2 * t^2)
	}
	den[lap_ind] <- phiU_lap(t[lap_ind]/2)


	phi_X <- num/den
	rehatphiX <- Re(phi_X)
	imhatphiX <- Im(phi_X)

	# Matrix of size length(t) x length(xx)
	xt <- outer(t / h, xx)
	longx <- length(xx)

	# Compute the DKDE estimator
	fXdecUK <- cos(xt) * kronecker(matrix(1, 1, longx), rehatphiX) + sin(xt) * 
			   kronecker(matrix(1, 1, longx), imhatphiX)
	fXdecUK <- apply(fXdecUK * kronecker(matrix(1, 1, longx), phiK(t)), 2, sum) / 
			   (2 * pi) * deltat / h
	fXdecUK[which(fXdecUK < 0)] <- 0

	if (rescale == 1){
		dx <- xx[2] - xx[1]
		fXdecUK <- fXdecUK / sum(fXdecUK) / dx
	}

	fXdecUK
}