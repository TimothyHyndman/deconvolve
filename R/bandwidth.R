#' Bandwidth Selectors for Deconvolution Kernel Density Estimation
#' 
#' Description
#' 
#' Details
#' 
#' @export

bandwidth <- function(W, errortype, sigU, phiU, phiK = phiK2, muK2 = 6, 
					  RK = 1024 / 3003 / pi, tt = seq(-1, 1, 2e-04)){
	n <- length(W)
	deltat <- tt[2] - tt[1]

	# # Error Known
	# PI_deconvUknownth4
	# # Error Heteroscedastic
	# PI_deconvUknownth4het
	# # Error estimated
	# PI_DeconvUEstTh4
	# #SIMEX
	# hSIMEXUknown
	# #CV
	# CVdeconv

}