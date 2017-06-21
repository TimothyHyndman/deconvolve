#' Deconvolution
#' 
#' Description
#' 
#' Details
#' 
#' @section Problems:
#' Currently doesn't take h into account for symmetric deconvolution
#' 
#' @export

deconvolve <- function(W, xx, errortype = NULL, sigU = NULL, phiU = NULL, 
					   h = NULL, rescale = FALSE, phiK = NULL, muK2 = 6, 
					   RK = 1024 / 3003 / pi, tt = seq(-1, 1, 2e-04)){

	# Decide on type of deconvolution ------------------------------------------
	if (is.null(errortype) & is.null(phiU)) {
		decon_type <- "symmetric"
	} else if (length(sigU) > 1  | length(phiU) > 1){
		decon_type <- "heteroscedastic"
	} else {
		decon_type <- "known"
	}

	# Calculate Bandwidth if not supplied --------------------------------------
	if (is.null(h) & (decon_type == "symmetric") == FALSE) {
		if (is.null(phiU)) {
			h <- bandwidth(W, errortype, sigU, phiK = phiK, muK2 = muK2, 
					   RK = RK, tt = tt)
		} else {
			h <- bandwidth(W, phiU = phiU, phiK = phiK, muK2 = muK2, 
					   RK = RK, tt = tt)
		}
	}
	
	# Use default PhiK if not supplied -----------------------------------------
	if(is.null(phiK)){
		phiK <- phiK2
	}

	# Perform appropriate deconvolution ----------------------------------------
	if (decon_type == "known"){
		output <- DeconErrKnownPdf(xx, W, h, errortype, sigU, phiU, rescale, 
								   phiK, muK2, RK, tt)
	}

	if (decon_type == "heteroscedastic"){
		output <- DeconErrKnownHetPdf(xx, W, h, errortype, sigU, phiU, rescale, 
								   phiK, muK2, RK, tt)
	}

	if (decon_type == "symmetric") {
		out <- DeconErrSymPmf(W)
		phi.W <- out$phi.W
		output <- DeconErrSymPmfToPdf(out, W, phi.W, xx, phiK, muK2, tt)
	}

	# Output PDF ---------------------------------------------------------------
	output
}
