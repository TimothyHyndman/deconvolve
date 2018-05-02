#' @export

hom_deconvolve <- function(W1, 
						   W2 = NULL,
						   x = seq(min(W1), max(W1), length.out = 100),
						   h = NULL,
						   kernel_type = c("default", "normal", "sinc"),
						   rescale = FALSE,
						   pmf = FALSE,
						   m = 20) {

	kernel_type <- match.arg(kernel_type)

	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	muK2 <- kernel_list$muk2
	RK <- kernel_list$rk
	tt <- kernel_list$tt
	deltat <- tt[2] - tt[1]

	n <- length(W1)

	if (!is.null(W2)) {
		errors <- "rep"
	} else {
		errors <- "sym"
	}

	if (is.null(h) & !(errors == "sym")) {
			h <- bandwidth(W1, W2, kernel_type = kernel_type)
	}

	if (errors == "rep") {
		t_search <- tt/h
		phiU_splined <- function(t){
			replicates_phiU(t, W1, W2, t_search)
		}
		W <- c(W1, W2)
		output <- hom_deconvolve_U_known(W, phiU_splined, h, x, kernel_type, rescale)
		output <- list("x" = x, "pdf" = output$pdf, "W1" = W1, "W2" = W2)
	}

	if (errors == "sym") {
		W <- W1
		out <- DeconErrSymPmf(W, m, kernel_type)
		h <- NULL

		if (!pmf) {
			phi_W <- out$phi.W
			pdf <- DeconErrSymPmfToPdf(out, W, phi_W, x, kernel_type, rescale, h)
			output <- list("x" = x, "pdf" = pdf, "support" = out$support, 
						   "probweights" = out$probweights, "W" = W)
		} else {
			output <- list("support" = out$support,
						   "probweights" = out$probweights,
						   "W" = W)
		}
	}


	class(output) <- c("deconvolve", "list")
	output
}