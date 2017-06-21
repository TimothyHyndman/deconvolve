#' Bandwidth Selectors for Deconvolution Kernel Density Estimation
#' 
#' Description
#' 
#' The function \code{bandwidth} chooses between three different algorithms
#' to calculate a bandwidth depending on the inputs.
#' 
#' @section Problems:
#' 
#' 
#' @export

bandwidth <- function(W, errortype, sigU, phiU, varX = NULL, algorithm = "PI", 
					  phiK = NULL, muK2 = 6, RK = 1024 / 3003 / pi, 
					  tt = seq(-1, 1, 2e-04)){
	
	n <- length(W)
	deltat <- tt[2] - tt[1]

	if(is.null(phiK)){
		phiK <- phiK2
	}

	# Determine Error Type Provided --------------------------------------------

	if (missing(phiU) & missing(sigU)) {
		errors <- "est"
	} else if (missing(phiU)) {
		if (length(sigU) > 1){
			errors <- "het"
			if ((length(sigU) == length(W)) == FALSE) {
				stop("sigU must be either length 1 for homoscedastic errors or have the same length as W for heteroscedastic errors.")
			}
		} else {
			errors <- "hom"
		}
	} else {
		if (length(phiU) > 1){
			errors <- "het"
			if ((length(phiU) == length(W)) == FALSE) {
				stop("phiU must be either length 1 for homoscedastic errors or have the same length as W for heteroscedastic errors.")
			}
		} else {
			errors <- "hom"
		}
	}

	# Check inputs -------------------------------------------------------------

	if ((algorithm == "CV" | algorithm == "PI") == FALSE) {
		stop("algorithm must be one of: 'PI', or 'CV'.")
	}

	if ((errors == "est") == FALSE) {
		if ((errortype == "norm" | errortype == "Lap") == FALSE) {
			stop("errortype must be one of: 'norm', or 'Lap'.")
		}
	}
	
	if (algorithm == "CV") {
		if (errors == "het") {
			stop("Algorithm type 'CV' can only be used with homoscedastic 
				 errors.")
		}
		if (errors == "est") {
			stop("You must supply error for algorithm 'CV'.")
		}
	}

	if (is.null(varX) & errors == "het") {
		stop("You must supply an estimate for the variance of X when the errors are heteroscedastic.")
	}

	# Perform appropriate bandwidth calculation --------------------------------

	if (algorithm == "CV"){
		output <- CVdeconv(n, W, errortype, sigU, phiU, phiK, muK2, RK, deltat, 
						   tt)
	}

	if (algorithm == "PI" & errors == "est") {
		stop("You must define the error distribution")
		# output <- PI_DeconvUEstTh4(W, phiU, hatvarU, phiK, muK2, tt)
	}

	if (algorithm == "PI" & errors == "het") {
		output <- PI_deconvUknownth4het(n, W, varX, errortype, sigU, phiU, phiK,
										muK2, RK, deltat, tt)
	}

	if (algorithm == "PI" & errors == "hom") {
		output <- PI_deconvUknownth4(n, W, errortype, sigU, phiU, phiK, muK2, 
									 RK, deltat, tt)
	}

	output
}