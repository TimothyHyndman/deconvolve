#' @export
#' 
create_phiU <- function(sd_U, 
						error_type = c("normal", "laplace")) {

	error_type <- match.arg(error_type)

	if (length(sd_U) == 1){
		errors = "hom"
	} else {
		errors = "het"
	}

	if(error_type == 'laplace' & errors == "hom") {
		phiU <- function(tt) {
			1 / (1 + sd_U^2 * tt^2 / 2)
		}
	}

	if(error_type == 'normal' & errors == "hom") {
		phiU <- function(tt) {
			exp(-sd_U^2 * tt^2 / 2)
		}
	}

	if(error_type == 'laplace' & errors == "het") {
		phiU <- c()
		for (sigUk in sd_U){
			phiUk <- function(tt) {
				1 / (1 + sigUk^2 * tt^2 / 2)
			}
			phiU <- c(phiU, phiUk)
		}
	}

	if(error_type == 'normal' & errors == "het") {
		phiU <- c()
		for (sigUk in sd_U){
			phiUk <- function(tt) {
				exp(-sigUk^2 * tt^2 / 2)
			}
			phiU <- c(phiU, phiUk)
		}
	}

	phiU
}