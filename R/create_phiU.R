create_phiU <- function(errors, errortype, sd_U){

	if(errortype == 'Lap' & errors == "hom") {
		phiU <- function(tt) {
			1 / (1 + sd_U^2 * tt^2 / 2)
		}
	}

	if(errortype == 'norm' & errors == "hom") {
		phiU <- function(tt) {
			exp(-sd_U^2 * tt^2 / 2)
		}
	}

	if(errortype == 'Lap' & errors == "het") {
		phiU <- c()
		for (sigUk in sd_U){
			phiUk <- function(tt) {
				1 / (1 + sigUk^2 * tt^2 / 2)
			}
			phiU <- c(phiU, phiUk)
		}
	}

	if(errortype == 'norm' & errors == "het") {
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