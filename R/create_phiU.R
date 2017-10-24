create_phiU <- function(errors, errortype, sigU){

	if(errortype == 'Lap' & errors == "hom") {
		phiU <- function(tt) {
			1 / (1 + sigU^2 * tt^2 / 2)
		}
	}

	if(errortype == 'norm' & errors == "hom") {
		phiU <- function(tt) {
			exp(-sigU^2 * tt^2 / 2)
		}
	}

	if(errortype == 'Lap' & errors == "het") {
		phiU <- c()
		for (sigUk in sigU){
			phiUk <- function(tt) {
				1 / (1 + sigUk^2 * tt^2 / 2)
			}
			phiU <- c(phiU, phiUk)
		}
	}

	if(errortype == 'norm' & errors == "het") {
		phiU <- c()
		for (sigUk in sigU){
			phiUk <- function(tt) {
				exp(-sigUk^2 * tt^2 / 2)
			}
			phiU <- c(phiU, phiUk)
		}
	}

	phiU
}