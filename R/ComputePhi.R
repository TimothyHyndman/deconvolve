#' @export
ComputePhiPmf <- function(theta, p, tt){
	# Returns phi calculated on tt as a complex vector
	i <- complex(imaginary = 1)
	phi <- p %*% exp( i * matrix(theta, ncol = 1) %*% tt )

	return(as.vector(phi))
}
#' @export
ComputePhiEmp <- function(W, tt){
	# Returns a list containing Phi and the Real, Imaginary, Normed parts of Phi
	n <- length( W )
	
	i <- complex(imaginary = 1)
	phi <- (1 / n) * 
		   colSums( exp( i * matrix(W, ncol = 1) %*% matrix(tt, nrow = 1) ) )

	return(list("complex" = phi, "re" = Re(phi), "im" = Im(phi), 
				"norm" = Mod(phi), "t.values" = tt))
}