#' @export
# Just for testing
ComputePhiPmf <- function(theta, p, tt){
	# Returns phi calculated on tt as a complex vector

	phi <- p %*% exp( complex(imaginary = 1) * matrix(theta, ncol = 1) %*% tt )

	return(as.vector(phi))
}

#' @export
# Just for testing
ComputePhiEmp <- function(W, tt){
	# Returns a list containing Phi and the Real, Imaginary, Normed parts of Phi
	n <- length( W )
	
	# Estimate empirical characteristic function of W on tt
	oo <- outerop( tt, W, '*' )
	re.hat.phi <- rowSums( cos( oo ) ) / n
	im.hat.phi <- rowSums( sin( oo ) ) / n
	norm.hat.phi <- sqrt( re.hat.phi^2 + im.hat.phi^2 )

	return(list("complex" = complex(real = re.hat.phi, imaginary = im.hat.phi),
				"re" = re.hat.phi, "im" = im.hat.phi, "norm" = norm.hat.phi, 
				"t.values" = tt))
}