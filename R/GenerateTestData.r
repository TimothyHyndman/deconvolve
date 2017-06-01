#' @export

GenerateTestData <- function(n = 500, NSR = 0.2){
	
	# True 
	df <- 3
	X <- stats::rchisq(n, df)
	var.X <- 2*df
	X <- X / sqrt(var.X)
	var.X <- 1

	# var.X.long <- 2 * df
	# c.var <- sqrt(var.X.long)
	# true.dens <- function(xx){
	# 	c.var * dchisq(xx * cvar, df)
	# }

	# Error
	sig.U = sqrt(NSR * var.X)
	U <- stats::rnorm(n, sd = sig.U)

	W <- X + U
	return(W)
}