#' @export

GenerateTestData <- function(n = 500, sig_X = 1, sig_U = 0.2, dist_type = "chi", 
							 error_type = "norm"){
	
	# Sample true data ---------------------------------------------------------
	if (dist_type == "chi"){
		df <- 3
		var_X <- 2*df
		X <- stats::rchisq(n, df)
		scaling_var <- (sqrt(var_X) / sig_X)
		X <- X / scaling_var
		var_X <- sig_X^2
		true_dens <- function(xx){
			scaling_var * stats::dchisq(xx * scaling_var, df)
		}
	}

	if (dist_type == "mix"){
		mu1 <- -3
		mu2 <- 2
		sig1 <- 1
		sig2 <- 1
		X <- stats::rnorm(n, mu1, sig1)
		X2 <- stats::rnorm(n, mu2, sig2)
		pmix <- 0.75
		tmp <- matrix(stats::runif(n, 0, 1), nrow=1, ncol=n, byrow=TRUE)
		X[which(tmp < pmix)] <- X2[which(tmp < pmix)]
		var_X <- stats::var(X)
		scaling_var <- (sqrt(var_X) / sig_X)
		X <- X / scaling_var
		var_X <- sig_X^2
		true_dens <- function(xx){
			xx <- xx * scaling_var
			(1 - pmix) * stats::dnorm(xx, mu1, sig1) + 
			pmix * stats::dnorm(xx, mu2, sig2)
		}
	}


	# Sample error -------------------------------------------------------------
	if (error_type == "norm"){
		if (length(sig_U) == 1){
			U <- stats::rnorm(n, sd = sig_U)
		} else {
			U <- numeric(n)
			for (i in 1:n) {
				U[i] <- stats::rnorm(1, 0, sig_U[i])
			}
		}
	}

	# Combine samples ----------------------------------------------------------
	W <- X + U

	W
}