# Just used for convenience in unit testing

GenerateTestData <- function(n, sd_X = 1, sd_U = 0.2, dist_type = "chi", 
							 error_type = "norm", create_Y = FALSE, 
							 replicates = FALSE){
	
	error_types <- c("normal", "laplace")
    error_type <- error_types[pmatch(tolower(error_type), error_types)]
    if (is.na(error_type)) {
            stop("Please provide a valid errortype.")
    }


	# Sample true data ---------------------------------------------------------
	if (dist_type == "chi"){
		df <- 3
		var_X <- 2*df
		X <- stats::rchisq(n, df)
		scaling_var <- (sqrt(var_X) / sd_X)
		X <- X / scaling_var
		var_X <- sd_X^2
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
		scaling_var <- (sqrt(var_X) / sd_X)
		X <- X / scaling_var
		var_X <- sd_X^2
		true_dens <- function(xx){
			xx <- xx * scaling_var
			(1 - pmix) * stats::dnorm(xx, mu1, sig1) + 
			pmix * stats::dnorm(xx, mu2, sig2)
		}
	}


	# Sample error -------------------------------------------------------------
	if (error_type == "normal"){
		if (length(sd_U) == 1){
			U <- stats::rnorm(n, sd = sd_U)
			U2 <- stats::rnorm(n, sd = sd_U)
		} else {
			U <- numeric(n)
			U2 <- numeric(n)
			for (i in 1:n) {
				U[i] <- stats::rnorm(1, 0, sd_U[i])
				U2[i] <- stats::rnorm(1, 0, sd_U[i])
			}
		}
	}

	if (error_type == "laplace") {
		if (length(sd_U) == 1) {
			sigLap <- sd_U / sqrt(2)
			U <- rlap(sigLap, 1, n)
			U2 <- rlap(sigLap, 1, n)
		} else {
			sigLap <- sd_U / sqrt(2)
			U <- numeric(n)
			U2 <- numeric(n)
			for (i in 1:n) {
				U[i] <- rlap(sigLap[i], 1, 1)
				U2[i] <- rlap(sigLap[i], 1, 1)
			}
		}
	}

	# Combine samples ----------------------------------------------------------
	W <- X + U
	W2 <- X + U2

	if (create_Y) {
		Y <- 2*X
		output <- list("W" = W, "Y" = Y)
	} else if (replicates) {
		output <- list("W1" = W, "W2" = W2)
	} else {
		output <- W
	}
	
	output
}
