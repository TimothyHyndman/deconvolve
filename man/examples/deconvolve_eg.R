\dontrun{
# Error estimated from replicates ----------------------------------------------
W1 <- (framingham$SBP21 + framingham$SBP22)/2
W2 <- (framingham$SBP31 + framingham$SBP32)/2

yy <- deconvolve(W1, W2)

# Symmetric Errors -------------------------------------------------------------
output <- deconvolve((framingham$SBP21 + framingham$SBP22)/2)

# Generate homoscedastic data --------------------------------------------------
n <- 50
X <- stats::rchisq(n, 3)

sd_U = 0.2
U <- stats::rnorm(n, sd = sd_U)

W <- X + U

# Homoscedastic Errors ---------------------------------------------------------
yy <- deconvolve(W, errortype = "norm", sd_U = sd_U)

# Generate heteroscedastic data ------------------------------------------------
n <- 50
X <- stats::rchisq(n, 3)

sd_U_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
U <- c()
for (sigUk in sd_U_vec){
	U <- c(U, stats::rnorm(1, 0, sigUk))
}

W <- X + U

# Heteroscedastic Errors -------------------------------------------------------
yy <- deconvolve(W, errortype = "norm", sd_U = sd_U_vec)

# Heteroscedastic Errors provided using a vector of phiUs ----------------------
phiU <- c()
for (sigUk in sd_U_vec){
	phiUk <- function(tt) {
		exp(-sigUk^2 * tt^2 / 2)
	}
	phiU <- c(phiU, phiUk)
}

yy <- deconvolve(W, sd_U = sd_U_vec, phiU = phiU)
}
