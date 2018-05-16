# Homoscedastic Errors ---------------------------------------------------------
n <- 50
sd_X <- 1
sd_U <- 0.2
W <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm")
yy <- deconvolve(W, errortype = "norm", sd_U = sd_U)

\dontrun{
# Heteroscedastic Errors -------------------------------------------------------
sd_U_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
W <- GenerateTestData(n, sd_X, sd_U_vec, dist_type = "mix", error_type = "norm")

yy <- deconvolve(W, errortype = "norm", sd_U = sd_U_vec)

# Heteroscedastic Errors provided using a vector of phiUs ----------------------
W <- GenerateTestData(n, sd_X, sd_U_vec, dist_type = "mix", error_type = "norm")

phiU_vec=c()
phiU <- function(k) {
	function(tt){
		exp(-sd_U_vec[k]^2 * tt^2 / 2)
	}
}
for(k in 1:n) {	
	phiU_vec <- c(phiU_vec, phiU(k))
}

yy <- deconvolve(W, sd_U = sd_U_vec, phiU = phiU_vec)

# Error estimated from replicates ----------------------------------------------
W1 <- (framingham$SBP21 + framingham$SBP22)/2
W2 <- (framingham$SBP31 + framingham$SBP32)/2

yy <- deconvolve(W1, W2)

# Symmetric Errors -------------------------------------------------------------
output <- deconvolve((framingham$SBP21 + framingham$SBP22)/2)
}

