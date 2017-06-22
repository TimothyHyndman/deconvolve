# Symmetric Errors -------------------------------------------------------------
n <- 200
W <- GenerateTestData(n, dist_type = "chi", error_type = "norm")
xx <- seq(min(W), max(W), length.out = 100)

yy <- deconvolve(W, xx)

# Homoscedastic Errors ---------------------------------------------------------
n <- 200
sigX <- 1
sigU <- 0.2
W <- GenerateTestData(n, sigX, sigU, dist_type = "mix", error_type = "norm")
xx <- seq(min(W), max(W), length.out = 100)

yy <- deconvolve(W, xx, errortype = "norm", sigU = sigU)

# Heteroscedastic Errors -------------------------------------------------------
n <- 200
sigX <- 1
sigU_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
W <- GenerateTestData(n, sigX, sigU_vec, dist_type = "mix", error_type = "norm")
xx <- seq(min(W), max(W), length.out = 100)
# Estimate the variance of X
varX <- mean(W^2) - (mean(W))^2 - sum(sigU_vec^2) / n

yy <- deconvolve(W, xx, errortype = "norm", sigU = sigU_vec, varX = varX)

# Heteroscedastic Errors provided using a vector of phiUs ----------------------
n <- 200
sigX <- 1
sigU <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
W <- GenerateTestData(n, sigX, sigU, dist_type = "mix", error_type = "norm")
xx <- seq(min(W), max(W), length.out = 100)
# Estimate the variance of X
varX <- mean(W^2) - (mean(W))^2 - sum(sigU^2) / n
# Calculate the characteristic functions of the errors
phiU_vec=c()
phiU <- function(k) {
	function(tt){
		exp(-sigU[k]^2 * tt^2 / 2)
	}
}
for(k in 1:n) {	
	# phiU <- function(tt,k) {
	# 	exp(-sigU[k]^2 * tt^2 / 2)
	# }
	phiU_vec <- c(phiU_vec, phiU(k))
}

yy <- deconvolve(W, xx, phiU = phiU_vec, varX = varX)
