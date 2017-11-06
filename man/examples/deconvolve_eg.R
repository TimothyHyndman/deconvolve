# Symmetric Errors -------------------------------------------------------------
n <- 200
W <- GenerateTestData(n, dist_type = "chi", error_type = "norm")
xx <- seq(min(W), max(W), length.out = 100)

d <- deconvolve(W, xx)
plot(d)
print(d)

# Symmetric Errors only returning PMF ------------------------------------------
d <- deconvolve(W, pmf = TRUE)
plot(d)
print(d)

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

yy <- deconvolve(W, xx, errortype = "norm", sigU = sigU_vec)

# Heteroscedastic Errors provided using a vector of phiUs ----------------------
n <- 200
sigX <- 1
sigU <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
W <- GenerateTestData(n, sigX, sigU, dist_type = "mix", error_type = "norm")
xx <- seq(min(W), max(W), length.out = 100)

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

yy <- deconvolve(W, xx, sigU = sigU, phiU = phiU_vec)
