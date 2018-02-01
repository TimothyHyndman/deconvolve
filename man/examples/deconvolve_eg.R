# Symmetric Errors -------------------------------------------------------------
n <- 200
W <- GenerateTestData(n, dist_type = "chi", error_type = "norm")
d <- deconvolve(W)
plot(d)
print(d)

# Symmetric Errors only returning PMF ------------------------------------------
d <- deconvolve(W, pmf = TRUE)
plot(d)
print(d)

# Homoscedastic Errors ---------------------------------------------------------
n <- 200
sd_X <- 1
sd_U <- 0.2
W <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm")

yy <- deconvolve(W, errortype = "norm", sd_U = sd_U)

# Heteroscedastic Errors -------------------------------------------------------
n <- 200
sd_X <- 1
sd_U_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
W <- GenerateTestData(n, sd_X, sd_U_vec, dist_type = "mix", error_type = "norm")

yy <- deconvolve(W, errortype = "norm", sd_U = sd_U_vec)

# Heteroscedastic Errors provided using a vector of phiUs ----------------------
n <- 200
sd_X <- 1
sd_U <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
W <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm")

phiU_vec=c()
phiU <- function(k) {
	function(tt){
		exp(-sd_U[k]^2 * tt^2 / 2)
	}
}
for(k in 1:n) {	
	phiU_vec <- c(phiU_vec, phiU(k))
}

yy <- deconvolve(W, sd_U = sd_U, phiU = phiU_vec)
