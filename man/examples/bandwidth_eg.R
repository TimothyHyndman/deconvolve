# CV bandwidth -----------------------------------------------------------------
n <- 200
sigX <- 1
sigU <- 0.2
W <- GenerateTestData(n, sigX, sigU, dist_type = "mix", error_type = "norm")

bw <- bandwidth(W, errortype = "norm", sigU = sigU, algorithm = "CV")


# PI bandwidth with heteroscedastic errors -------------------------------------
n <- 200
sigX <- 1
sigU_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
W <- GenerateTestData(n, sigX, sigU_vec, dist_type = "mix", error_type = "norm")
# Estimate the variance of X
varX <- mean(W^2) - (mean(W))^2 - sum(sigU_vec^2) / n

bw <- bandwidth(W, errortype = "norm", sigU = sigU_vec, varX = varX)

# SIMEX bandwidth --------------------------------------------------------------
n <- 200
sigX <- 1
sigU <- 0.2
data <- GenerateTestData(n, sigX, sigU, dist_type = "mix", error_type = "norm", 
						 create_Y = TRUE)
output <- bandwidth(data$W, errortype = "norm", sigU = sigU, Y = data$Y, 
					algorithm = "SIMEX", n_cores = 2)
bw <- output$h

# PI bandwidth with heteroscedastic errors supplied using phiU -----------------
n <- 200
sigX <- 1
sigU_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
phiU <- c()
for (sigUk in sigU_vec){
	phiUk <- function(tt) {
		exp(-sigUk^2 * tt^2 / 2)
	}
	phiU <- c(phiU, phiUk)
}

W <- GenerateTestData(n, sigX, sigU_vec, dist_type = "mix", error_type = "norm")

# Estimate the variance of X
varX <- mean(W^2) - (mean(W))^2 - sum(sigU_vec^2) / n

bw <- bandwidth(W, phiU = phiU, varX = varX)
