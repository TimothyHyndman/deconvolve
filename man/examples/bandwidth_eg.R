\dontrun{
# PI bandwidth with error estimated from replicates ----------------------------
W1 <- (framingham$SBP21 + framingham$SBP22)/2
W2 <- (framingham$SBP31 + framingham$SBP32)/2

bw <- bandwidth(W1, W2)

# CV bandwidth -----------------------------------------------------------------
n <- 50
sd_X <- 1
sd_U <- 0.2
W <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm")

bw <- bandwidth(W, errortype = "norm", sd_U = sd_U, algorithm = "CV")

# PI bandwidth with heteroscedastic errors -------------------------------------
sd_U_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
W <- GenerateTestData(n, sd_X, sd_U_vec, dist_type = "mix", error_type = "norm")

bw <- bandwidth(W, errortype = "norm", sd_U = sd_U_vec)

# SIMEX bandwidth --------------------------------------------------------------
data <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm", 
						 create_Y = TRUE)
output <- bandwidth(data$W, errortype = "norm", sd_U = sd_U, Y = data$Y, 
					algorithm = "SIMEX", n_cores = 2)
bw <- output$h

# PI bandwidth with heteroscedastic errors supplied using phiU -----------------
sd_U_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
phiU <- c()
for (sigUk in sd_U_vec){
	phiUk <- function(tt) {
		exp(-sigUk^2 * tt^2 / 2)
	}
	phiU <- c(phiU, phiUk)
}

W <- GenerateTestData(n, sd_X, sd_U_vec, dist_type = "mix", error_type = "norm")

bw <- bandwidth(W, sd_U = sd_U_vec, phiU = phiU)
}
