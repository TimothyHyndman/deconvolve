\dontrun{
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

# Error estimated from replicates
data <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm", replicates = TRUE)
bw <- bandwidth(data$W1, data$W2)
}
