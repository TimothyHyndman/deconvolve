\dontrun{
# PI bandwidth with error estimated from replicates ----------------------------
W1 <- (framingham$SBP21 + framingham$SBP22)/2
W2 <- (framingham$SBP31 + framingham$SBP32)/2

bw <- bandwidth(W1, W2)


# Generate homoscedastic data --------------------------------------------------
n <- 50
X <- stats::rchisq(n, 3)

sd_U = 0.2
U <- stats::rnorm(n, sd = sd_U)

W <- X + U

# CV bandwidth -----------------------------------------------------------------
bw <- bandwidth(W, errortype = "norm", sd_U = sd_U, algorithm = "CV")

# SIMEX bandwidth --------------------------------------------------------------
Y <- 2*X

output <- bandwidth(W, errortype = "norm", sd_U = sd_U, Y = Y, 
					algorithm = "SIMEX", n_cores = 2)
bw <- output$h
rho <- output$rho

# Generate heteroscedastic data ------------------------------------------------
n <- 50
X <- stats::rchisq(n, 3)

sd_U_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
U <- c()
for (sigUk in sd_U_vec){
	U <- c(U, stats::rnorm(1, 0, sigUk))
}

W <- X + U

# PI bandwidth with heteroscedastic errors -------------------------------------
bw <- bandwidth(W, errortype = "norm", sd_U = sd_U_vec)

# PI bandwidth with heteroscedastic errors supplied using phiU -----------------
phiU <- c()
for (sigUk in sd_U_vec){
	phiUk <- function(tt) {
		exp(-sigUk^2 * tt^2 / 2)
	}
	phiU <- c(phiU, phiUk)
}

bw <- bandwidth(W, sd_U = sd_U_vec, phiU = phiU)
}
