\dontrun{
# Error from replicates --------------------------------------------------------
W1 <- (framingham$SBP21 + framingham$SBP22)/2
W2 <- (framingham$SBP31 + framingham$SBP32)/2
Y <- framingham$FIRSTCHD
h <- 1.120537 #Precalculated using SIMEX option from bandwidth()
rho <- 0.0103959 #Precalculated using SIMEX option from bandwidth()
output <- reg_deconvolve(Y, W1, W2, bw = h, rho = rho)

# Error known ------------------------------------------------------------------
n <- 50
X <- stats::rchisq(n, 3)
Y <- 2*X

sd_U = 0.2
U <- stats::rnorm(n, sd = sd_U)

W <- X + U

output <- reg_deconvolve(W, Y, errortype = "norm", sd_U = 0.2, n_cores = 2)
}
