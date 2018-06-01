\dontrun{
# Error from replicates --------------------------------------------------------
W1 <- (framingham$SBP21 + framingham$SBP22)/2
W2 <- (framingham$SBP31 + framingham$SBP32)/2
Y <- framingham$FIRSTCHD
h <- 2
rho <- 1
output <- reg_deconvolve(Y, W1, W2, bw = h, rho = rho)

# Error known ------------------------------------------------------------------
n <- 50
sd_X <- 1
sd_U <- 0.2
data <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm", create_Y = TRUE)
output <- reg_deconvolve(data$W, data$Y, errortype = "norm", sd_U = 0.2, n_cores = 2)
}
