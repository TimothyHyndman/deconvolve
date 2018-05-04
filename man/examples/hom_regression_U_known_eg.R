\dontrun{
# Generate some data to work with
n <- 100
sd_U <- 0.2
X <- stats::rchisq(n, df = 3)
U <- stats::rnorm(n, sd = sd_U)
W <- X + U
Y <- 2 * X

# Create phiU
phiU <- create_phiU(sd_U, "norm")

# Calculate a bandwidth
output <- bandwidth(W, sd_U = sd_U, errortype = "norm", Y = Y, algorithm = "SIMEX")
h <- output$h
rho <- output$rho

# Perform regression
output <- hom_regression_U_known(W, Y, phiU, h, rho)
}
