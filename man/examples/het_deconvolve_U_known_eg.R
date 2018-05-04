# Generate some data to work with
n <- 100
sd_U <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
X <- stats::rchisq(n, df = 3)
for (sd in sd_U) {
	U[i] <- stats::rnorm(1, 0, sd)
}
W <- X + U

# Create phiUs
phiU <- create_phiU(sd_U, "norm")

# Calculate a bandwidth
h <- bandwidth(W, sd_U = sd_U, phiU = phiU)

# Deconvolve
yy <- deconvolve(W, phiU, h)
