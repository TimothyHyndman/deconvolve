# Generate some data to work with
n <- 100
sd_U <- 0.2
X <- stats::rchisq(n, df = 3)
U <- stats::rnorm(n, sd = sd_U)
W <- X + U

# Create phiU
phiU <- create_phiU(sd_U, "norm")

# Calculate a bandwidth
h <- bandwidth(W, sd_U = sd_U, phiU = phiU)

# Perform deconvolution
y <- hom_deconvolve_U_known(W, phiU, h)
