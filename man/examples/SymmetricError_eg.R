# Generate Test Data

n <- 500
W <- GenerateTestData(n)

# Deconvolve
out <- DeconvolveSymmetricErrorPmf(W, m = 10)
theta <- out$support
p <- out$probweights
phi.W <- out$phi.W

# Convert to Probability Density
fX <- SymmetricErrorPmfToPdf(theta, p, W, phi.W)

# Plot results
par(mfcol = c(1,2))

plot(theta, p)
plot(fX$x, fX$y,"l")
