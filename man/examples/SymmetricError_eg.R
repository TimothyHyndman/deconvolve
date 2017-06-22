# Generate Test Data
W <- GenerateTestData(n = 500)

# Deconvolve
out <- DeconErrSymPmf(W)

phi.W <- out$phi.W

# Convert to Probability Density
xx <- seq(min(W), max(W), length.out = 100)
fX <- DeconErrSymPmfToPdf(out, W, phi.W, xx)
