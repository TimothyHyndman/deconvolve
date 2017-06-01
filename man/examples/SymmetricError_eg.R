# Generate Test Data
W <- GenerateTestData(n = 500)

# Deconvolve
out <- DeconErrSymPmf(W)

phi.W <- out$phi.W

# Convert to Probability Density
fX <- DeconErrSymPmfToPdf(out, W, phi.W)

# Plot results
PlotPmfAndPdf(out, fX)
