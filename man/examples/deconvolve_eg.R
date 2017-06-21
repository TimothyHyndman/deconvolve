# Generate some data
W <- GenerateTestData()

# Define a vector on which to calculate the density
xx <- seq(min(W), max(W), length.out = 100)

# Calculate the DKDE assuming the errors are symmetric
yy <- deconvolve(W, xx)
