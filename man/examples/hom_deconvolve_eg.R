\dontrun{
# Symmetric Errors -------------------------------------------------------------
output <- hom_deconvolve(framingham$SBP21)

# Error estimated from replicates ----------------------------------------------
W1 <- framingham$SBP21
W2 <- framingham$SBP22

yy <- hom_deconvolve(W1, W2)
}
