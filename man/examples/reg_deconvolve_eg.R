n <- 200
sigX <- 1
sigU <- 0.2
data <- GenerateTestData(n, sigX, sigU, dist_type = "mix", error_type = "norm", create_Y = TRUE)
xx <- seq(-2, 8, 0.1)
output <- reg_deconvolve(xx = xx, W = data$W, Y = data$Y, error_type = "norm", sigU = 0.2)
