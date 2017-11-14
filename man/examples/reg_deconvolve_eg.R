n <- 200
sigX <- 1
sigU <- 0.2
data <- GenerateTestData(n, sigX, sigU, dist_type = "mix", error_type = "norm", create_Y = TRUE)
xx <- seq(-2, 8, 0.1)
output <- reg_deconvolve(data$W, data$Y, xx, errortype = "norm", sigU = 0.2, n_cores = 2)
