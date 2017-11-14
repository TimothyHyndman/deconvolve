n <- 200
sd_X <- 1
sd_U <- 0.2
data <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm", create_Y = TRUE)
xx <- seq(-2, 8, 0.1)
output <- reg_deconvolve(data$W, data$Y, xx, errortype = "norm", sd_U = 0.2, n_cores = 2)
