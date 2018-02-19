library(deconvolve)
context("reg_deconvolve")

set.seed(1)
n <- 50
sd_X <- 1
sd_U <- 0.2
data <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm", create_Y = TRUE)


test_that("reg_deconvolve gives expected result", {
	output_test <- reg_deconvolve(data$W, data$Y, errortype = "norm", sd_U = 0.2, n_cores = 2)
	load("reg_deconvolve_test_result.RData")
	expect_equal(output_test, output)
})