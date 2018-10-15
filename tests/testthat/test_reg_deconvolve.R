library(deconvolve)
context("reg_deconvolve")
Sys.unsetenv("R_TESTS")

set.seed(1)
n <- 50
sd_X <- 1
sd_U <- 0.2
data <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm", create_Y = TRUE)


test_that("reg_deconvolve gives expected result", {
	skip_on_cran()
	output_test <- reg_deconvolve(data$Y, data$W, errortype = "norm", sd_U = 0.2, bw = 0.06, rho = 0.09)
	load("reg_deconvolve_test_result.RData")
	expect_equal(output_test, output)
})

test_that("replicates reg_deconvolve", {
	skip_on_cran()
	set.seed(1)
	data <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm", replicates = TRUE)
	Y <- 2*data$W1
	output_test <- reg_deconvolve(Y, data$W1, W2 = data$W2, bw = 0.06, rho = 0.09)
	load("reg_deconvolve_rep_test_result.RData")
	expect_equal(output_test, output)
})