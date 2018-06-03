library(deconvolve)
context("SIMEX")

test_that("replicates SIMEX method 1", {
	skip_on_cran()
	set.seed(1)
	n <- 50
	sd_X <- 1
	sd_U <- 0.2
	data <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm", replicates = TRUE)
	Y <- 2*data$W1
	output_test <- bandwidth(data$W1, data$W2, algorithm = "SIMEX", Y = Y, seed = 100, n_cores = 2)
	expect_equal(output_test$h, 0.07305357, tolerance = 0.0000001)
	expect_equal(output_test$rho, 0.09636439, tolerance = 0.0000001)
})

test_that("replicates SIMEX method 2", {
	skip_on_cran()
	set.seed(1)
	n <- 50
	sd_X <- 1
	sd_U <- 0.2
	data <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm", replicates = TRUE)
	Y <- 2*data$W1
	output_test <- bandwidth(data$W1, data$W2, algorithm = "SIMEX", Y = Y, seed = 100, use_alt_SIMEX_rep_opt = TRUE, n_cores = 2)
	expect_equal(output_test$h, 0.06711634, tolerance = 0.0000001)
	expect_equal(output_test$rho, 0.004956175, tolerance = 0.0000001)
})