library(deconvolve)
context("Bandwidth")

# load("sym_error_test_result.RData")

set.seed(1)
n <- 50
sd_X <- 1
sd_U <- 0.2
W <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm")

set.seed(1)
test_that("homoscedastic errors case gives expected result", {
	skip_on_cran()
	expect_equal(bandwidth(W, errortype = "norm", sd_U = sd_U), 0.1278101, tolerance = 0.0000001)
})

set.seed(1)
sd_U_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
W <- GenerateTestData(n, sd_X, sd_U_vec, dist_type = "mix", error_type = "norm")
test_that("heteroscedastic error case gives expected result", {
	skip_on_cran()
	expect_equal(bandwidth(W, errortype = "norm", sd_U = sd_U_vec), 0.2312728, tolerance = 0.0000001)
})

set.seed(1)
data <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm", replicates = TRUE)
test_that("replicates case gives expected result", {
	skip_on_cran()
	expect_equal(bandwidth(data$W1, data$W2), 0.09959522, tolerance = 0.0000001)
})

set.seed(2)
n <- 500
test_that("no error case gives expected result", {
	skip_on_cran()
	W <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm")
	expect_equal(bandwidth(W), 0.1104239, tolerance = 0.0000001)
})

set.seed(1)
n <- 50
W <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm")
test_that("CV case gives expected result", {
	skip_on_cran()
	expect_equal(bandwidth(W, errortype = "norm", sd_U = sd_U, algorithm = "CV"), 0.1387812, tolerance = 0.0000001)
})

test_that("replicates SIMEX method 1", {
	set.seed(1)
	n <- 50
	data <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm", replicates = TRUE)
	Y <- 2*data$W1
	output_test <- bandwidth(data$W1, data$W2, algorithm = "SIMEX", Y = Y, seed = 100, n_cores = 2)
	expect_equal(output_test$h, 0.07305357, tolerance = 0.0000001)
	expect_equal(output_test$rho, 0.09636439, tolerance = 0.0000001)
})

test_that("replicates SIMEX method 2", {
	set.seed(1)
	n <- 50
	data <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm", replicates = TRUE)
	Y <- 2*data$W1
	output_test <- bandwidth(data$W1, data$W2, algorithm = "SIMEX", Y = Y, seed = 100, use_alt_SIMEX_rep_opt = TRUE, n_cores = 2)
	expect_equal(output_test$h, 0.06711634, tolerance = 0.0000001)
	expect_equal(output_test$rho, 0.004956175, tolerance = 0.0000001)
})
