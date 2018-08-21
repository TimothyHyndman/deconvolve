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

	# NlcOptim works different on different OS (I think) so we need multiple tests
	if (Sys.info()['sysname'] == "Windows"){
		expect_equal(bandwidth(W), 0.1104239, tolerance = 0.0000001)
	} else if (Sys.info()['sysname'] == "Darwin"){
		expect_equal(bandwidth(W), 0.1104239, tolerance = 0.0000001)
	} else if (Sys.info()['sysname'] == "Linux"){
		expect_equal(bandwidth(W), 0.1104239, tolerance = 0.0000001)
	} else {
		warning("OS wasn't one of Windows, Darwin, or Linux and so test didn't work.")
	}
	
})

set.seed(1)
n <- 50
W <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm")
test_that("CV case gives expected result", {
	skip_on_cran()
	expect_equal(bandwidth(W, errortype = "norm", sd_U = sd_U, algorithm = "CV"), 0.1387812, tolerance = 0.0000001)
})


