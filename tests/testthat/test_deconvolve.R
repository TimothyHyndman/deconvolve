library(deconvolve)
context("Deconvolution")

set.seed(1)
n <- 50
W <- GenerateTestData(n, dist_type = "chi", error_type = "norm")

test_that("symmetric error case gives expected result", {
	d_test <- deconvolve(W)
	load("sym_error_test_result.RData")
	expect_equal(d_test, d)

	set.seed(1)
	d_test <- deconvolve(W, pmf = TRUE)
	load("sym_error_pmf_test_result.RData")
	expect_equal(d_test, d)
})

set.seed(1)
sd_X <- 1
sd_U <- 0.2
W <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm")

test_that("hom error gives expected result", {
	yy_test <- deconvolve(W, errortype = "norm", sd_U = sd_U)
	load("hom_error_test_result.RData")
	expect_equal(yy_test, yy)
})

set.seed(1)
sd_U_vec <- 0.6 * sqrt(1 + (1:n) / n) * sqrt(0.5)
W <- GenerateTestData(n, sd_X, sd_U_vec, dist_type = "mix", error_type = "norm")
phiU_vec=c()
phiU <- function(k) {
	function(tt){
		exp(-sd_U_vec[k]^2 * tt^2 / 2)
	}
}
for(k in 1:n) {	
	phiU_vec <- c(phiU_vec, phiU(k))
}


test_that("het error gives expected result", {
	yy_test <- deconvolve(W, sd_U = sd_U_vec, phiU = phiU_vec)
	load("het_error_test_result.RData")
	expect_equal(yy_test, yy)
})