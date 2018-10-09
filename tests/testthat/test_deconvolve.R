library(deconvolve)
context("Deconvolution")

set.seed(1)
n <- 50
W <- GenerateTestData(n, dist_type = "chi", error_type = "norm")

test_that("symmetric error pmftopdf gives expected result", {
	skip_on_cran()
	load('sym_error_pmf_input.RData')
	d_test <- DeconErrSymPmfToPdf(X_pmf, W, phi_W, seq(min(W), max(W), length.out = 100), 'default', 1, NULL)
	load("sym_error_pmftopdf_test_result.RData")
	expect_equal(d_test, d)
})

test_that("symmetric error case gives reasonable result", {
	d_test <- deconvolve(W, m = 2)
	expect_equal(length(d_test$pdf), 100)
})

test_that("replicates gives expected result", {
  skip_on_cran()
  set.seed(1)
  sd_X <- 1
  sd_U <- 0.2
  n <- 50
  data <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm", 
    replicates = TRUE)
  yy_test <- deconvolve(data$W1, data$W2)
  load("hom_error_rep_test_result.RData")
  expect_equal(yy_test, yy)
})

test_that("het_replicates", {
	skip_on_cran()
	set.seed(1)
	sd_X <- 1
	sd_U <- 0.2
	n <- 50
	data <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm", 
	replicates = TRUE)
	yy_test <- deconvolve(data$W1, data$W2, het_replicates = TRUE)
	load("het_error_rep_test_result.RData")
	expect_equal(yy_test, yy)
})


set.seed(1)
sd_X <- 1
sd_U <- 0.2
W <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm")

test_that("hom error gives expected result", {
	skip_on_cran()
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
	skip_on_cran()
	yy_test <- deconvolve(W, sd_U = sd_U_vec, phiU = phiU_vec)
	load("het_error_test_result.RData")
	expect_equal(yy_test, yy)
})