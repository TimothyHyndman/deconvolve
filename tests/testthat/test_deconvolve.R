library(deconvolve)
context("Deconvolution")

set.seed(1)
n <- 50
W <- GenerateTestData(n, dist_type = "chi", error_type = "norm")

# test_that("symmetric error case gives expected result", {
# 	skip_on_cran()
# 	d_test <- hom_deconvolve(W)
# 	load("sym_error_test_result.RData")
# 	expect_equal(d_test, d)

# 	set.seed(1)
# 	d_test <- hom_deconvolve(W, pmf = TRUE)
# 	load("sym_error_pmf_test_result.RData")
# 	expect_equal(d_test, d)
# })

test_that("replicates gives expected result", {
	skip_on_cran()
	set.seed(1)
	sd_X <- 1
	sd_U <- 0.2
	n <- 50
	data <- GenerateTestData(n, sd_X, sd_U, dist_type = "chi", error_type = "norm", 
		replicates = TRUE)
	yy_test <- hom_deconvolve(data$W1, data$W2)
	load("hom_error_rep_test_result.RData")
	expect_equal(yy_test, yy)
})

set.seed(1)
sd_X <- 1
sd_U <- 0.2
W <- GenerateTestData(n, sd_X, sd_U, dist_type = "mix", error_type = "norm")

test_that("hom error gives expected result", {
	skip_on_cran()
	phiU <- create_phiU(sd_U, "norm")
	h <- bandwidth(W, NULL, "norm", sd_U, phiU)
	yy_test <- hom_deconvolve_U_known(W, phiU, h)
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
	h <- bandwidth(W, NULL, NULL, sd_U_vec, phiU_vec)
	yy_test <- het_deconvolve_U_known(W, phiU_vec, h)
	load("het_error_test_result.RData")
	expect_equal(yy_test, yy)
})