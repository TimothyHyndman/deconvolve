DeconErrSymPmf <- function(W, m, kernel_type, n_tp_iter = 20, n_var_iter = 2, 
						   show_diagnostics = FALSE){

	diagnostic <- function(message){
		print_diagnostic(message, show_diagnostics)
	}

	kernel_list <- kernel(kernel_type)
	phiK <- kernel_list$phik
	mu_K2 <- kernel_list$muk2
	RK <- kernel_list$rk

	if (m < 2){
		stop("m must be at least 2")
	}
	
	n <- length(W)

	#--------------------------------------------------------------------------#
	# Pre-calculate PhiW and weight w(t)
	#--------------------------------------------------------------------------#
	# Calculate phi_W on [-8,8] so we can find t*
	tt_length <- 100
	hnaive <- ((8 * sqrt(pi) * RK/3/mu_K2^2)^0.2) * sqrt(stats::var(W)) * 
		n^(-1/5)
	hmin <- hnaive/3
	tt <- seq(-1/hmin, 1/hmin, length.out = tt_length)
	phi_W <- ComputePhiEmp(W, tt)

	# Calculate t*
	t_star <- find_t_cutoff(phi_W$norm, tt, n)

	# Calculate phi_W on [-t*,t*]
	tt_new_length <- 100
	tt_new <- seq(-t_star, t_star, length.out = tt_new_length)
	phi_W <- ComputePhiEmp(W, tt_new)
	sqrt_psi_W <- compute_sqrt_psi_W(tt_new, W)

	# Calculate weight w(t) on [-t*, t*]
	weight <- KernelWeight(tt_new)

	#--------------------------------------------------------------------------#
	# Solve optimization problem to find PMF
	#--------------------------------------------------------------------------#

	matrices <- create_bound_matrices(W, m)
	A <- matrices$A
	B <- matrices$B

	# ----------------------
	# Min T(p)
	# ----------------------

	tp_objective_NLC <- function(x) {
		tp_objective(x, phi_W, sqrt_psi_W, weight)
	}
	diagnostic("Minimizing T(p)")

	tp_min <- Inf
	for (i in 1:n_tp_iter){
		looping <- TRUE
		while (looping) {
			# NlcOptim can't handle it if the Generated QP problem is infeasible
			# However, this error only pops up some of the time. For now, a 
			# dirty solution is to keep running the code with different starting
			# values and ignore everything until we get no error
			looping <- FALSE

			theta0 <- sort(stats::runif(m, min = min(W), max = max(W)))
			p0 <- stats::runif(m, min = 0, max = 1)
			p0 <- p0 / sum(p0)

			x0 <- pmf_to_x(theta0, p0)

			# test_min_tp_sol <- stats::constrOptim(x0, tp_objective, NULL, A_tp, B_tp, 
			# 							 		  phi_W = phi_W, weight = weight, outer.iterations = 10000)
			test_min_tp_sol <- tryCatch({
				test_min_tp_sol <- NlcOptim::solnl(x0, tp_objective_NLC, NULL, A, 
													B, 
													maxIter = 400, 
													maxnFun = 100*(2*m-1), 
													tolX = 1e-6)
			}, error = function(e){
					error_message <- paste("ERROR :", conditionMessage(e))
					diagnostic(error_message)
					test_min_tp_sol <- NULL						
			})

			if (is.null(test_min_tp_sol)) {
					looping <- TRUE
			}
		}

		# if (test_min_tp_sol$value < tp_min){
		# 	tp_min <- test_min_tp_sol$value
		# 	min_tp_sol <- test_min_tp_sol
		# 	diagnostic(tp_min)
		# }
		if (test_min_tp_sol$fn < tp_min){
			tp_min <- test_min_tp_sol$fn
			min_tp_sol <- test_min_tp_sol
			diagnostic(tp_min)
		}
	}

	# con_fun(x) <= 0
	# diagnostic(tp_min)
	# diagnostic(min_tp_sol)

	X_pmf <- x_to_pmf(min_tp_sol$par)
	tt <- phi_W$t.values
	# Calculate phi_X
	phi_X <- ComputePhiPmf(X_pmf$support, X_pmf$prob_weights, tt)

	tp_max <- calculate_tp(phi_X, phi_W, sqrt_psi_W, weight)
	diagnostic(paste("tp_max = ", tp_max))

	penalty_tolerance_scale = 0.00
	penalties_max <- calculate_penalties(phi_X, phi_W) * (1 + penalty_tolerance_scale)
	diagnostic(paste("penalties = ", penalties_max))

	con_fun <- function(x){
		# constraints(x, phi_W, weight, min_tp_sol$value)
		constraints(x, phi_W, sqrt_psi_W, weight, tp_max, penalties_max)
	}
	
	diagnostic(is_feasible(X_pmf$support, X_pmf$prob_weights, A, B, con_fun))

	diagnostic("Minimizing Variance")

	var_min <- Inf
	for (i in 1:n_var_iter){
		looping <- TRUE
		while (looping){			
			# NlcOptim can't handle it if the Generated QP problem is infeasible
			# However, this error only pops up some of the time. For now, a 
			# dirty solution is to keep running the code with different starting
			# values and ignore everything until we get no error
			looping <- FALSE

			theta0 <- sort(stats::runif(m, min = min(W), max = max(W)))
			p0 <- stats::runif(m, min = 0, max = 1)
			p0 <- p0 / sum(p0)

			x0 <- pmf_to_x(theta0, p0)

			test_min_var_sol <- tryCatch({
	    		test_min_var_sol <- NlcOptim::solnl(x0, var_objective, con_fun, 
	    											A, B, 
	    											maxIter = 400, 
													maxnFun = 100*(2*m-1), 
													tolX = 1e-6)
			}, error = function(e){
				error_message <- paste("ERROR :", conditionMessage(e))
				diagnostic(error_message)
				test_min_var_sol <- NULL						
			})

			if (is.null(test_min_var_sol)) {
				looping <- TRUE
			} else {
				X_pmf <- x_to_pmf(test_min_var_sol$par)
				if (is_feasible(X_pmf$support, X_pmf$prob_weights, A, B, con_fun)) {
					diagnostic("feasible solution returned")
				} else {
					diagnostic("non feasible solution returned")
					# looping <- TRUE
				}
			}
		}


		
		if (test_min_var_sol$fn < var_min){
			var_min <- test_min_var_sol$fn
			min_var_sol <- test_min_var_sol
			diagnostic(var_min)
		}
	}

	diagnostic(paste("Initial Variance was", var_objective(min_tp_sol$par)))
	diagnostic(paste("Final Variance is", var_min))
	# tp_diff <- tp_min - tp_objective(min_var_sol$par, phi_W, weight)
	# diagnostic(paste("T(p) decreased by", tp_diff, "while minimizing variance"))

	X_pmf <- x_to_pmf(min_var_sol$par)
	# Calculate phi_X
	tt <- phi_W$t.values
	phi_X <- ComputePhiPmf(X_pmf$support, X_pmf$prob_weights, tt)

	diagnostic(paste("T(p) =", calculate_tp(phi_X, phi_W, sqrt_psi_W, weight)))
	diagnostic(paste("Penalties", calculate_penalties(phi_X, phi_W)))

	
	#--------------------------------------------------------------------------#
	# Convert to nice formats and return results
	#--------------------------------------------------------------------------#

	# Convert back to normal vectors
	x_sol <- min_var_sol$par
	p_sol <- c( x_sol[1:m-1], 1 - sum(x_sol[1:m-1]))
	theta_sol <- x_sol[m:(2 * m - 1)]
	# simple_sol <- simplify_pmf(theta_sol, p_sol)
	# p_sol <- simple_sol$ProbWeights
	# theta_sol <- simple_sol$Support

	x_min_tp <- min_tp_sol$par
	p_min_tp <- c( x_min_tp[1:m-1], 1 - sum(x_min_tp[1:m-1]))
	theta_min_tp <- x_min_tp[m:(2 * m - 1)]

	# simple_min_tp <- simplify_pmf(theta_min_tp, p_min_tp)
	# p_min_tp <- simple_min_tp$ProbWeights
	# theta_min_tp <- simple_min_tp$Support

	list("support" = theta_sol, 
		 "probweights" = p_sol, 
		 "support_min_tp" = theta_min_tp,
		 "probweights_min_tp" = p_min_tp,
		 "phi_W" = phi_W,
		 "var_opt_results" = min_var_sol,
		 "tp_opt_results" = min_tp_sol)
}

x_to_pmf <- function(x) {
	m <- (length(x) + 1) / 2
	probweights <- c( x[ 1:m - 1 ], 1 - sum( x[ 1:m - 1 ] ) )
	support <- x[ m:(2 * m - 1) ]

	list(support = support, prob_weights = probweights)
}

pmf_to_x <- function(theta, p) {
	m <- length(theta)

	if (m == 1){
		x <- theta
	} else {
		x <- c(p[1:(m-1)], theta)	
	}
	
	x
}

is_feasible <- function(xj, pj, A, B, nonlcon) {
	flag <- TRUE
	x <- pmf_to_x(xj, pj)

	lin_tol <- 1e-5
	if (any(A %*% x - B > lin_tol)) {
		flag <- FALSE
	}

	non_lin_tol <- 1e-5
	cons <- nonlcon(x)
	if (any(cons$c > non_lin_tol)) {
		flag <- FALSE
	}

	flag
}

tp_objective <- function(x, phi_W, sqrt_psi_W, weight) {
	X_pmf <- x_to_pmf(x)

	# Calculate phi_X
	tt <- phi_W$t.values
	phi_X <- ComputePhiPmf(X_pmf$support, X_pmf$prob_weights, tt)

	tp <- calculate_tp(phi_X, phi_W, sqrt_psi_W, weight)
	penalties <- calculate_penalties(phi_X, phi_W)

	penalty_scale <- 500
	tp + penalty_scale * sum(penalties)
}

calculate_tp <- function(phi_X, phi_W, sqrt_psi_W, weight){
	
	tt <- phi_W$t.values
	dt <- tt[2] - tt[1]
	# integrand <- abs(phi_W$complex - sqrt_psi_W * phi_X / Mod(phi_X))^2 * weight

	integrand <- abs(phi_W$complex * Mod(phi_X) - sqrt_psi_W * phi_X)^2 * weight

	tp <- dt * sum(integrand)

	tp
}

calculate_penalties <- function(phi_X, phi_W) {
	penalty1 = sum(abs(Re(phi_X) * phi_W$im - Im(phi_X) * phi_W$re))

	mod_phi_U = phi_W$norm / Mod(phi_X)
	penalty2 = sum(mod_phi_U[mod_phi_U > 1] - 1)

	c(penalty1, penalty2)
}

var_objective <- function(x){
	X_pmf <- x_to_pmf(x)
	p <- X_pmf$prob_weights
	theta <- X_pmf$support
	mean <- sum(p*theta)
	var <- sum(p*(theta - mean)^2)
	return(var)
}

constraints <- function(x, phi_W, sqrt_psi_W, weight, tp_max, penalties_max){
	X_pmf <- x_to_pmf(x)
	probweights <- X_pmf$prob_weights
	support <- X_pmf$support

	tt <- phi_W$t.values
	# Calculate phi_X
	phi_X <- ComputePhiPmf(support, probweights, tt)

	tp <- calculate_tp(phi_X, phi_W, sqrt_psi_W, weight)
	const1 <- tp - tp_max
	const23 <- calculate_penalties(phi_X, phi_W) - penalties_max
	
	list(ceq = NULL, c = c(const1, const23))
}

simplify_pmf <- function(theta, p, zero_tol = 1e-3, adj_tol = 1e-3){
	
	# Remove ps that are too small
	theta <- theta[p > zero_tol]
	p <- p[p > zero_tol]
	p <- p / sum(p)

	# Combine thetas that are too close together
	if (length(p) > 1){
		looping <- TRUE
	} else {
		looping <- FALSE
	}

	i <- 1
	while (looping){
		if (theta[i+1] - theta[i] > adj_tol){
			i <- i + 1
		} else {
			theta[i] <- (p[i] * theta[i] + p[i + 1] * theta[i + 1]) / 
						(p[i] + p[i + 1])
			theta <- theta[- (i + 1)]
			p[i] <- p[i] + p[i + 1]
			p <- p[- (i + 1)]
		}

		if (i >= length(p)){
			looping <- FALSE
		}
	}

	return(list("Support" = theta, "ProbWeights" = p))
}

compute_sqrt_psi_W <- function(tt, W){
	n <- length(W)

	WW <- outer(W, W, "-")
	WW <- c(WW[!(col(WW) == row(WW))])

	OO <- WW %o% tt

	re_psi_W <- colSums(cos(OO)) / (n*(n-1))
	im_psi_W <- colSums(sin(OO)) / (n*(n-1))

	sqrt_psi_W <- sqrt(sqrt(re_psi_W^2 + im_psi_W^2))

	sqrt_psi_W
}

print_diagnostic <- function(message, show_diagnostics){
	if (show_diagnostics){
		print(message)
	}
}

create_bound_matrices <- function(W, m){
	# pj non-negative
	A <- matrix(0, nrow = 2*m+1, ncol = 2*m-1)
	A[1:m-1, 1:m-1] <- -diag(m - 1)
	B <- numeric(2 * m - 1)
	# pj sum to less than 1
	A[m, 1:m-1] = matrix(1, nrow=1, ncol = m-1)
	B[m] = 1
	# thetaj are increasing
	if (m > 1){
		for (i in 1:(m-1)){
			A[m+i, (m+i-1):(m+i)] <- c(1, -1)	
		}
	}
	# min(W) < thetaj < max(W)
	A[2*m, m] <- -1
	A[2*m+1, 2*m-1] <- 1
	B[2*m] <- -min(W)
	B[2*m+1] <- max(W)

	list(A = A, B = B)
}

