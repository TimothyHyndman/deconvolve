DeconErrSymPmf <- function(W, m, kernel_type, n_tp_iter = 10, n_var_iter = 10, 
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
	t_star <- find_t_cutoff(phi_W$norm, tt)

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

	matrices <- create_bound_matrix(W, m)

	# ------------------
	# Min T(p)
	# ------------------
	tp_obj_min <- Inf
	for (i in 1:n_tp_iter) {
		theta0 <- sort(stats::runif(m, min = min(W), max = max(W)))
		p0 <- stats::runif(m, min = 0, max = 1)
		p0 <- p0 / sum(p0)
		x0 <- theta_p_to_x(theta0, p0)

		control <- list(maxit = 100000000)
		optim_result <- optim(x0,
							  tp_objective, 
							  control = control,
							  method = "Nelder-Mead",
							  phi_W = phi_W, 
							  sqrt_psi_W = sqrt_psi_W, 
							  weight = weight,
							  A = matrices$A, 
							  B = matrices$B)

		diagnostic(optim_result$convergence)

		if (optim_result$value < tp_obj_min) {
			optim_result_min <- optim_result
			tp_obj_min <- optim_result$value

			diagnostic(optim_result$value)
		}
	}

	tp_optim_result_min <- optim_result_min
	x_sol <- tp_optim_result_min$par
	theta_sol <- x_to_theta(x_sol)
	p_sol <- x_to_p(x_sol)
	
	# ------------------
	# Min Var
	# ------------------
	phi_X <- ComputePhiPmf(theta_sol, p_sol, tt_new)
	tp_max <- calculate_tp(phi_X, phi_W, sqrt_psi_W, weight)
	penalties_max <- calculate_penalties(phi_X, phi_W)

	diagnostic(tp_max)
	diagnostic(penalties_max)

	optim_result <- optim(x_sol,
						  var_objective, 
						  control = control,
						  method = "Nelder-Mead",
						  phi_W = phi_W, 
						  tp_max = tp_max,
						  penalties_max = penalties_max,
						  sqrt_psi_W = sqrt_psi_W, 
						  weight = weight,
						  A = matrices$A, 
						  B = matrices$B)

	var_obj_min <- optim_result$value
	diagnostic(var_obj_min)

	for (i in 1:n_var_iter) {
		theta0 <- sort(stats::runif(m, min = min(W), max = max(W)))
		p0 <- stats::runif(m, min = 0, max = 1)
		p0 <- p0 / sum(p0)
		x0 <- theta_p_to_x(theta0, p0)

		optim_result <- optim(x0,
							  var_objective, 
							  control = control,
							  method = "Nelder-Mead",
							  phi_W = phi_W, 
							  tp_max = tp_max,
							  penalties_max = penalties_max,
							  sqrt_psi_W = sqrt_psi_W, 
							  weight = weight,
							  A = matrices$A, 
							  B = matrices$B)

		diagnostic(optim_result$convergence)

		if (optim_result$value < var_obj_min) {
			optim_result_min <- optim_result
			var_obj_min <- optim_result$value

			diagnostic(optim_result$value)
		}
	}

	x_sol <- optim_result_min$par
	theta_sol <- x_to_theta(x_sol)
	p_sol <- x_to_p(x_sol)

	#--------------------------------------------------------------------------#
	# Convert to nice formats and return results
	#--------------------------------------------------------------------------#

	list("support" = theta_sol, 
		 "probweights" = p_sol, 
		 "phi_W" = phi_W,
		 "var_opt_results" = optim_result_min,
		 "tp_opt_results" = tp_optim_result_min)
}

x_to_theta <- function(x) {
	m <- (length(x) + 1) / 2
	
	x[ m:(2 * m - 1) ]
}

x_to_p <- function(x) {
	m <- (length(x) + 1) / 2
	
	c( x[ 1:m - 1 ], 1 - sum( x[ 1:m - 1 ] ) )
}

theta_p_to_x <- function(theta, p) {
	m <- length(theta)
	c(p[1:m-1], theta)
}

tp_objective <- function(x, phi_W, sqrt_psi_W, weight, A, B) {
	theta <- x_to_theta(x)
	p <- x_to_p(x)

	# Calculate phi_X
	tt <- phi_W$t.values
	phi_X <- ComputePhiPmf(theta, p, tt)

	tp <- calculate_tp(phi_X, phi_W, sqrt_psi_W, weight)
	penalties <- calculate_penalties(phi_X, phi_W)

	cliff = 0
	if (!is_valid_pmf(x, A, B)) {
		cliff = 1e20
	}

	tp + sum(penalties) + cliff
}

calculate_tp <- function(phi_X, phi_W, sqrt_psi_W, weight){
	
	tt <- phi_W$t.values
	dt <- tt[2] - tt[1]
	integrand <- abs(phi_W$complex - sqrt_psi_W * phi_X / Mod(phi_X))^2 * weight
	tp <- dt * sum(integrand)

	tp
}

calculate_penalties <- function(phi_X, phi_W) {
	penalty1 = sum(abs(Re(phi_X) * phi_W$im - Im(phi_X) * phi_W$re))

	mod_phi_U = phi_W$norm / Mod(phi_X)
	penalty2 = sum(mod_phi_U[mod_phi_U > 1])

	c(penalty1, penalty2)
}

var_objective <- function(x, phi_W, tp_max, penalties_max, sqrt_psi_W, weight, A, B){
	p <- x_to_p(x)
	theta <- x_to_theta(x)

	mean <- sum(p*theta)
	var <- sum(p*(theta - mean)^2)

	tt <- phi_W$t.values
	phi_X <- ComputePhiPmf(theta, p, tt)
	penalties <- max(0, calculate_penalties(phi_X, phi_W) - penalties_max)
	tp_penalty <- max(0, calculate_tp(phi_X, phi_W, sqrt_psi_W, weight) - tp_max)

	cliff <- 0

	if (!is_valid_pmf(x, A, B)) {
		cliff <- 1e20
	}

	if (any(c(tp_penalty, penalties) > 0)) {
		penalty_scale <- 1e3
		penalty_height <- 1e3
		cliff <- cliff + penalty_height + penalty_scale*(tp_penalty + sum(penalties))
	}

	var + cliff
}

is_valid_pmf <- function(x, A, B) {
	flag = TRUE

	if (any(A %*% x - B > 0)) {
		flag = FALSE
	}

	flag
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

create_bound_matrix <- function(W, m) {
	# Set up constraints in form Ax \geq B

	# pj non-negative
	A <- matrix(0, nrow = 2*m+1, ncol = 2*m-1)
	A[1:m-1, 1:m-1] <- -diag(m - 1)
	B <- numeric(2 * m - 1)
	# pj sum to less than 1
	A[m, 1:m-1] = matrix(1, nrow=1, ncol = m-1)
	B[m] = 1
	# thetaj are increasing
	for (i in 1:(m-1)){
		A[m+i, (m+i-1):(m+i)] <- c(1, -1)	
	}
	# min(W) < thetaj < max(W)
	A[2*m, m] <- -1
	A[2*m+1, 2*m-1] <- 1
	B[2*m] <- -min(W)
	B[2*m+1] <- max(W)

	list(A = A, B = B)
}