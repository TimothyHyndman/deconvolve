DeconErrSymPmf <- function(W, m, kernel_type, n_tp_iter = 5, n_var_iter = 2, 
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

	# ------------------
	# Min T(p)
	# ------------------

	for (i in 1:n_tp_iter) {
		theta0 <- sort(stats::runif(m, min = min(W), max = max(W)))
		p0 <- stats::runif(m, min = 0, max = 1)
		p0 <- p0 / sum(p0)
		x0 <- theta_p_to_x(theta0, p0)

		optim_result <- optim(x0,
							  tp_objective, 
							  phi_W = phi_W, 
							  sqrt_psi_W = sqrt_psi_W, 
							  weight = weight,
							  W = W)

		if (optim_result$value < tp_min) {
			tp_min <- optim_result$value
			x_sol <- optim_result$par

			diagnostic(optim_result$value)
			p_sol <- x_to_p(x_sol)
			theta_sol <- x_to_theta(x_sol)

			diagnostic(p_sol)
			diagnostic(theta_sol)
		}
	}
	

	# ------------------
	# Min Var
	# ------------------

	
	#--------------------------------------------------------------------------#
	# Convert to nice formats and return results
	#--------------------------------------------------------------------------#

	# Convert back to normal vectors
	x_sol <- min_var_sol$par
	p_sol <- c( x_sol[1:m-1], 1 - sum(x_sol[1:m-1]))
	theta_sol <- x_sol[m:(2 * m - 1)]
	simple_sol <- simplify_pmf(theta_sol, p_sol)
	p_sol <- simple_sol$ProbWeights
	theta_sol <- simple_sol$Support

	x_min_tp <- min_tp_sol$par
	p_min_tp <- c( x_min_tp[1:m-1], 1 - sum(x_min_tp[1:m-1]))
	theta_min_tp <- x_min_tp[m:(2 * m - 1)]

	simple_min_tp <- simplify_pmf(theta_min_tp, p_min_tp)
	p_min_tp <- simple_min_tp$ProbWeights
	theta_min_tp <- simple_min_tp$Support

	list("support" = theta_sol, 
		 "probweights" = p_sol, 
		 "support_min_tp" = theta_min_tp,
		 "probweights_min_tp" = p_min_tp,
		 "phi_W" = phi_W,
		 "var_opt_results" = min_var_sol,
		 "tp_opt_results" = min_tp_sol)
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

tp_objective <- function(x, phi_W, sqrt_psi_W, weight, W) {
	theta <- x_to_theta(x)
	p <- x_to_p(x)

	# Calculate phi_X
	tt <- phi_W$t.values
	phi_X <- ComputePhiPmf(theta, p, tt)

	tp <- calculate_tp(phi_X, phi_W, sqrt_psi_W, weight)
	penalties <- calculate_penalties(phi_X, phi_W)

	cliff = 0
	if (!is_valid_pmf(theta, p, W)) {
		cliff = 1e20
	}
	tp + sum(penalties) + cliff
}

is_valid_pmf <- function(theta, p, W) {
	flag = TRUE

	if (any(p < 0)) {
		flag = FALSE
	}

	if (any(theta < min(W))) {
		flag = FALSE
	}

	if (any(theta > max(W))) {
		flag = FALSE
	}

	flag
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

var_objective <- function(x){
	m <- (length(x) + 1) / 2
	p <- c( x[ 1:(m - 1) ], 1 - sum( x[ 1:(m - 1) ] ) )
	theta <- x[ m:(2 * m - 1) ]
	mean <- sum(p*theta)
	var <- sum(p*(theta - mean)^2)
	return(var)
}

constraints <- function(x, phi_W, sqrt_psi_W, weight, tp_max, penalties_max){
	m <- (length(x) + 1) / 2
	probweights <- c( x[ 1:m - 1 ], 1 - sum( x[ 1:m - 1 ] ) )
	support <- x[ m:(2 * m - 1) ]
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
