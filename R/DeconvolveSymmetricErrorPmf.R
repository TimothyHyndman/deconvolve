#' Estimate Discrete Deconvolution of Data With Symmetric Error
#' 
#' Estimates deconvolution of data \eqn{W = X + U}, where \eqn{U} is symmetric, 
#' with a probability mass function.
#' 
#' Details here
#' 
#' @param W A vector of the contaminated data
#' @param m The number of support points to use in finding the Pmf
#' 
#' @return A list with components:
#' 	\item{support}{The support of the deconvolved distribution}
#' 	\item{probweights}{The sizes of the probability mass at each respective 
#' 					   point of support.}
#' 	\item{phi.W}{The empirical distribution of W and the t values on which it is 
#' 				 calculated}
#' 
#' @example man/examples/SymmetricError_eg.R
#' 
#' @export

# There are two ways we can solve this optimization problem. 
# 1. Find T* = min(T(p)) then minimize variance under constraint that T(p) = T*
#
# 2. Minimize var + lambda*T(p) for some suitable lambda.
#
# Option 1 is preferable I think, but harder to implement since T(p) = T* is a 
# very tight constraint

DeconvolveSymmetricErrorPmf <- function(W, m = 10){

	n <- length(W)

	#--------------------------------------------------------------------------#
	# Pre-calculate PhiW and weight w(t)
	#--------------------------------------------------------------------------#
	# Calculate phi.W on [-8,8] so we can find t*
	tt.length <- 100
	a <- -8
	b <- 8
	tt <- seq( a, b, length.out = tt.length )
	phi.W <- ComputePhiEmp(W, tt)

	# Calculate t*
	tmp <- tt[phi.W$norm < n^(-0.25)]
	if ( length( tmp[ tmp > 0] ) == 0 ){
		t.star <- max( tt )
	} else {
		t.star <- min( tmp[ tmp > 0 ] )
	}

	# Calculate phi.W on [-t*,t*]
	tt.new.length <- 100
	tt.new <- seq( -t.star, t.star, length.out = tt.new.length )
	phi.W <- ComputePhiEmp(W, tt.new)

	# Calculate weight w(t) on [-t*, t*]
	weight <- KernelWeight(tt.new)

	#--------------------------------------------------------------------------#
	# Solve optimization problem to find PMF
	#--------------------------------------------------------------------------#

	# Initial point for optimization routine
	theta0 <- seq( min(W), max(W), length.out = m )
	p0 <- numeric( m - 1 ) + 1/m
	x0 <- c(p0, theta0)

	# Set up constraints in form Ax \geq b
	A <- matrix(0, nrow = m, ncol = 2*m-1)
	A[1:m-1, 1:m-1] <- diag( m - 1 )   #pj non-negative
	b = numeric(m)
	A[m, 1:m-1] = matrix(-1, nrow=1, ncol = m-1) #pj sum to less than 1
	b[m] = -1
	

	#----------------------------------#
	# Option 1
	#----------------------------------#

	# Minimize T(p)

	# min.tp.sol <- stats::constrOptim(x0, CalculateTp, NULL, A, b, phi.W = phi.W, 
	# 						  weight = weight)
	
	# Minimize the variance subject to T(theta,p) not increasing
	# This is currently not working well as optim can't find any feasible 
	# points other than the initial.

	# min.var.sol <- stats::constrOptim(min.tp.sol$par, CalculateVarwithConstr, NULL, A, 
	# 						   b, max.tp = min.tp.sol$value, phi.W = phi.W, 
	# 						   weight = weight)

	# min.var.sol <- stats::constrOptim(min.tp.sol$par, CombinedObjective, NULL, 
	# 								  A, b, phi.W = phi.W, weight = weight)

	# x.sol <- min.var.sol$par
	# p.sol <- c( x.sol[1:m-1], 1 - sum(x.sol[1:m-1]))
	# theta.sol <- x.sol[m:(2 * m - 1)]

	# return(list("support" = theta.sol, "probweights" = p.sol))


	#------------------------------------#
	# Option 2
	#------------------------------------#

	# min.comb.sol <- stats::constrOptim(x0, CombinedObjective, NULL, A, b, 
	# 							phi.W = phi.W, weight = weight)

	# x.sol <- min.comb.sol$par
	# p.sol <- c( x.sol[1:m-1], 1 - sum(x.sol[1:m-1]))
	# theta.sol <- x.sol[m:(2 * m - 1)]

	#------------------------------------#
	# Option 3
	#------------------------------------#

	min.tp.sol <- stats::constrOptim(x0, CalculateTp, NULL, A, b, phi.W = phi.W, 
							  weight = weight)


	lb <- numeric(2*m-1)
	lb[m:(2 * m - 1)] <- min(W)
	ub <- numeric(2*m-1) + 1
	ub[m:(2 * m - 1)] <- max(W)

	opts <- list("algorithm"="NLOPT_LD_SLSQP", "maxeval" = 1e4)
	min.var.sol <- nloptr::nloptr(min.tp.sol$par, eval_f = CalculateVar, 
								  eval_grad_f = CalculateVarGrad, lb = lb, 
								  ub = ub, eval_g_ineq = Constraints, 
								  phi.W = phi.W, weight = weight,
								  tp.max = min.tp.sol$value, opts = opts)

	#--------------------------------------------------------------------------#
	# Return
	#--------------------------------------------------------------------------#

	return(list("support" = theta.sol, "probweights" = p.sol, "phi.W" = phi.W))
}

CalculateTp <- function(x, phi.W, weight){
	m <- (length(x) + 1) / 2

	p <- c( x[ 1:m - 1 ], 1 - sum( x[ 1:m - 1 ] ) )
	theta <- x[ m:(2 * m - 1) ]
	tt <- phi.W$t.values

	# Calculate phi.ptheta
	phi.ptheta <- ComputePhiPmf(theta, p, tt)	

	# Calculate integral
	fred <- phi.W$complex - phi.W$norm * phi.ptheta / Mod(phi.ptheta)
	integrand <- abs(fred)^2 * weight
	dt <- tt[2] - tt[1]
	tp <- dt * sum(integrand)

	return(tp)
}

CalculateVar <- function(x, phi.W, weight, tp.max){
	m <- (length(x) + 1) / 2
	p <- c( x[ 1:(m - 1) ], 1 - sum( x[ 1:(m - 1) ] ) )
	theta <- x[ m:(2 * m - 1) ]
	mean <- sum(p*theta)
	var <- sum(p*(theta - mean)^2)
	return(var)
}

CalculateVarGrad <- function(x, phi.W, weight, tp.max){
	var.x <- CalculateVar(x, phi.W, weight, tp.max)
	dx <- 0.01
	grad <- numeric(length(x))

	for (i in 1:length(x)){
		y <- x
		y[i] <- x[i] + dx
		var.x1 <- CalculateVar(x1, phi.W, weight, tp.max)
		grad[i] <- (var.x1 - var.x) / dx	
	}
	
	return(grad)
}

# CalculateVarwithConstr <- function(x, max.tp, phi.W, weight){
# 	#Check that tp hasn't gotten any bigger
# 	tp <- CalculateTp(x, phi.W, weight)
# 	if (tp > max.tp){
# 		return(Inf)
# 	}

# 	# All constraints are met so calculate variance
# 	var <- CalculateVar(x)
# 	return(var)
# }

# CombinedObjective <- function(x, phi.W, weight, lambda = 10000){
# 	CalculateVar(x) + lambda * CalculateTp(x, phi.W, weight)
# }

Constraints <- function(x, phi.W, weight, tp.max){
	tp <- CalculateTp(x, phi.W, weight)
	const1 <- tp - tp.max

	m <- (length(x) + 1) / 2
	p <- c( x[ 1:m - 1 ], 1 - sum( x[ 1: (m - 1) ] ) )
	theta <- x[ m:(2 * m - 1) ]

	const2 <- sum(p < 0)
	const3 <- sum(p > 1)

	return(c(const1, const2, const3))
}

ConstraintsGrad <- function(x, phi.W, weight, tp.max){
	tp.x <- CalculateTp(x, phi.W, weight, tp.max)
	dx <- 0.01
	grad <- numeric(length(x))
}