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
	theta0 <- seq( min(W) + 0.01, max(W) - 0.01, length.out = m )
	p0 <- numeric( m - 1 ) + 1/m
	x0 <- c(p0, theta0)

	# Set up constraints in form Ax \geq b
	A <- matrix(0, nrow = m-1 + 1 + 2*m, ncol = 2*m-1)
	A[1:m-1, 1:m-1] <- diag( m - 1 )   # pj non-negative
	b = numeric(m + 2 * m)
	A[m, 1:m-1] = matrix(-1, nrow=1, ncol = m-1) # pj sum to less than 1
	b[m] = -1
	# thetaj greater than min(W)
	A[(m + 1):(2 * m), m:(2 * m - 1)] <- diag(m)
	b[(m+1):(2*m)] <- min(W)
	# thetaj less than max(W)
	A[(2 * m + 1):(3 * m), m:(2 * m - 1)] <- -diag(m)
	b[(2*m+1):(3*m)] <- -max(W)


	min.tp.sol <- stats::constrOptim(x0, CalculateTp, NULL, A, b, phi.W = phi.W, 
							  weight = weight)


	lb <- numeric(2*m-1)
	lb[m:(2 * m - 1)] <- min(W)
	ub <- numeric(2*m-1) + 1
	ub[m:(2 * m - 1)] <- max(W)

	# opts <- list("algorithm"="NLOPT_LD_SLSQP", "maxeval" = 1e3, 
	# 			 "xtol_rel" = 1e-4)
	local.opts <- list( "algorithm" = "NLOPT_LD_MMA", 
					  	"xtol_rel" = 1e-4 )
	opts <- list( "algorithm"="NLOPT_LD_AUGLAG", 
				  "maxeval" = 1e4, 
				  "xtol_rel" = 1e-4, 
				  "local_opts" = local.opts )
	
	min.var.sol <- nloptr::nloptr(min.tp.sol$par, eval_f = CalculateVar, 
								  eval_grad_f = CalculateVarGrad, lb = lb, 
								  ub = ub, eval_g_ineq = Constraints, 
								  eval_jac_g_ineq = ConstraintsGrad,
								  opts = opts, phi.W = phi.W, weight = weight,
								  tp.max = min.tp.sol$value)

	x.sol <- min.var.sol$solution
	p.sol <- c( x.sol[1:m-1], 1 - sum(x.sol[1:m-1]))
	theta.sol <- x.sol[m:(2 * m - 1)]
	plot(theta.sol, p.sol)

	return(list("support" = theta.sol, "probweights" = p.sol, "phi.W" = phi.W))
}
#' @export
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
#' @export
CalculateVar <- function(x, phi.W, weight, tp.max){
	m <- (length(x) + 1) / 2
	p <- c( x[ 1:(m - 1) ], 1 - sum( x[ 1:(m - 1) ] ) )
	theta <- x[ m:(2 * m - 1) ]
	mean <- sum(p*theta)
	var <- sum(p*(theta - mean)^2)
	return(var)
}
#' @export
CalculateVarGrad <- function(x, phi.W, weight, tp.max){
	numDeriv::grad(CalculateVar, x, method = "simple", phi.W = phi.W, weight = weight, 
				   tp.max = tp.max)
}

#' @export
Constraints <- function(x, phi.W, weight, tp.max){
	tp <- CalculateTp(x, phi.W, weight)
	const1 <- tp - tp.max	# Does giving it a little slack help?

	m <- (length(x) + 1) / 2
	p <- c( x[ 1:m - 1 ], 1 - sum( x[ 1: (m - 1) ] ) )
	theta <- x[ m:(2 * m - 1) ]

	const2 <- -sum(p[p < 0])
	const3 <- sum(p[p > 1] - 1)

	return(c(const1, const2, const3))
}
#' @export
ConstraintsGrad <- function(x, phi.W, weight, tp.max){
	numDeriv::jacobian(Constraints, x, method = "simple", phi.W = phi.W, weight = weight,
					   tp.max = tp.max)
}