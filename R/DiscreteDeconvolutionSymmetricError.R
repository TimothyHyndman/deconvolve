#' Estimate Discrete Deconvolution of Data With Symmetric Error
#' 
#' Estimates deconvolution of data \eqn{W = X + U}, where \eqn{U} is symmetric, 
#' with a probability mass function.
#' 
#' Details here
#' 
#' @param W A vector of the contaminated data
#' 
#' @return A list with components:
#' 	\item{support}{The support of the deconvolved distribution}
#' 	\item{probweights}{The sizes of the probability mass at each respective 
#' 					   point of support.}
#' 
#'
#' @export

# There are two ways we can solve this optimization problem. 
# 1. Find T* = min(T(p)) then minimize variance under constraint that T(p) = T*
#
# 2. Minimize var + lambda*T(p) for some suitable lambda.
#
# Option 1 is preferable I think, but harder to implement.

DiscreteDeconvolutionSymmetricError <- function(W){

	n <- length(W)

	# Pre-calculate PhiW and weight w(t)

	tt.length <- 100
	a <- -8
	b <- 8
	tt <- seq( a, b, length.out = tt.length )

	tt.new.length <- 100
	phi.W <- ComputePhi(W, tt, tt.new.length )	#Currently not calculating sqrt 
												#psi hat W
	weight <- KernelWeight(phi.W$t.values)

	#Choose theta and p to minimize T(theta,p) = int |psi.hat(t) - blahblah | 
	# w(t) dt under the constraint that blah
	m <- 10

	theta0 <- seq( min(W), max(W), length.out = m )
	p0 <- numeric( m - 1 ) + 1/m
	x0 <- c(p0, theta0)

	A <- matrix(0, nrow = m, ncol = 2*m-1)
	A[1:m-1, 1:m-1] <- diag( m - 1 )   #pj non-negative
	ci = numeric(m)
	A[m, 1:m-1] = matrix(-1, nrow=1, ncol = m-1) #pj sum to less than 1
	ci[m] = -1
	

	#--------------------------------------------------------------------------#
	# Option 1
	#--------------------------------------------------------------------------#

	# Minimize T(p)

	# min.tp.sol <- constrOptim(x0, CalculateTp, NULL, A, ci, phi.W = phi.W, 
	# 						  weight = weight)
	
	# Minimize the variance subject to T(theta,p) not increasing
	# This is currently not working well as optim can't find any feasible 
	# points other than the initial.

	# min.var.sol <- constrOptim(min.tp.sol$par, CalculateVarwithConstr, NULL, A, 
	# 						   ci, max.tp = min.tp.sol$value, phi.W = phi.W, 
	# 						   weight = weight)

	# x.sol <- min.var.sol$par
	# p.sol <- c( x.sol[1:m-1], 1 - sum(x.sol[1:m-1]))
	# theta.sol <- x.sol[m:(2 * m - 1)]

	# return(list("support" = theta.sol, "probweights" = p.sol))


	#--------------------------------------------------------------------------#
	# Option 2
	#--------------------------------------------------------------------------#

	min.comb.sol <- constrOptim(x0, CombinedObjective, NULL, A, ci, 
								phi.W = phi.W, weight = weight, lambda = 10000)

	x.sol <- min.comb.sol$par
	p.sol <- c( x.sol[1:m-1], 1 - sum(x.sol[1:m-1]))
	theta.sol <- x.sol[m:(2 * m - 1)]

	return(list("support" = theta.sol, "probweights" = p.sol))
}

CalculateTp <- function(x, phi.W, weight){
	m <- (length(x) + 1) / 2

	p <- c( x[ 1:m - 1 ], 1 - sum( x[ 1:m - 1 ] ) )
	theta <- x[ m:(2 * m - 1) ]
	tt <- phi.W$t.values

	# Compute phi.ptheta
	phi.ptheta <- p %*% exp( complex(imaginary = 1) * matrix(theta, ncol = 1) %*% tt )

	sqrt.psi.hat.W <- phi.W$norm 	# This is calculated differently in the paper, fix later
	hat.phi.W <- complex(real = phi.W$re, imaginary = phi.W$im)	#Check this is how it's calculated in the paper

	fred <- hat.phi.W - sqrt.psi.hat.W * phi.ptheta / abs(phi.ptheta)

	integrand <- abs(fred)^2 * weight
	dt <- tt[2] - tt[1]
	tp <- dt * sum(integrand)

	return(tp)
}

CalculateVarwithConstr <- function(x, max.tp, phi.W, weight){
	#Check that tp hasn't gotten any bigger
	tp <- CalculateTp(x, phi.W, weight)
	if (tp > max.tp){
		return(Inf)
	}

	# All constraints are met so calculate variance
	var <- CalculateVar(x)
	return(var)
}

CalculateVar <- function(x){
	m <- (length(x) + 1) / 2
	p <- c( x[ 1:m - 1 ], 1 - sum( x[ 1:m - 1 ] ) )
	theta <- x[ m:(2 * m - 1) ]
	mean <- sum(p*theta)
	var <- sum(p*(theta - mean)^2)
	return(var)
}

CombinedObjective <- function(x, phi.W, weight, lambda = 10000){
	CalculateVar(x) + lambda * CalculateTp(x, phi.W, weight)
}