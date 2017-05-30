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

	# ---------
	# Min T(p)
	# ---------

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

	# Initial point for optimization routine
	theta0 <- seq( min(W) + 0.01, max(W) - 0.01, length.out = m )
	p0 <- numeric( m - 1 ) + 1/m
	x0 <- c(p0, theta0)

	min.tp.sol <- stats::constrOptim(x0, CalculateTp, NULL, A, b, phi.W = phi.W, 
							  weight = weight)

	# ---------
	# Min Var
	# ---------

	# Ax <= B
	# pj non-negative
	A <- matrix(0, nrow = 2*m-1, ncol = 2*m-1)
	A[1:m-1, 1:m-1] <- - diag(m - 1)
	B = numeric(2 * m - 1)
	# pj sum to less than 1
	A[m, 1:m-1] = matrix(1, nrow=1, ncol = m-1)
	B[m] = 1
	#thetaj are increasing
	for (i in 1:(m-1)){
		A[m+i, (m+i-1):(m+i)] <- c(1, -1)	
	}
	B <- matrix(B, ncol = 1)
	# lb <= x <= ub
	lb <- numeric(2*m-1)
	ub <- numeric(2*m-1) + 1
	lb[m:(2 * m - 1)] <- min(W)
	ub[m:(2 * m - 1)] <- max(W)
	lb <- matrix(lb, ncol = 1)
	ub <- matrix(lb, ncol = 1)
	ConFun <- function(x){
		Constraints(x, phi.W, weight, min.tp.sol$value)
	}

	x0 <- min.tp.sol$par
	min.var.sol <- NlcOptim::solnl(x0, CalculateVar, ConFun, A, B, lb, ub)	#lb and ub 
	# aren't working yet :(
	
	#--------------------------------------------------------------------------#
	# Convert to nice formats and return results
	#--------------------------------------------------------------------------#

	# Convert back to normal vectors
	x.sol <- min.var.sol$solution
	p.sol <- c( x.sol[1:m-1], 1 - sum(x.sol[1:m-1]))
	theta.sol <- x.sol[m:(2 * m - 1)]
	simple.sol <- SimplifyPmf(theta.sol, p.sol)
	p.sol <- simple.sol$ProbWeights
	theta.sol <- simple.sol$Support

	# See effect of min variance.
	x.min.tp <- min.tp.sol$par
	p.min.tp <- c( x.min.tp[1:m-1], 1 - sum(x.min.tp[1:m-1]))
	theta.min.tp <- x.min.tp[m:(2 * m - 1)]
	simple.min.tp <- SimplifyPmf(theta.min.tp, p.min.tp)
	p.min.tp <- simple.min.tp$ProbWeights
	theta.min.tp <- simple.min.tp$Support



	# Plot for diagnostics
	df.sol <- data.frame(p.sol, theta.sol)
	df.min.tp <- data.frame(p.min.tp, theta.min.tp)
	plot <- ggplot() + geom_point(data = df.min.tp, 
								  aes(theta.min.tp, p.min.tp), 
								  color = "magenta") + 
				 	   geom_point(data = df.sol, 
				 	   			  aes(theta.sol, p.sol), 
				 	   			  color = "blue")

	return(list("plot" = plot, 
				"support" = theta.sol, 
				"probweights" = p.sol, 
				"phi.W" = phi.W,
				"var.opt.results" = min.var.sol,
				"tp.opt.results" = min.tp.sol))
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
CalculateVar <- function(x){
	m <- (length(x) + 1) / 2
	p <- c( x[ 1:(m - 1) ], 1 - sum( x[ 1:(m - 1) ] ) )
	theta <- x[ m:(2 * m - 1) ]
	mean <- sum(p*theta)
	var <- sum(p*(theta - mean)^2)
	return(var)
}

#' @export
Constraints <- function(x, phi.W, weight, tp.max){
	tp <- CalculateTp(x, phi.W, weight)
	const1 <- tp - tp.max
	return(list(ceq = NULL, c = const1))
}

#' @export
SimplifyPmf <- function(theta, p, zero.tol = 1e-3, adj.tol = 1e-3){
	
	# Remove ps that are too small
	theta <- theta[p > zero.tol]
	p <- p[p > zero.tol]
	p <- p / sum(p)

	# Combine thetas that are too close together
	if (length(p) > 1){
		looping <- T
	} else {
		looping <- F
	}

	i <- 1
	while (looping){
		if (theta[i+1] - theta[i] > adj.tol){
			i <- i + 1
		} else {
			theta[i] <- (p[i] * theta[i] + p[i + 1] * theta[i + 1]) / 
						(p[i] + p[i + 1])
			theta <- theta[- (i + 1)]
			p[i] <- p[i] + p[i + 1]
			p <- p[- (i + 1)]
		}

		if (i >= length(p)){
			looping <- F
		}
	}

	return(list("Support" = theta, "ProbWeights" = p))
}