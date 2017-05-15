#' @export

DiscreteDeconvolutionSymmetricError <- function(W, weight.type = 'Epanechnikov'){

	n <- length(W)

	# Pre-calculate PhiW

	tt.length <- 100
	a <- -8
	b <- 8
	tt <- seq( a, b, length.out = tt.length )

	tt.new.length <- 100
	phi.W <- ComputePhi(W, tt, tt.new.length )	#Currently not calculating sqrt psi hat W

	# Choose kernel Weight
	weight <- KernelWeight(weight.type, phi.W$t.values)

	#Choose theta and p to minimize T(theta,p) = int |psi.hat(t) - blahblah | w(t) dt
	# under the constraint that blah
	m <- 5

	theta0 <- seq( min(W), max(W), length.out = m )
	p0 <- numeric( m - 1 ) + 1/m
	x0 <- c(p0, theta0)

	A <- matrix(0, nrow = m, ncol = 2*m-1)
	A[1:m-1, 1:m-1] <- diag( m - 1 )   #pj non-negative
	ci = numeric(m)
	A[m, 1:m-1] = matrix(-1, nrow=1, ncol = m-1) #pj sum to less than 1
	ci[m] = -1
	
	min.sol <- constrOptim(x0, CalculateTp, NULL, A, ci, phi.W = phi.W, weight = weight)
	x.sol <- min.sol$par

	p.sol <- c( x.sol[1:m-1], 1 - sum(x.sol[1:m-1]))
	theta.sol <- x.sol[m:(2 * m - 1)]

	return( list( "ProbWeights" = p.sol, "Support" = theta.sol) )

	# Minimize the variance subject to T(theta,p) not increasing

}

CalculateTp <- function(x, phi.W, weight){
	m <- (length(x) + 1) / 2

	p <- c( x[ 1:m - 1 ], 1 - sum( x[ 1:m - 1 ] ) )
	theta <- x[ m:(2 * m - 1) ]
	tt <- phi.W$t.values

	# Compute phi.ptheta
	phi.ptheta <- p %*% exp(complex(imaginary = 1) * matrix(theta,ncol = 1) %*% tt )

	sqrt.psi.hat.W <- phi.W$norm 	# This is calculated differently in the paper, fix later
	hat.phi.W <- complex(real = phi.W$re, imaginary = phi.W$im)	#Check this is how it's calculated in the paper

	fred <- hat.phi.W - sqrt.psi.hat.W * phi.ptheta / abs(phi.ptheta)

	integrand <- abs(fred)^2 * weight
	dt <- tt[2] - tt[1]
	tp <- dt * sum(integrand)

	return(tp)
}