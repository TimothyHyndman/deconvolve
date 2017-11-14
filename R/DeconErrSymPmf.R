DeconErrSymPmf <- function(W, m, n.iter.tp = 5, n.iter.var = 2, 
						   show.diagnostics = F){

	Diagnostic <- function(message){
		PrintDiagnostic(message, show.diagnostics)
	}

	if (m < 2){
		stop("m must be at least 2")
	}
	
	n <- length(W)

	#--------------------------------------------------------------------------#
	# Pre-calculate PhiW and weight w(t)
	#--------------------------------------------------------------------------#
	# Calculate phi.W on [-8,8] so we can find t*
	tt.length <- 100
	tt <- seq(-1, 8, length.out = tt.length)
	phi.W <- ComputePhiEmp(W, tt)

	# Calculate t*
	tmp <- tt[phi.W$norm < n^(-0.25)]
	if ( length( tmp[ tmp > 0] ) == 0 ){
		t.star <- max(tt)
	} else {
		t.star <- min( tmp[ tmp > 0 ] )
	}

	# Calculate phi.W on [-t*,t*]
	tt.new.length <- 100
	tt.new <- seq(-t.star, t.star, length.out = tt.new.length)
	phi.W <- ComputePhiEmp(W, tt.new)

	# Calculate weight w(t) on [-t*, t*]
	weight <- KernelWeight(tt.new)

	#--------------------------------------------------------------------------#
	# Solve optimization problem to find PMF
	#--------------------------------------------------------------------------#

	# ------------------
	# Setup for Min T(p)
	# ------------------
	# Set up constraints in form Ax \geq B
	
	# pj non-negative
	A.tp <- matrix(0, nrow = 2*m+1, ncol = 2*m-1)
	A.tp[1:m-1, 1:m-1] <- diag( m - 1 )   
	B.tp = numeric(2 * m - 1)
	# pj sum to less than 1
	A.tp[m, 1:m-1] = matrix(-1, nrow=1, ncol = m-1) 
	B.tp[m] = -1
	# thetaj are increasing
	for (i in 1:(m-1)){
		A.tp[m+i, (m+i-1):(m+i)] <- c(-1, 1)	
	}
	# min(W) < thetaj < max(W)
	A.tp[2*m, m] <- 1
	A.tp[2*m+1, 2*m-1] <- -1
	B.tp[2*m] <- min(W)
	B.tp[2*m+1] <- -max(W)

	# ------------------
	# Setup for Min Var
	# ------------------
	# Set up constraints in form Ax \leq B

	# pj non-negative
	A.var <- matrix(0, nrow = 2*m+1, ncol = 2*m-1)
	A.var[1:m-1, 1:m-1] <- -diag(m - 1)
	B.var <- numeric(2 * m - 1)
	# pj sum to less than 1
	A.var[m, 1:m-1] = matrix(1, nrow=1, ncol = m-1)
	B.var[m] = 1
	# thetaj are increasing
	for (i in 1:(m-1)){
		A.var[m+i, (m+i-1):(m+i)] <- c(1, -1)	
	}
	# min(W) < thetaj < max(W)
	A.var[2*m, m] <- -1
	A.var[2*m+1, 2*m-1] <- 1
	B.var[2*m] <- -min(W)
	B.var[2*m+1] <- max(W)

	# ----------------------
	# Perform Minimizations
	# ----------------------
	Diagnostic("Minimizing T(p)")

	tp.min <- Inf
	for (i in 1:n.iter.tp){
		theta0 <- sort(stats::runif(m, min = min(W), max = max(W)))
		p0 <- stats::runif(m, min = 0, max = 1)
		p0 <- p0 / sum(p0)

		x0 <- c(p0[1:(m-1)], theta0)

		min.tp.sol.test <- stats::constrOptim(x0, CalculateTp, NULL, A.tp, B.tp, 
									 		  phi.W = phi.W, weight = weight)
		if (min.tp.sol.test$value < tp.min){
			tp.min <- min.tp.sol.test$value
			min.tp.sol <- min.tp.sol.test
			Diagnostic(tp.min)
		}
	}

	# ConFun(x) <= 0
	ConFun <- function(x){
		Constraints(x, phi.W, weight, min.tp.sol$value)
	}
	
	Diagnostic("Minimizing Variance")

	var.min <- Inf
	for (i in 1:n.iter.var){
		looping <- T
		while (looping){			
			# NlcOptim can't handle it if the Generated QP problem is infeasible
			# However, this error only pops up some of the time. For now, a 
			# dirty solution is to keep running the code with different starting
			# values and ignore everything until we get no error
			looping <- F

			theta0 <- sort(stats::runif(m, min = min(W), max = max(W)))
			p0 <- stats::runif(m, min = 0, max = 1)
			p0 <- p0 / sum(p0)
			x0 <- c(p0[1:(m-1)], theta0)

			min.var.sol.test <- tryCatch({
	    		min.var.sol.test <- NlcOptim::solnl(x0, CalculateVar, ConFun, 
	    											A.var, B.var, 
	    											maxIter = 400, 
													maxnFun = 100*(2*m-1), 
													tolX = 1e-6)
			}, error = function(e){
				err.mess <- paste("ERROR :", conditionMessage(e), sep = " ")
				Diagnostic(err.mess)
				min.var.sol.test <- NULL						
			})

			if (is.null(min.var.sol.test)) {
				looping <- T
			}
		}
		
		if (min.var.sol.test$fn < var.min){
			var.min <- min.var.sol.test$fn
			min.var.sol <- min.var.sol.test
			Diagnostic(var.min)
		}
	}

	Diagnostic(paste("Initial Variance was", CalculateVar(min.tp.sol$par), 
				 	 sep  = " "))
	Diagnostic(paste("Final Variance is", var.min, sep  = " "))

	tp.diff <- tp.min - CalculateTp(min.var.sol$par, phi.W, weight)
	Diagnostic(paste("T(p) decreased by", tp.diff, "while minimizing variance", 
					 sep = " "))
	
	#--------------------------------------------------------------------------#
	# Convert to nice formats and return results
	#--------------------------------------------------------------------------#

	# Convert back to normal vectors
	x.sol <- min.var.sol$par
	p.sol <- c( x.sol[1:m-1], 1 - sum(x.sol[1:m-1]))
	theta.sol <- x.sol[m:(2 * m - 1)]
	simple.sol <- SimplifyPmf(theta.sol, p.sol)
	p.sol <- simple.sol$ProbWeights
	theta.sol <- simple.sol$Support

	x.min.tp <- min.tp.sol$par
	p.min.tp <- c( x.min.tp[1:m-1], 1 - sum(x.min.tp[1:m-1]))
	theta.min.tp <- x.min.tp[m:(2 * m - 1)]

	simple.min.tp <- SimplifyPmf(theta.min.tp, p.min.tp)
	p.min.tp <- simple.min.tp$ProbWeights
	theta.min.tp <- simple.min.tp$Support

	return(list("support" = theta.sol, 
				"probweights" = p.sol, 
				"support.min.tp" = theta.min.tp,
				"probweights.min.tp" = p.min.tp,
				"phi.W" = phi.W,
				"var.opt.results" = min.var.sol,
				"tp.opt.results" = min.tp.sol))
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

CalculateVar <- function(x){
	m <- (length(x) + 1) / 2
	p <- c( x[ 1:(m - 1) ], 1 - sum( x[ 1:(m - 1) ] ) )
	theta <- x[ m:(2 * m - 1) ]
	mean <- sum(p*theta)
	var <- sum(p*(theta - mean)^2)
	return(var)
}

Constraints <- function(x, phi.W, weight, tp.max){
	tp <- CalculateTp(x, phi.W, weight)
	lambda <- 1.0
	const1 <- tp - lambda*tp.max
	return(list(ceq = NULL, c = const1))
}

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

PrintDiagnostic <- function(message, show.diagnostics){
	if (show.diagnostics){
		print(message)
	}
}