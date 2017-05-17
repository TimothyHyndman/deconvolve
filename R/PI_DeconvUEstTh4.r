PI_DeconvUEstTh4 <- function(W, t.limits, phi.U, hat.var.U, phi.U, tt){
	n <- length(W)
	PhiK <- function(t){
		(1 - t^2)^3
	}
	mu.K.2 <- 6
	

	st.dev.X <- max( !is.na(sqrt( stats::var(W) - hat.var.U )), 1/n )
	th4 <- (st.dev.X^(-9) * factorial(8) ) / (2^9 * factorial(4) * sqrt(pi))

	# Grid of h values on which to search for a solution
	max.h <- ( max(W) - min(W) ) / 10
	h.naive <- 1.06 * sqrt( stats::var(W) ) * n^(-1/5)
	h.grid <- seq(h.naive / 3, max.h, length.out = 100)
	# length.h.grid <- length(h.grid)

	dt <- 0.0002
	t <- seq(-1, 1, by = dt)

	t <- matrix(t, ncol = 1)
	h.grid <- matrix(h.grid, nrow = 1)
	t.over.h <- t %*% (1 / h.grid)

	phi.U.2 <- PhiUSpline(t.over.h, hat.var.U, t.limits, phi.U, tt)
	phi.K.2 <- (PhiK(t))^2

	# Find h3 for th3
	rr <- 3
	term1 <- -h.grid^2 * mu.K.2 * th4
	term2 <- matrix( rep( t^(2*rr) * phi.K.2, length(h.grid) ), 
					ncol = length(h.grid) ) 
	term2 <- term2 / phi.U.2
	term2 <- colSums(term2) * dt / ( 2 * pi * n * h.grid^(2 * rr + 1) )

	A.bias.2 <- (term1 + term2)^2
	ind.h3 <- which(A.bias.2 == min(A.bias.2))
	ind.h3 <- ind.h3[1]
	h3 <- h.grid[ind.h3]

	# Estimate empirical characterstic function of W
	phi.W <- ComputePhiEmp(W, t/h3)
	th3 <- t^(2 * rr) * phi.W$norm * phi.K.2 / phi.U.2[,ind.h3]
	th3 <- sum(th3) * dt / ( 2 * pi * h3^(2 * rr + 1) )

	# Find h2 for th2
	rr <- 2
	

}

PhiUSpline <- function(t.over.h, hat.var.U, t.limits, phi.U, tt){
	ind1 <- ( t.over.h >= t.limits[1] ) & ( t.over.h <= t.limits[2] )
	ind2 <- ( t.over.h <  t.limits[1] ) | ( t.over.h >  t.limits[2] )

	PhiULap <- function(t){
		1 / ( 1 + hat.var.U / 2 * t^2)
	}

	y <- 0*t.over.h

	y[ind1] <- spline( tt, phi.U, xout = t.over.h[ind1] )$yinch #What method is this using? Compare to MATLAB code.
	y[ind2] <- PhiULap(t.over.h[ind2])

	return(y)
}