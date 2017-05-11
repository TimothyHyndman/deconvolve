
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


	# Minimize the variance subject to T(theta,p) not increasing
}
