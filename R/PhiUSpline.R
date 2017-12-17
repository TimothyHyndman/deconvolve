PhiUSpline <- function(t.over.h, hat.var.U, phi.U, tt){
	t.limits <- c(min(tt), max(tt))
	ind1 <- ( t.over.h >= t.limits[1] ) & ( t.over.h <= t.limits[2] )
	ind2 <- ( t.over.h <  t.limits[1] ) | ( t.over.h >  t.limits[2] )

	PhiULap <- function(t){
		1 / ( 1 + hat.var.U / 2 * t^2)
	}

	y <- 0*t.over.h

	y[ind1] <- stats::spline( tt, phi.U, xout = t.over.h[ind1] )$y 
	#What method is this using? Compare to MATLAB code.
	y[ind2] <- PhiULap(t.over.h[ind2])

	return(y)
}
