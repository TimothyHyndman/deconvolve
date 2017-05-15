#' @export
# Just exporting for testing, will remove later

ComputePhi <- function(W, tt, longt){
# W is vector of data for which we calculate PhiW
# tt is coarse vector on which we initially calculate PhiW
# longt is length of fine vector we want to calculate PhiW on
#
# Returns a list containing Real, Imaginary, Normed parts of Phi and the T 
# values on which they are calculated

	n <- length( W )
	
	# Estimate empirical characteristic function of W on tt
	oo <- outerop( tt, W, '*' )
	re.hat.phi <- rowSums( cos( oo ) ) / n
	im.hat.phi <- rowSums( sin( oo ) ) / n
	norm.hat.phi <- sqrt( re.hat.phi^2 + im.hat.phi^2 )

	# Calculate t^* where t^* is smallest t>0 such that |hat.phi(t^*)| <= n^(-0.25)
	tmp <- tt[norm.hat.phi < n^(-0.25)]

	if ( length( tmp[ tmp > 0] ) == 0 ){
		t.star <- max( tt )
	} else {
		t.star <- min( tmp[ tmp > 0 ] )
	}

	# Estimate empirical characteristic function of W on [-t.star, t.star]
	tt.new <- seq( -t.star, t.star, length.out = longt )
	
	oo <- outerop( tt.new, W, "*" )
	re.hat.phi <- rowSums( cos( oo ) ) / n
	im.hat.phi <- rowSums( sin( oo ) ) / n
	norm.hat.phi <- sqrt( re.hat.phi^2 + im.hat.phi^2 )

	return( list( "norm" = norm.hat.phi, "im" = im.hat.phi, "re" = re.hat.phi,
				  "t.values" = tt.new) )
}