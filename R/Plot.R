#' @export
PlotPmf <- function(pmf){
	# Plot for diagnostics
	p <- pmf$probweights
	theta <- pmf$support

	p.min.tp <- pmf$probweights.min.tp
	theta.min.tp <- pmf$support.min.tp

	df <- data.frame(p, theta)
	df.min.tp <- data.frame(p.min.tp, theta.min.tp)
	plot <- ggplot2::ggplot() + 
			ggplot2::geom_point(data = df.min.tp, 
			  					ggplot2::aes(theta.min.tp, p.min.tp), 
			  					color = "magenta") + 
			ggplot2::geom_point(data = df, 
			  					ggplot2::aes(theta, p), 
			  					color = "blue")	

	plot
}

#' @export
PlotPdf <- function(pdf){

	x <- pdf$x
	y <- pdf$y

	df <- data.frame(x,y)
	
	ggplot2::ggplot(data = df, ggplot2::aes(x,y)) + 
			ggplot2::geom_path(color = "blue", size = 1)
			# ggplot2::geom_path() + 
			# ggplot2::geom_point()
}

#' @export
PlotPmfAndPdf <- function(pmf, pdf){
	p <- pmf$probweights
	theta <- pmf$support

	p.min.tp <- pmf$probweights.min.tp
	theta.min.tp <- pmf$support.min.tp

	df.sol <- data.frame(p, theta)
	df.min.tp <- data.frame(p.min.tp, theta.min.tp)

	x <- pdf$x
	y <- pdf$y
	df <- data.frame(x,y)
	
	ggplot2::ggplot(data = df, ggplot2::aes(x,y)) + 
			ggplot2::geom_path(color = "blue", size = 1) +
			ggplot2::geom_point(data = df.min.tp, 
			  					ggplot2::aes(theta.min.tp, p.min.tp), 
			  					color = "magenta") + 
			ggplot2::geom_point(data = df.sol, 
			  					ggplot2::aes(theta, p), 
			  					color = "blue")	
}