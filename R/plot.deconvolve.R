# S3 method for plot
#' @export
plot.deconvolve <- function(x, ...){
	
	final_plot <- ggplot2::ggplot()
	text_col = "darkslategrey"
	dens_col = "black"
	# Plot distribution --------------------------------------------------------
	if ("pdf" %in% names(x)) {
		xvec <- x$x
		yvec <- x$pdf
		df_pdf <- data.frame(xvec,yvec)
		plot_pdf <- ggplot2::geom_path(data = df_pdf, ggplot2::aes(xvec,yvec), 
								   color = dens_col, size = 1)
		final_plot <- final_plot + plot_pdf
	}

	if ("probweights" %in% names(x)) {
		p <- x$probweights
		theta <- x$support
		df_pmf_sol <- data.frame(p, theta)
		plot_pmf_sol <- ggplot2::geom_point(data = df_pmf_sol, 
											ggplot2::aes(theta, p), 
											color = dens_col)
		final_plot <- final_plot + plot_pmf_sol
	}

	if ("probweights_mintp" %in% names(x)) {
		p_mintp <- x$probweights_mintp
		theta_mintp <- x$support_mintp
		df_mintp <- data.frame(p_mintp, theta_mintp)
		plot_pmf_mintp <- ggplot2::geom_point(data = df_mintp, 
			  							ggplot2::aes(theta_mintp, p_mintp), 
			  							color = "magenta")
		final_plot <- final_plot + plot_pmf_mintp
	}

	# Add title ----------------------------------------------------------------
	if ("pdf" %in% names(x)){
		title <- 'Deconvolved Density'
	} else {
		title <- 'Deconvolved Distribution'
	}
	final_plot <- final_plot + ggplot2::ggtitle(title) + 
				  ggplot2::theme(plot.title = ggplot2::element_text(size=20, face="bold",
				  				 hjust = 0.5, color = text_col))
	
	# X and Y labels
	final_plot <- final_plot + ggplot2::labs(x = "x", y = "f(x)") + 
				  ggplot2::theme(
				  	axis.title.x = ggplot2::element_text(color = text_col, vjust=-0.35), 
				  	axis.title.y = ggplot2::element_text(color = text_col, vjust=0.35))

	# Add legend ---------------------------------------------------------------
	

	final_plot
}
