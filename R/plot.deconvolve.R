# S3 method for plot
#' @export
plot.deconvolve <- function(decon_obj){
	
	final_plot <- ggplot2::ggplot()

	if ("pdf" %in% names(decon_obj)) {
		x <- decon_obj$x
		y <- decon_obj$pdf
		df_pdf <- data.frame(x,y)
		plot_pdf <- ggplot2::geom_path(data = df_pdf, ggplot2::aes(x,y), 
								   color = "blue", size = 1)
		final_plot <- final_plot + plot_pdf
	}

	if ("probweights" %in% names(decon_obj)) {
		p <- decon_obj$probweights
		theta <- decon_obj$support
		df_pmf_sol <- data.frame(p, theta)
		plot_pmf_sol <- ggplot2::geom_point(data = df_pmf_sol, 
											ggplot2::aes(theta, p), 
											color = "blue")
		final_plot <- final_plot + plot_pmf_sol
	}

	if ("probweights_mintp" %in% names(decon_obj)) {
		p_mintp <- decon_obj$probweights_mintp
		theta_mintp <- decon_obj$support_mintp
		df_mintp <- data.frame(p_mintp, theta_mintp)
		plot_pmf_mintp <- ggplot2::geom_point(data = df_mintp, 
			  							ggplot2::aes(theta_mintp, p_mintp), 
			  							color = "magenta")
		final_plot <- final_plot + plot_pmf_mintp
	}

	final_plot
}