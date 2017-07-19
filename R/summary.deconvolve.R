# S3 summary method for objects of class deconvolve
#' @export

summary.deconvolve <- function(decon_obj){
	cat("Object of class ", class(decon_obj), ".\n", sep = "")

	obj_names <- names(decon_obj)

	for (name in obj_names){
		cat('\n', name, ':\n', sep = "")

		if (name == 'support' | name == 'probweights'){
			print(decon_obj[[name]])
		} else if (name == 'x'){
			x <- decon_obj[[name]]
			Min. <- min(x)
			Max. <- max(x)
			if (!any(!(x == seq(Min., Max., length.out = length(x))))){
				Deltax. = x[2] - x[1]
				print(data.frame(Min., Deltax., Max.), row.names = FALSE)
			} else {
				print(data.frame(Min., Max.), row.names = FALSE)	
			}
		} else {
			print(summary(decon_obj[[name]]))	
		}
	}
}