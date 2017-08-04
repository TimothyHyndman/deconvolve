# S3 summary method for objects of class deconvolve
#' @export

summary.deconvolve <- function(object, ...){
	cat("Object of class ", class(object), ".\n", sep = "")

	obj_names <- names(object)

	for (name in obj_names){
		cat('\n', name, ':\n', sep = "")

		if (name == 'support' | name == 'probweights' | name == 'pdf'){
			print(object[[name]])
		} else if (name == 'x'){
			x <- object[[name]]
			Min. <- min(x)
			Max. <- max(x)
			if (!any(!(x == seq(Min., Max., length.out = length(x))))){
				Deltax. = x[2] - x[1]
				print(data.frame(Min., Deltax., Max.), row.names = FALSE)
			} else {
				print(data.frame(Min., Max.), row.names = FALSE)	
			}
		} else {
			print(summary(object[[name]]))	
		}
	}
}