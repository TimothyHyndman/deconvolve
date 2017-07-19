# S3 print method for objects of class deconvolve
#' @export

print.deconvolve <- function(decon_obj){
	cat("Object of class ", class(decon_obj), ".\n \n", sep = "")
	print(unclass(decon_obj))

	# obj_names <- names(decon_obj)

	# for (name in obj_names){
	# 	cat('$', name, '\n', sep = "")

	# 	if (length(yy[[name]]) > 6){
	# 		cat(head(decon_obj[[name]]), '...', '\n')
	# 	} else {
	# 		print(decon_obj[[name]])
	# 	}
		
	# 	cat('\n')
	# }
}