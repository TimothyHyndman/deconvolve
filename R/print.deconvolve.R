# S3 print method for objects of class deconvolve
#' @export

print.deconvolve <- function(x, ...){
	cat("Object of class ", class(x), ".\n \n", sep = "")

	if (is.null(x$pdf)){
		cat("$support \n")
		print(x$support)
		cat("$probweights \n")
		print(x$probweights)
	} else {
		cat("$pdf \n")
		print(x$pdf)
	}
}