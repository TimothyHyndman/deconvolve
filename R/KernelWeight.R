#' @export
# Just exporting for testing, will remove later

KernelWeight <- function(x, weight.type = "Epanechnikov"){
# x is vector of values (symmetric about zero) on which to calculate weight
# weight.type tells us which type of weight we're calcualting

    scale <- x[length(x)]
    weight <- switch( weight.type,
                    Epanechnikov = (0.75/scale)*(1 - (x/scale)^2),
                    Uniform = numeric(length(x)) + 1/(2*scale),
                    Triangular = ( 1 / scale ) * ( 1 - abs( x / scale ) )
                    )

    return(weight)
}