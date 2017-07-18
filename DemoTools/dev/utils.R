
# Author: tim
###############################################################################

#' shift a vector left or right
#' 
#' @description a simple function to move the elements of a vector to earlier or 
#' later positions depending on the value of \code{shift}. Auxiliary to other functions,
#'  rather multipurpose.

#' @param x a vector
#' @param shift an integer value \code{< length(x)}
#' @param fill what value do we put in the newly created positions? Possibly values like \code{FALSE}, \code{NA}, or  \code{0}.

#' @details Nothing fancy here. This is used for example in \code{Noumbissi()} to match denominator ranges to numerator positions using logical vectors.
#' @return the vector x, shifted left or right.
#' @export
shift.vector <- function(x,shift = 0, fill = FALSE){
	n <- length(x)
	
	if (shift > 0){
		x <- c(rep(fill,shift),x[-((n-shift+1):n)])
	}
	if (shift < 0){
		x <- c(x[-(1:abs(shift))],rep(fill,abs(shift)))
	}
	x
}


