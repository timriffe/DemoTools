
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



#' a simple centered moving average function

#' @description This function is defined based on a code chunk found here \url{https://stackoverflow.com/questions/743812/calculating-moving-average}. This is a centered moving average of arbitrary width. 
#' 
#' @param x numeric vector to produce a moving average of
#' @param n integer the width of the moving avarage. Default 5 time steps (years)
#' 
#' @return numeric vector of same length as \code{x}.
#' 
#' @details \code{NA} values
#' are used as padding on the left and right to make the returned vector of equal length.
#' 
#' @export
#' @examples 
#' x  <- runif(100)
#' xx <- cumsum(x)
#' plot(x)
#' lines(ma(x))
#' lines(ma(x),9)


ma <- function(x,n=5){
	stats::filter(x, rep(1 / n, n), sides = 2)
}


#' rescale a vector proportionally to a new sum

#' @description THis is a frequently needed operation, so it's being added as a general utility.
#' 
#' @param x numeric vector
#' @param scale what you want the vector to sum to. Default 1.
#' 
#' @details For a distribution, use \code{scale = 1}. For percentages, use \code{scale = 100}, etc.
#' 
#' @return the vector, rescaled
#' @examples 
#' x <- runif(10)
#' sum(x)
#' xx <- rescale.vector(x,100)
#' sum(xx)
#' @export

rescale.vector <- function(x, scale = 1){
	scale * x / sum(x, na.rm = TRUE)
}