
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


#' group single ages into equal age groups of arbitrary width
#' 
#' @description This can be useful to check constrained sums, or as an intermediate step for smoothing.
#' 
#' @param Value numeric vector of single age counts
#' @param Age integer vector of lower bounds of single age groups
#' @param N integer. Default 5. The desired width of resulting age groups
#' @param shiftdown integer. Default 0. Optionally shift age groupings down by single ages 
#' 
#' @return vector of counts in N-year age groups
#' 
#' @details If you shift the groupings, then the first age groups may have a negative lower bound 
#' (for example of -5). These counts would be discarded for the osculatory version of Sprague smoothing,
#' for example, but they are preserved in this function. The important thing to know is that if you shift 
#'  the groups, the first and last groups won't be N years wide. For example if \code{shiftdown} is 1, the first age group is 4-ages wide. The ages themselves are not returned, 
#' but they are the name attribute of the output count vector.

#' @examples
#' India1991males <- c(9544406,7471790,11590109,11881844,11872503,12968350,11993151,10033918
#' 		,14312222,8111523,15311047,6861510,13305117,7454575,9015381,10325432,9055588,5519173,12546779,4784102,13365429,4630254,9595545,4727963
#' 		,5195032,15061479,5467392,4011539,8033850,1972327,17396266,1647397,6539557,2233521,2101024,16768198,3211834,1923169,4472854,1182245
#' 		,15874081,1017752,3673865,1247304,1029243,12619050,1499847,1250321,2862148,723195,12396632,733501,2186678,777379,810700,7298270
#' 		,1116032,650402,1465209,411834,9478824,429296,1190060,446290,362767,4998209,388753,334629,593906,178133,4560342,179460
#' 		,481230,159087,155831,1606147,166763,93569,182238,53567,1715697,127486,150782,52332,48664,456387,46978,34448
#' 		,44015,19172,329149,48004,28574,9200,7003,75195,13140,5889,18915,21221,72373)
#' Age <- 0:100
#' groupAges(India1991males, N = 5)
#' groupAges(India1991males, N = 5, shift = 1)
#' groupAges(India1991males, N = 5, shift = 2)
#' groupAges(India1991males, N = 5, shift = 3)
#' groupAges(India1991males, N = 5, shift = 4)
groupAges <- function(Value, Age = 1:length(Value) - 1, N = 5, shiftdown = 0){
	shift <- abs(shift)
	stopifnot(shift < N)
	
	ageN  <- (Age + shift) - (Age + shift) %% N 
	tapply(Value, ageN, sum)
}