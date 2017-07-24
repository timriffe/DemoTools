
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
	as.vector(stats::filter(x, rep(1 / n, n), sides = 2))
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

#' @title determine whether a year is a leap year. 
#' 
#' @description In order to remove \code{lubridate} dependency, we self-detect leap years and adjust February accordingly. Code inherited from HMD.
#' 
#' @param Year integer of year to query
#' 
#' @return logical is the Year a leap year or not
#' 
#' @export
#' @author Carl Boe

isLeapYear <- function (Year){      # CB: mostly good algorithm from wikipedia
	ifelse(
			( (Year %% 4) == 0  &  (Year %% 100) != 0   ) | ( (Year %% 400) == 0 ),
			TRUE, FALSE )
}

#' @title determine the proportion of a year passed as of a particular date
#' 
#' @description The fraction returned by this is used e.g. for intercensal estimates.
#' 
#' @param Year 4-digit year (string or integer)
#' @param Month month digits (string or integer, 1 or 2 characters)
#' @param Day Day of month digits (string or integer, 1 or 2 characters)
#' @param detect.mid.year logical. if \code{TRUE}, June 30 or July 1 will always return .5.
#' @param detect.start.end logical. default \code{TRUE}. Should Jan 1 always be 0 and Dec 31 always be 1?
#' @details Code inherited from HMD, slightly modified to remove matlab legacy bit.
#' @export
#' @examples
#' ypart(2001,2,14) # general use
#' ypart(2001,6,30) # mid year options, default detection
#' ypart(2001,7,1)  # also
#' ypart(2000,6,30) # no change for leap year, still detected as mid year
#' ypart(2000,6,30, FALSE) # exact measure
#' ypart(2000,7,1, FALSE)  # July 1 leap year
#' ypart(2001,7,1, FALSE)  # July 1 not leap year
#' ypart(2002,12,31, detect.start.end = FALSE) # assumes end of day by default.
#' ypart(2002,1,1, detect.start.end = FALSE) # end of day year fraction
#' ypart(2002,1,1, detect.start.end = TRUE)  # assume very begining of year


ypart <- function(Year, Month, Day, detect.mid.year = TRUE, detect.start.end = TRUE){
	M <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)

	if (detect.mid.year){
		.d <- as.integer(Day)
		.m <- as.integer(Month)
		if ((.d == 30 & .m == 6) | (.d == 1 & .m == 7)){
			return(.5)
		}
	}
	
	if (detect.start.end){
		.d <- as.integer(Day)
		.m <- as.integer(Month)
		if (.d == 1 & .m == 1){
			return(0)
		}
		if(.d == 31 & .m == 12){
			return(1)
		}
	}
	
	monthdur    <- diff(c(M,365))
	monthdur[2] <- monthdur[2] + isLeapYear(Year)
	M           <- cumsum(monthdur) - 31
	return((M[Month] + Day) / sum(monthdur))
	
}


#' a convert date to decimal year fraction
#'  
#' @description Convert a character or date class to decimal, taking into account leap years. 

#' @details This makes use of two HMD functions, \code{ypart()}, and \code{isLeapYear()} to compute. If the date is numeric, it is returned as such. If it is \code{"character"}, we try to coerce to \code{"Date"} class, ergo, it is best to specify a character string in an unambiguous \code{"YYYY-MM-DD"} format. If \code{date} is given in a \code{"Date"} class it is dealt with accordingly.
#' @param date either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}.
#' 
#' @return numeric expression of the date, year plus the fraction of the year passed as of the date.
#' @export
dec.date <- function(date){
	if (class(date) == "numeric"){
		return(date)
	}
	if (class(date) == "character"){
		date <- as.Date(date)
	}
	day 	<- as.numeric(format(date,'%d'))
	month 	<- as.numeric(format(date,'%m'))
	year 	<- as.numeric(format(date,'%Y'))
	frac    <- ypart(Year = year, 
		                  Month = month, 
		                  Day = day,
		                  detect.mid.year = TRUE,
						  detect.start.end = TRUE)
	year + frac
}