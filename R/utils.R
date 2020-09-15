

# Author: tim
###############################################################################

#' Shift a vector left or right.
#'
#' @description Simple function to move the elements of a vector to earlier or
#' later positions depending on the value of \code{shift}. Auxiliary to other functions,
#'  rather multipurpose.

#' @param x vector.
#' @param shift integer. Value \code{< length(x)}. Default 0.
#' @param fill Values to fill newly created positions, e.g. \code{FALSE}, \code{NA}, or  \code{0}.

#' @details Nothing fancy here. This is used for example in \code{Noumbissi()} to match denominator ranges to numerator positions using logical vectors.
#' @return The vector x, shifted left or right.
#' @export
shift.vector <- function(x, shift = 0, fill = FALSE) {
  n          <- length(x)

  if (shift > 0) {
    x        <- c(rep(fill, shift), x[-((n - shift + 1):n)])
  }
  if (shift < 0) {
    x        <- c(x[-(1:abs(shift))], rep(fill, abs(shift)))
  }
  x
}

#' Logging that does not cause jams.
#' @description Zeros happen. Logging zeros leads to errors in downstream code.
#' For the sake of robustness, we replace zeros with a small positive number to be
#' able to continue with calculations.
#' @details Use judiciously, since errors are good too.
#' @param x  numeric or complex vector.
#' @param base positive or complex number. The base with respect to which
#'           logarithms are computed.  Defaults to \code{e=exp(1)}.
#' @export
rlog <- function(x, base = exp(1)) {
  x <- pmax(x, .Machine$double.xmin)
  log(x, base = base)
}

#' A simple centered moving average function.

#' @description This function is defined based on a code chunk found \href{https://stackoverflow.com/questions/743812/calculating-moving-average}{here}.
#' This is a centered moving average of arbitrary width.
#'
#' @param x numeric. Vector to produce a moving average of.
#' @param n integer. The width of the moving average. Default 5 time steps.
#'
#' @return Numeric vector of same length as \code{x}.
#'
#' @details \code{NA} values
#' are used as padding on the left and right to make the returned vector of equal length.
#'
#' @export
#' @examples
#' x  <- runif(100)
#' xx <- cumsum(x)
#' \dontrun{
#' plot(x)
#' lines(ma(x))
#' lines(ma(x),9)
#' }


ma <- function(x, n = 5) {
  as.vector(stats::filter(x, rep(1 / n, n), sides = 2))
}


#' Rescale a vector proportionally to a new sum.

#' @description Function used to rescale a vector to a given value. This is a frequently needed operation.
#'
#' @param x numeric vector.
#' @param scale numeric. Value the vector should sum to. Default 1.
#'
#' @details For a distribution, use \code{scale = 1}. For percentages, use \code{scale = 100}, etc.
#'
#' @return The vector rescaled.
#' @examples
#' x <- runif(10)
#' sum(x)
#' xx <- rescale_vector(x,100)
#' sum(xx)
#' @export

rescale_vector <- function(x, scale = 1) {
   scale * x / sum(x, na.rm = TRUE)
}

#' @title Determine the proportion of a year passed as of a particular date.
#'
#' @description The fraction returned by this is used e.g. for intercensal estimates.
#'
#' @param Year string or integer. 4-digit year.
#' @param Month string or integer. Month digits, 1 or 2 characters.
#' @param Day string or integer. Day of month digits, 1 or 2 characters.
#' @param detect.mid.year logical. If \code{TRUE}, June 30 or July 1 will always return .5.
#' @param detect.start.end logical. Whether or not Jan 1 always be 0 and Dec 31 always be 1. Default \code{TRUE}.
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


ypart           <-
  function(Year,
           Month,
           Day,
           detect.mid.year = TRUE,
           detect.start.end = TRUE) {
    M           <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)

    if (detect.mid.year) {
      .d        <- as.integer(Day)
      .m        <- as.integer(Month)
      if ((.d == 30 & .m == 6) | (.d == 1 & .m == 7)) {
        return(.5)
      }
    }

    if (detect.start.end) {
      .d        <- as.integer(Day)
      .m        <- as.integer(Month)
      if (.d == 1 & .m == 1) {
        return(0)
      }
      if (.d == 31 & .m == 12) {
        return(1)
      }
    }

    monthdur    <- diff(c(M, 365))
    monthdur[2] <- monthdur[2] + lubridate::leap_year(Year)
    M           <- cumsum(monthdur) - 31
    return((M[Month] + Day) / sum(monthdur))

  }


#' Convert date to decimal year fraction.
#'
#' @description Convert a character or date class to decimal, taking into account leap years.
#' @details This makes use of the \code{ypart()} HMD function to compute. If the date is numeric, it is returned as such. If it is \code{"character"}, we try to coerce to \code{"Date"} class, ergo, it is best to specify a character string in an unambiguous \code{"YYYY-MM-DD"} format.  If \code{date} is given in a \code{"Date"} class it is dealt with accordingly.
#'
#' @param date Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}.
#'
#' @return Numeric expression of the date, year plus the fraction of the year passed as of the date.
#' @export
dec.date  <- function(date) {
  if (class(date) == "numeric") {
    return(date)
  }
  if (class(date) == "character") {
    date   <- as.Date(date)
  }
  day 	   <- as.numeric(format(date, '%d'))
  month 	 <- as.numeric(format(date, '%m'))
  year 	   <- as.numeric(format(date, '%Y'))
  frac     <- ypart(
    Year = year,
    Month = month,
    Day = day,
    detect.mid.year = TRUE,
    detect.start.end = TRUE
  )
  year + frac
}


#' Take consecutive ratios of a vector.
#'
#' @description This can be used, for example to take survival ratios.
#' @details Behavior similar to \code{diff()}, in that returned vector is \code{k} elements
#' shorter than the given vector \code{fx}.
#'
#' @param fx numeric. Vector of \code{length > k}.
#' @param k integer. The size of the lag in elements of \code{fx}.
#'
#' @export
#' @examples
#' fx <- 1:10
#' ratx(fx)
#' ratx(fx,-1)
#' ratx(fx,0)
ratx   <- function(fx, k = 1) {
  k    <- as.integer(k)
  N    <- length(fx)
  m    <- N - abs(k)
  if (k > 0) {
    fx <- fx[-c(1:k)] / fx[1:m]
  }
  if (k < 0) {
    fx <- fx[1:m] / fx[-c(1:abs(k))]
  }
  fx
}




#' Wrapper to provide a single location to reference all model life tables.
#' @description Still in the works.
#' @param ModelName string naming the life table to return. Can be "coale-demeny west".
#' @param Sex string indicating which sex should be returned. Can be either "m" or "f".
#' @return list of life tables
#' @details More model families can easily be added.
#' @importFrom demogR cdmltw
#' @export
getModelLifeTable <- function(ModelName, Sex) {
  Sex             <- toupper(substr(Sex, 1, 1))

  stopifnot(Sex %in% c("M", "F"))
  stopifnot(ModelName == "coale-demeny west")

  outputLT        <- demogR::cdmltw(sex = Sex)

  return(outputLT)
}

#' calculate average of vector elements adjacent to and excluding the index element
#' @description Calculate average of vector elements adjacent to and excluding the index element. For example, the second element of the result is the average of the first and third elements of the input vector \code{x}. Used by \code{smooth_age_5_zigzag_inner()}, and possibly useful elsewhere.
#' @details Tails are given a value of \code{NA}.
#' @param x numeric vector
#' @return numeric vector the same length as \code{x}.
#' @export
#' @examples
#' x <- 1:10
#' all(avg_adj(x) == x, na.rm = TRUE)

avg_adj <- function(x) {
  n <- length(x)
  (shift.vector(x,-1, NA) + shift.vector(x, 1, NA)) / 2
}

#' convert strings to concatenation of lower case alphabet
#' @description all characters are converteusd to lower case, and all non-alphabet characters are removed. This is useful for reducing checks on names of user-specified strings, like \code{method} arguments. For example, \code{"Carrier-Farrag"}, \code{"Carrier Farrag"}, or \code{"carrier_farrag"} all get converted to \code{"carrierfarrag"}.
#' @param string character string
#' @return character string
#' @export
#' @examples
#' simplify.text(c("carrier_farrag","CarrierFarrag","Carrier-Farrag"))

simplify.text <- function(string) {
  lc <- tolower(string)
  #	stringr::str_replace_all(lc, "[^[:lower:]]", "")
  gsub("[^a-z]", "" , lc)
}

# deprecated functions

# TR: deprecated 20 July, 2018. Use graduate_uniform() instead
##' Convert arbitrary age groupings into single years of age
##'
##' @description Splits aggregate counts for a vector of age groups of a common width into single year age groups.
##'
##' @param Value numeric vector of age group counts
##' @param Age numeric vector of ages corresponding to the lower integer bound of the age range.
##' @param OAG boolean argument that determines whether the final age group (assumed open ended) is kept as it is or has the same length as the rest of the age groups. Default is FALSE, i.e. use the same length for the final age group.
##'
##' @return numeric vector of counts for single year age groups.
##'
##' @details Assumes that all age groups are equal width. (Default is 5-year age groups.) Also, assumes that the population is uniformly distributed across each age interval. If there is a split under 5 age group (0 and 1-4 age groups), use \code{groupAges()} to consolidate before using this function. The default setting is to assume that the final age group is the same width as the other age groups.
##'
##' @export
##' @examples
# MalePop <- c(9544406,7471790,11590109,11881844,11872503,12968350,
# 11993151,10033918,14312222,8111523,15311047,6861510,13305117,7454575,
#              9015381,10325432,9055588,5519173)
# Ages <- seq(0,85, by=5)
##' splitToSingleAges(MalePop, Ages)
##' splitToSingleAges(MalePop, Ages, OAG = TRUE)
#splitToSingleAges <- function(Value, Age, OAG = FALSE){
#	ageMax <- max(Age)             # the lower bound of the largest age group
#	N      <- length(Value)        # number of age groups from the census.
#	M      <- ageMax / (N - 1)     # length of each age group from the census.
#
#	ageGroupBirths   <- Value / M  # vector of the number of births in a single year for each age group assuming uniformity.
#
#	if (OAG){
#		ageGroupBirths <- ageGroupBirths[0:(length(ageGroupBirths)-1)] # remove final age group
#		singleAgeGroupBirths <- rep(ageGroupBirths, each = M) # vector of the single year births for all ages except final age group
#		singleAgeGroupBirths[length(singleAgeGroupBirths)+1] <- Value[length(Value)]
#	}
#	else{
#		singleAgeGroupBirths  <- rep(ageGroupBirths, each = M)  # vector of the single year births for all ages
#	}
#
#	return(singleAgeGroupBirths)
#}
