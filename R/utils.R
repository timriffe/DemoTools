

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
#' # some de facto unit tests:
#' x <- 1:10
#' d2 <- ma(x,2) - x
#' d3 <- ma(x,3) - x
#' d4 <- ma(x,4) - x
#' d5 <- ma(x,5) - x
#' # pop <- sample.int(10, 5, replace = T)
#' # all should give same for linear data.
#' stopifnot(all(abs(d2) < 1e-10, na.rm = TRUE))
#' stopifnot(all(abs(d3) < 1e-10, na.rm = TRUE))
#' stopifnot(all(abs(d4) < 1e-10, na.rm = TRUE))
#' stopifnot(all(abs(d5) < 1e-10, na.rm = TRUE))

ma <- function(x, n = 5) {
  # This is PJ's fix from issue #127
  if (n %% 2 == 1) {   # odd as before
    out <- as.vector(stats::filter(x, rep(1 / n, n), sides = 2))
  } else { # even...
    temp <- as.vector(stats::filter(x, rep(1 / (2 * n), n), sides = 1))
    out  <- shift.vector(temp, 
                         shift = -n / 2 + 1, 
                         fill = NA) + 
            shift.vector(temp, 
                         shift = -n / 2, fill = NA)
  }
  out
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

#' Convert date to decimal year fraction.
#'
#' @description Convert a character or date class to decimal, taking into account leap years.
#' @details This makes use of the \code{lubridate::decimal_date} to compute the proportion of the year that has passed. If the date is numeric, it is returned as such. If it is \code{"character"}, we try to coerce to date through \code{lubridate::ymd}, ergo, it is best to specify a character string in an unambiguous \code{"YYYY-MM-DD"} format.  If \code{date} is given in a \code{"Date"} class it is dealt with accordingly.
#'
#' @param date Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}.
#'
#' @return Numeric expression of the date, year plus the fraction of the year passed as of the date.
#'
#' @export
dec.date  <- function(date) {

  if (inherits(date, "numeric")) {
    return(date)
  }

  res <- lubridate::decimal_date(lubridate::ymd(date))
  res
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
  out <- (shift.vector(x,-1, NA) + shift.vector(x, 1, NA)) / 2
  names(out) <- names(x)
  out
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

#' Convert single age groups to five-year group abridged
#' @description Convert single age groups to five-year group abridged
#' @param Age a numeric vector with counts for age groups. These are assumed to begin at zero.
#' @return A numeric vector with five-year age groups abridged
#'
#' @export
#' @examples
#' pop_male_counts <-
#'   c(
#'     `0` = 11684,
#'     `1` = 11473,
#'     `2` = 11647,
#'     `3` = 11939,
#'     `4` = 11680,
#'     `5` = 10600,
#'     `6` = 11100,
#'     `7` = 11157,
#'     `8` = 11238,
#'     `9` = 11544,
#'     `10` = 7216,
#'     `11` = 7407,
#'     `12` = 7461,
#'     `13` = 7656,
#'     `14` = 7774,
#'     `15` = 5709,
#'     `16` = 5629,
#'     `17` = 5745,
#'     `18` = 6056,
#'     `19` = 6259,
#'     `20` = 5303,
#'     `21` = 5423,
#'     `22` = 5497,
#'     `23` = 5547,
#'     `24` = 5417,
#'     `25` = 5441,
#'     `26` = 5466,
#'     `27` = 5500,
#'     `28` = 5668,
#'     `29` = 5694,
#'     `30` = 4365,
#'     `31` = 4252,
#'     `32` = 4122,
#'     `33` = 4142,
#'     `34` = 4039,
#'     `35` = 3210,
#'     `36` = 3222,
#'     `37` = 3258,
#'     `38` = 3413,
#'     `39` = 3871,
#'     `40` = 2684,
#'     `41` = 2844,
#'     `42` = 3052,
#'     `43` = 3182,
#'     `44` = 3237,
#'     `45` = 2263,
#'     `46` = 2298,
#'     `47` = 2318,
#'     `48` = 2257,
#'     `49` = 2194,
#'     `50` = 2231,
#'     `51` = 2172,
#'     `52` = 2072,
#'     `53` = 2008,
#'     `54` = 1932,
#'     `55` = 1301,
#'     `56` = 1262,
#'     `57` = 1213,
#'     `58` = 1197,
#'     `59` = 1191,
#'     `60` = 1601,
#'     `61` = 1593,
#'     `62` = 1490,
#'     `63` = 1348,
#'     `64` = 1299,
#'     `65` = 568,
#'     `66` = 745,
#'     `67` = 843,
#'     `68` = 801,
#'     `69` = 925,
#'     `70` = 806,
#'     `71` = 883,
#'     `72` = 796,
#'     `73` = 725,
#'     `74` = 672,
#'     `75` = 470,
#'     `76` = 441,
#'     `77` = 340,
#'     `78` = 300,
#'     `79` = 289,
#'     `80` = 4200
#'   )
#' pop_female_counts <-
#'   c(
#'     `0` = 11673,
#'     `1` = 11474,
#'     `2` = 11670,
#'     `3` = 11934,
#'     `4` = 11614,
#'     `5` = 10603,
#'     `6` = 11144,
#'     `7` = 11179,
#'     `8` = 11269,
#'     `9` = 11617,
#'     `10` = 6772,
#'     `11` = 6948,
#'     `12` = 7030,
#'     `13` = 7211,
#'     `14` = 7306,
#'     `15` = 6531,
#'     `16` = 6443,
#'     `17` = 6535,
#'     `18` = 6951,
#'     `19` = 7213,
#'     `20` = 6096,
#'     `21` = 6234,
#'     `22` = 6327,
#'     `23` = 6410,
#'     `24` = 6285,
#'     `25` = 6464,
#'     `26` = 6492,
#'     `27` = 6549,
#'     `28` = 6739,
#'     `29` = 6795,
#'     `30` = 5013,
#'     `31` = 4888,
#'     `32` = 4735,
#'     `33` = 4747,
#'     `34` = 4646,
#'     `35` = 3040,
#'     `36` = 3068,
#'     `37` = 3107,
#'     `38` = 3246,
#'     `39` = 3658,
#'     `40` = 2650,
#'     `41` = 2788,
#'     `42` = 2977,
#'     `43` = 3108,
#'     `44` = 3156,
#'     `45` = 1756,
#'     `46` = 1784,
#'     `47` = 1802,
#'     `48` = 1764,
#'     `49` = 1724,
#'     `50` = 1982,
#'     `51` = 1935,
#'     `52` = 1846,
#'     `53` = 1795,
#'     `54` = 1731,
#'     `55` = 863,
#'     `56` = 850,
#'     `57` = 825,
#'     `58` = 819,
#'     `59` = 816,
#'     `60` = 1348,
#'     `61` = 1342,
#'     `62` = 1246,
#'     `63` = 1138,
#'     `64` = 1101,
#'     `65` = 391,
#'     `66` = 520,
#'     `67` = 585,
#'     `68` = 560,
#'     `69` = 659,
#'     `70` = 670,
#'     `71` = 750,
#'     `72` = 686,
#'     `73` = 634,
#'     `74` = 604,
#'     `75` = 353,
#'     `76` = 340,
#'     `77` = 270,
#'     `78` = 246,
#'     `79` = 247,
#'     `80` = 4143
#'   )
#'
#' single2abridged(pop_male_counts)
#' 
single2abridged <- function(Age) {
  age_names <- as.integer(seq_along(Age) - 1)
  abridged_index <- calcAgeAbr(age_names)
  tapply(Age, abridged_index, sum)
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
