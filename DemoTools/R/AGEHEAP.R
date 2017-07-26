# TODO: Add comment
# 
# Author: tim
###############################################################################



#' calculate Bachi's index of age heaping

#' @description Two Implementations: one following the PASEX spreadsheet SINGAGE (\code{pasex = TRUE}), with ages hard-coded, and another with flexible upper and lower age bounds, but that does not match the PASEX implementation.

#' @param Value a vector of demographic counts by single age
#' @param Age a vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin the lowest age included in calcs. Default 30
#' @param ageMax the upper age bound used for calcs. Default 80
#' @param pasex logical default \code{FALSE}. Do we want to reproduce the specific age weightings in the PASEX spreadsheet?
#' 
#' @details \code{ageMax} is the hard upper bound, treated as interval. If you want ages
#' 20 to 89, then give \code{ageMin = 20} and \code{ageMax = 90}, not 89. These are only heeded if \code{pasex = FALSE}.

#' @export 
#' @examples
#' Value <- c(80626,95823,104315,115813,100796,105086,97266,116328,
#' 		75984,89982,95525,56973,78767,65672,53438,85014,
#' 		47600,64363,42195,42262,73221,30080,34391,29072,
#' 		20531,66171,24029,44227,24128,23599,82088,16454,
#' 		22628,17108,12531,57325,17220,28425,16206,17532,
#' 		65976,11593,15828,13541,8133,44696,11165,18543,
#' 		12614,12041,55798,9324,10772,10453,6773,28358,
#' 		9916,13348,8039,7583,42470,5288,5317,6582,
#' 		3361,17949,3650,5873,3279,3336,27368,1965,
#' 		2495,2319,1335,12022,1401,1668,1360,1185,
#' 		9167,424,568,462,282,6206,343,409,333,291,4137,133,169,157,89,2068,68,81,66,57)
#' Age <- 0:99
#' 
#' Bachi(Value, Age, ageMin = 20, ageMax = 80, pasex = TRUE) # reproduces PASEX SINGAGE
#' Bachi(Value, Age, ageMin = 20, ageMax = 80) # default simpler

Bachi <- function(Value, Age, ageMin = 30, ageMax = 80, pasex = FALSE){
	
	
	# make a matrix for numerators
	w1           <- matrix(0, nrow = length(Age), ncol = 10)
	w2           <- matrix(0, nrow = length(Age), ncol = 10)
	
	if (pasex){
		w1[Age %in% seq(30,70,10),1]  <- 1
		w1[Age %in% seq(31,71,10),2]  <- 1
		w1[Age %in% seq(32,72,10),3]  <- 1
		w1[Age %in% seq(33,63,10),4]  <- 1
		w1[Age %in% c(23,73),4]       <- .5
		w1[Age %in% seq(34,64,10),5]  <- 1
		w1[Age %in% c(24,74),5]       <- .5
		w1[Age %in% seq(35,65,10),6]  <- 1
		w1[Age %in% c(25,75),6]       <- .5
		w1[Age %in% seq(36,66,10),7]  <- 1
		w1[Age %in% c(26,76),7]       <- .5
		w1[Age %in% seq(37,67,10),8]  <- 1
		w1[Age %in% c(27,77),8]       <- .5
		w1[Age %in% seq(28,68,10),9]  <- 1
		w1[Age %in% seq(29,69,10),10] <- 1
		
		# more quirky ranges
		w2[Age > 25 & Age < 75,1]     <- 1
		w2[Age %in% c(25,75),1]       <- .5
		w2[Age > 26 & Age < 76,2]     <- 1
		w2[Age %in% c(26,76),2]       <- .5
		w2[Age > 27 & Age < 77,3]     <- 1
		w2[Age %in% c(27,77),3]       <- .5
		w2[Age > 23 & Age < 73,4]     <- 1
		w2[Age %in% c(23,73),4]       <- .5
		w2[Age > 24 & Age < 74,5]     <- 1
		w2[Age %in% c(24,74),5]       <- .5
		w2[Age > 25 & Age < 75,6]     <- 1
		w2[Age %in% c(25,75),6]       <- .5
		w2[Age > 26 & Age < 76,7]     <- 1
		w2[Age %in% c(26,76),7]       <- .5
		w2[Age > 27 & Age < 77,8]     <- 1
		w2[Age %in% c(27,77),8]       <- .5
		w2[Age > 23 & Age < 73,9]     <- 1
		w2[Age %in% c(23,73),9]       <- .5
		w2[Age > 24 & Age < 74,10]    <- 1
		w2[Age %in% c(24,74),10]      <- .5
	} else {
		markers      <- row(w1) - col(w1)
		w1[markers %% 10 == 0 & markers >= ageMin & markers < ageMax]   <- 1
		w2[markers == ageMin - 5 | markers == ageMax - 5] <- .5
		w2[markers > ageMin - 5 & markers < ageMax - 5]   <- 1
	}
	
	numerators   <- colSums(Value * w1)
	denominators <- colSums(Value * w2)
	
	ratio   <- 100 * numerators / denominators
	ratioeq <- ratio - 10
	sum(abs(ratioeq))  / 2
}

# ----------------------------
# Coale, A. and S. Li (1991) The effect of age misreporting in China on 
# the calculation of mortality rates at very high ages. Demography 28(2)
# ----------------------------

#' Coale-Li age heaping index
#' 
#' @description This implementation is based largely on the sparse verbal description given in
#' Coale and Li (1991): calculate a two-stage 5-term moving average as a reference pattern, then
#' take ratios with respect to this. Ratios for a given terminal digit can then be averaged to produce
#' an index. This procedure was used in that paper for ages 65-100 for mortality rates. 
#' It is probably better suited to rates than counts, but that is not a hard rule.
#' 
#' @param Value a vector of demographic rates or counts by single age
#' @param Age a vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin the lowest age included in calcs. Default 20
#' @param ageMax the upper age bound used for calcs. Default 65
#' @param digit any digit 0-9. Default \code{0}.

#' @details \code{digit} could also be a vector of digits, but the more digits one includes (excepting 0 and 5) the closer the index will get to 1. It is therefore recommended for single digits, or else \code{c(0,5)}
#' 
#' @return the index value 
#' 
#' @references 
#' Coale, A. and S. Li (1991) The effect of age misreporting in China on the calculation of mortality rates at very high ages. Demography 28(2)
#' @export
#' @examples 
#' Value <- c(80626,95823,104315,115813,100796,105086,97266,116328,
#'75984,89982,95525,56973,78767,65672,53438,85014,
#'47600,64363,42195,42262,73221,30080,34391,29072,
#'20531,66171,24029,44227,24128,23599,82088,16454,
#'22628,17108,12531,57325,17220,28425,16206,17532,
#'65976,11593,15828,13541,8133,44696,11165,18543,
#'12614,12041,55798,9324,10772,10453,6773,28358,
#'9916,13348,8039,7583,42470,5288,5317,6582,
#'3361,17949,3650,5873,3279,3336,27368,1965,
#'2495,2319,1335,12022,1401,1668,1360,1185,
#'9167,424,568,462,282,6206,343,409,333,291,4137,133,169,157,89,2068,68,81,66,57)
#' Age <- 0:99
#' CoaleLi(Value, Age, 65, 95, 5, 0) # 3.7
#' CoaleLi(Value, Age, 65, 95, 5, 5) # 3.5 almost just as high
#' Whipple(Value, Age, 65, 95, 0)    # 3.4
#' Whipple(Value, Age, 65, 95, 5)    # 3.2
#' Noumbissi(Value, Age, 65, 95, 0)  # 3.6
#' Noumbissi(Value, Age, 65, 95, 5)  # 3.0

CoaleLi <- function(Value, Age, minAge = 60, maxAge = max(Age), terms = 5, digit = 0){
	
	reference <- ma(ma(Value, n = terms), n = terms)
	
	ratio     <- Value / reference
	
# deprecated: make sure tested age range is divisible by 10.. 
#	agerange  <- maxAge - minAge
#	agerange  <- 10 * floor(agerange / 10)
#	# adjust maxage if necessary
#	maxAge    <- minAge + agerange
	
	ind       <- Age >= minAge & Age <= maxAge
	
	ages      <- max(c(min(Age),minAge)):min(c(maxAge, max(Age)))
	avgRatios <- tapply(ratio[ind], ages %% 10, mean, na.rm = TRUE)
	
	# return avg deviation for specified digit(s)
	mean(avgRatios[as.character(digit)])
}

