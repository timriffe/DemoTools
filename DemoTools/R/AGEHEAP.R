
# Author: tim
###############################################################################


#' calculate Whipple's index of age heaping

#' @description Implementation following the PASEX spreadsheet SINGAGE, with some extra options for more digit checking.

#' @param Value a vector of demographic counts by single age
#' @param Age a vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin the lowest age included in calcs. Default 25
#' @param ageMax the upper age bound used for calcs. Default 65
#' @param digit any digit 0-9. Default \code{c(0,5)}. Otherwise it needs to be a single digit.
#'  
#' 
#' @details \code{ageMax} is the hard upper bound, treated as interval. If you want ages
#' 20 to 89, then give \code{ageMin = 20} and \code{ageMax = 90}, not 89. You can get 
#' arbitrary W(i) indicies by specifying other digits. Note you can only do pairs of digits
#' if it's 0,5. Otherwise just one digit at a time. 

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
#' Whipple(Value, Age, 25, 60, digit = c(0,5)) # 2.34,replicates SINGAGE males
#' 
#' # implements formula from GDA_1981_Structures-par-age-et-sexe-en-Afrique_[IREDA]
#' # p 148
#' Whipple(Value, Age, 25, 60, digit = 0)
#' Whipple(Value, Age, 25, 60, digit = 5) 
#' 
#' # Whipple types
#' Whipple(Value, Age, 25, 60, digit = 3) 

Whipple <- function(Value, Age, ageMin = 25, ageMax = 65, digit = c(0,5)){
	stopifnot(length(digit) <= 2)
	stopifnot(length(Value) == length(Age))
	
	numeratorind   <- Age >= ageMin & Age <= ageMax & Age %% 10 %in% digit
	# if we are checking just one digit, go down 7 up two, so that the right nr
	# of counts in denom. This per the French formulas.
	denominatorind <- Age >= (ageMin - ifelse(all(digit %in% c(0,5)),2,7)) & Age <= (ageMax + 2)
	
	whip           <- ifelse(length(digit) == 2,5,10) * sum(Value[numeratorind]) / sum(Value[denominatorind])
	
	return(whip)
}


#' calculate Myer's blended index of age heaping

#' @description Implementation following the PASEX spreadsheet SINGAGE

#' @param Value a vector of demographic counts by single age
#' @param Age a vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin the lowest age included in calcs. Default 10
#' @param ageMax the upper age bound used for calcs. Default 90

#' @details \code{ageMax} is the hard upper bound, treated as interval. If you want ages
#' 20 to 89, then give \code{ageMin = 20} and \code{ageMax = 90}, not 89.

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
#' Myers(Value, Age, 10, 90) * 2 #47.46, replicates SINGAGE males

Myers <- function(Value, Age, ageMin = 10, ageMax = 90){
	
	# hard code period to 10 for digits
	period <- 10
	
	# must be of same length for indexing
	stopifnot(length(Value) == length(Age))
	stopifnot(ageMin %% period == 0 & ageMax %% period == 0)
	
	
	# select out ages, place into matrix for summing over digits
	ind     <- Age >= ageMin & Age < ageMax
	# a row corresponds to a digit
	VA      <- matrix(Value[ind], nrow = period, dimnames = list(0:(period-1), NULL))
	
	# sum staggered, once without the youngest group but with the oldest one (tab2)
	# and once with the youngest and without the oldest
	tab1    <- rowSums(VA) # differs from other implementations, but matches PASEX
	
	tab2    <- rowSums(VA[, - 1])
	
	# weighted tabulation
	TAB     <- tab1 * 1:period + tab2 * c(period:1 - 1)
	# interpret as % that would need to be redistributed...
	my      <- sum(abs(TAB / sum(TAB) - 1 / period)) * 50
	return(my)
}


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
#' @param terms integer (default 5). How wide shall the (centered) moving average be?
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

CoaleLi <- function(Value, Age, ageMin = 60, ageMax = max(Age), terms = 5, digit = 0){
	
	reference <- ma(ma(Value, n = terms), n = terms)
	
	ratio     <- Value / reference
	
# deprecated: make sure tested age range is divisible by 10.. 
#	agerange  <- maxAge - minAge
#	agerange  <- 10 * floor(agerange / 10)
#	# adjust maxage if necessary
#	maxAge    <- minAge + agerange
	
	ind       <- Age >= ageMin & Age <= ageMax
	
	ages      <- max(c(min(Age),ageMin)):min(c(ageMax, max(Age)))
	avgRatios <- tapply(ratio[ind], ages %% 10, mean, na.rm = TRUE)
	
	# return avg deviation for specified digit(s)
	mean(avgRatios[as.character(digit)])
}


#' calculate Noubmissi's digit heaping index
#' 
#' @description this method compares single terminal digit numerators to denominators 
#' consisting in 5-year age groups centered on the digit in question. As seen in Spoorenberg (2007)

#' @param Value a vector of demographic counts by single age
#' @param Age a vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin the lowest age included in calcs. Default 20
#' @param ageMax the upper age bound used for calcs. Default 65
#' @param digit which terminal digit do we calculate for?
#' @details  \code{ageMin} and \code{ageMax} are applied to numerator ages, not denominators.
#'  Denominators are always 5-year age groups centered on the digit in question,
#' and these therefore stretch into ages a bit higher or lower than the numerator ages.
#' @return the index value 

#' @export
#' @examples
#' Value <- c(80626,95823,104315,115813,100796,105086,97266,116328,
#' 75984,89982,95525,56973,78767,65672,53438,85014,
#' 47600,64363,42195,42262,73221,30080,34391,29072,
#' 20531,66171,24029,44227,24128,23599,82088,16454,
#' 22628,17108,12531,57325,17220,28425,16206,17532,
#' 65976,11593,15828,13541,8133,44696,11165,18543,
#' 12614,12041,55798,9324,10772,10453,6773,28358,
#' 9916,13348,8039,7583,42470,5288,5317,6582,
#' 361,17949,3650,5873,3279,3336,27368,1965,
#' 2495,2319,1335,12022,1401,1668,1360,1185,
#' 9167,424,568,462,282,6206,343,409,333,291,4137,133,169,157,89,2068,68,81,66,57)
#' Age <- 0:99
#' Noumbissi(Value, Age, digit = 0) # 2.32
#' Noumbissi(Value, Age, digit = 1) # 0.55
#' Noumbissi(Value, Age, digit = 2) # 0.73
#' Noumbissi(Value, Age, digit = 3) # 0.76
#' Noumbissi(Value, Age, digit = 4) # 0.49
#' Noumbissi(Value, Age, digit = 5) # 2.08
#' Noumbissi(Value, Age, digit = 6) # 0.66
#' Noumbissi(Value, Age, digit = 7) # 1.08  7 looks good!
#' Noumbissi(Value, Age, digit = 8) # 0.57
#' Noumbissi(Value, Age, digit = 9) # 0.59
Noumbissi <- function(Value, Age, ageMin = 20, ageMax = 65, digit = 0){
	stopifnot(length(Age) == length(Value))
	stopifnot(length(digit) == 1)
	
	numi   <- Age >= ageMin & Age <= ageMax & Age %% 10 %in% digit
	denomi <- shift.vector(numi,-2) | 
			shift.vector(numi,-1) |
			numi |
			shift.vector(numi,1) |
			shift.vector(numi,2)
	5 * sum(Value[numi]) / sum(Value[denomi])
}

#' Spoorenberg's total modified Whipple index

#' @description This index consists in calculating the Noumbissi index, Wi, for each
#' terminal digit, then taking the sum of the absolute deviations of these from unity.

#' @param Value a vector of demographic counts by single age
#' @param Age a vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin the lowest age included in calcs. Default 20
#' @param ageMax the upper age bound used for calcs. Default 65
#' @return the index value 
#' 
#' @details  \code{ageMin} and \code{ageMax} are applied to numerator ages, not denominators.
#'  Denominators are always 5-year age groups centered on the digit in question,
#' and these therefore stretch into ages a bit higher or lower than the numerator ages.
#' 
#' @export 
#' @examples
#' Value <- c(80626,95823,104315,115813,100796,105086,97266,116328,
#' 75984,89982,95525,56973,78767,65672,53438,85014,
#' 47600,64363,42195,42262,73221,30080,34391,29072,
#' 20531,66171,24029,44227,24128,23599,82088,16454,
#' 22628,17108,12531,57325,17220,28425,16206,17532,
#' 65976,11593,15828,13541,8133,44696,11165,18543,
#' 12614,12041,55798,9324,10772,10453,6773,28358,
#' 9916,13348,8039,7583,42470,5288,5317,6582,
#' 361,17949,3650,5873,3279,3336,27368,1965,
#' 2495,2319,1335,12022,1401,1668,1360,1185,
#' 9167,424,568,462,282,6206,343,409,333,291,4137,133,169,157,89,2068,68,81,66,57)
#' Age <- 0:99
#' Spoorenberg(Value, Age)

Spoorenberg <- function(Value, Age, ageMin = 20, ageMax = 65){
	
	digits <- 0:9
	Wi   <- sapply(digits, 
			Noumbissi, 
			Value = Value, 
			Age = Age, 
			ageMin = ageMin, 
			ageMax = ageMax)
	Wtot <- sum(abs(1 - Wi))
	return(Wtot)
}
