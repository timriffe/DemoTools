
# Author: tim
###############################################################################


#' Calculate Whipple's index of age heaping

#' @description Implementation following the PASEX spreadsheet SINGAGE, with some extra options for more digit checking. 
#' Whipple's index is a summary measure of age heaping on ages ending usually in the digits 0 or 5 used to determine variability in the
#' quality of age reporting between regions or countries and its evolution over time. It assumes a
#' linear distribution of ages in each five-year age range. Digits refer to terminal age digits. For example, 9, 19, 29, etc are all of digit 9.

#' @param Value numeric. A vector of demographic counts by single age.
#' @param Age numeric. A vector of ages corresponding to the lower integer bound of the counts.
#' @param ageMin integer. The lowest age included in calculations. Default 25.
#' @param ageMax integer. The upper age bound used for calculations. Default 65.
#' @param digit integer. Any digit between 0 and 9. Default \code{c(0,5)}. Otherwise it needs to be a single digit.
#'  
#' 
#' @details \code{ageMin} and \code{ageMax} refer to the bounds of the numerator, where \code{ageMax} is inclusive. 
#' The denominator looks 7 ages lower and 2 ages higher, so these ages must be available. You can get 
#' arbitrary W(i) indices by specifying other digits. Note you can only do pairs of digits
#' they are 0 and 5. Otherwise just one digit at a time. 
#' @return The value of the index.
#' @references 
#' \insertRef{GDA1981IREDA}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#' @export 
#' @examples 
#'  Age <- 0:99
#'  
#'  (w05 <- Whipple(pop1m_pasex, Age, 25, 60, digit = c(0,5))) # 2.34,replicates SINGAGE males
#'  
#'  # implements formula from Roger et al. (1981, p. 148)
#'  (w0 <- Whipple(pop1m_pasex, Age, 25, 60, digit = 0))
#'  (w5 <- Whipple(pop1m_pasex, Age, 25, 60, digit = 5)) 
#'  
#'  # Whipple types
#'  Whipple(pop1m_pasex, Age, 25, 60, digit = 3) 
Whipple <- function(Value, Age, ageMin = 25, ageMax = 65, digit = c(0,5)){
	
	stopifnot(is_single(Age))
	stopifnot(length(digit) == 1 || (all(c(0, 5) %in% digit) & length(digit == 2)))
	stopifnot(length(Value) == length(Age))
	
	numeratorind   <- Age >= ageMin & Age <= ageMax & Age %% 10 %in% digit
	# if we are checking just one digit, go down 7 up two, so that the right nr
	# of counts in denom. This per the French formulas.
	# TR: sufficient to check length here, since hard check done above
	denominatorind <- Age >= (ageMin - ifelse(length(digit) == 2,2,7)) & Age <= (ageMax + 2)
	
	whip           <- ifelse(length(digit) == 2,5,10) * sum(Value[numeratorind]) / sum(Value[denominatorind])
	
	return(whip)
}

# test to see if can override inheritParams
#' Calculate Myer's blended index of age heaping

#' @description Implementation following the PASEX spreadsheet SINGAGE. Myers' measures preferences for each of the ten possible digits
#' as a blended index. It is based on the principle that in the absence of age heaping, the aggregate population of each age ending in one of the digits
#'  0 to 9 should represent 10 percent of the total population.
#' @inheritParams Whipple
#' @param ageMax integer. The upper age bound used for calculations. Default 89.
#' 
#' @details \code{ageMax} is an inclusive upper bound, treated as interval. If you want ages
#' 20 to 89, then give \code{ageMin = 20} and \code{ageMax = 89}, not 90. \code{ageMax} may be
#' internally rounded down if necessary so that \code{ageMax - ageMin + 1} is evenly divisible by 10.
#' @return The value of the index. 
#' @references 
#' \insertRef{myers1954accuracy}{DemoTools}
#' \insertRef{spoorenberg2007quality}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#' @export 
#' @examples
#' Age <- 0:99
#' Myers(pop1m_pasex, Age, 10, 89) 
#' 
#' Myers(pop1m_pasex, Age, 10, 90) * 2 #47.46, replicates SINGAGE males

Myers <- function(Value, Age, ageMin = 10, ageMax = 89){
	stopifnot(is_single(Age))
	# hard code period to 10 for digits
	period  <- 10
	
	# must be of same length for indexing
	stopifnot(length(Value) == length(Age))
	#stopifnot(ageMin %% period == 0 & ageMax %% period == 0)
	
	# ageMax dynamically rounded down 
	# to a total span divisible by period
	Diff    <- ageMax - ageMin + 1
	AgeMax  <- ageMin + Diff - Diff %% period 
	
	# may as well be certain here
	stopifnot(ageMax <= max(Age))
	
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

#' @description Bachi's index involves applying the Whipple method repeatedly to determine the extent of preference for each final digit.
#' Similarly to Myers', it equals the sum of the positive deviations from 10 percent. It has a theoretical range from 0 to 90, 
#' and 10 is the expected value for each digit.
#' Two Implementations: one following the PASEX spreadsheet SINGAGE (\code{pasex = TRUE}), with ages hard-coded, and another with flexible upper and lower age bounds, but that does not match the PASEX implementation.

#' @inheritParams Whipple
#' @param pasex logical. Whether or not reproduce the specific age weightings in the PASEX spreadsheet. Default \code{FALSE}.
#' 
#' @details \code{ageMax} is an inclusive upper bound, treated as interval. If you want ages
#' 30 to 79, then give \code{ageMin = 30} and \code{ageMax = 79}, not 80. These are only heeded if \code{pasex = FALSE}.
#' @return The value of the index.
#' @references 
#' \insertRef{PAS}{DemoTools}
#' \insertRef{bachi1951tendency}{DemoTools}
#' \insertRef{shryock1973methods}{DemoTools}
#' @export 
#' @examples
#' Age <- 0:99
#' 
#' Bachi(pop1m_pasex, Age, ageMin = 20, ageMax = 79, pasex = TRUE) # reproduces PASEX SINGAGE
#' Bachi(pop1m_pasex, Age, ageMin = 20, ageMax = 79) # default simpler

Bachi <- function(Value, Age, ageMin = 30, ageMax = 79, pasex = FALSE){
	stopifnot(length(Age) == length(Value))
	stopifnot(is_single(Age))
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
		# cheap way to make upper bound inclusive
		ageMax       <- ageMax + 1
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
#' @inheritParams Whipple
#' @param terms integer. Length of the (centered) moving average be. Default 5.
#' @param digit integer. Any digit 0-9. Default 0.


#' @details \code{digit} could also be a vector of digits, but the more digits one includes (excepting 0 and 5) the closer the index will get to 1. 
#' It is therefore recommended for single digits, or else \code{c(0,5)}. 
#' \code{ageMax} is an inclusive upper bound, treated as interval.
#'  If you want ages 20 to 89, then give \code{ageMin = 20} and \code{ageMax = 89}, not 90. By default all available ages greater than or equal to \code{ageMin} are used. 
#' 
#' @return The value of the index.
#' 
#' @references 
#' \insertRef{coale1991effect}{DemoTools}
#' @export
#' @examples 
#' Age <- 0:99
#' CoaleLi(pop1m_pasex, Age, 65, 95, 5, 0) # 3.7
#' CoaleLi(pop1m_pasex, Age, 65, 95, 5, 5) # 3.5 almost just as high

CoaleLi <- function(Value, Age, ageMin = 60, ageMax = max(Age), terms = 5, digit = 0){
	stopifnot(is_single(Age))
	stopifnot(length(Age) == length(Value))
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


#' calculate Noumbissi's digit heaping index
#' 
#' @description Noumbissi's method improves on Whipple's method by extending its basic principle to all ten digits. It compares single terminal digit numerators to denominators 
#' consisting in 5-year age groups centered on the terminal digit of age in question. For example, 9, 19, 29, etc are all of digit 9.

#' @param Value numeric. A vector of demographic counts by single age.
#' @param Age numeric. A vector of ages corresponding to the lower integer bound of the counts.
#' @param ageMin integer. The lowest age included in calculations. Default 20.
#' @param ageMax integer. The upper age bound used for calculations. Default 64.
#' @param digit integer. Any digit 0-9. Default 0.
#' @details  \code{ageMin} and \code{ageMax} are applied to numerator ages, not denominators.
#'  Denominators are always 5-year age groups centered on the digit in question,
#' and these therefore stretch into ages a bit higher or lower than the numerator ages. \code{ageMax} is an inclusive upper bound, treated as interval. 
#' If you want ages 20 to 64, then give \code{ageMin = 20} and \code{ageMax = 64}, not 65. 
#' @return The value of the index. 
#' @references 
#' \insertRef{noumbissi1992indice}{DemoTools}
#' \insertRef{spoorenberg2007quality}{DemoTools}
#' @export
#' @examples
#' Age <- 0:99
#' Noumbissi(pop1m_pasex, Age, digit = 0) # 2.32
#' Noumbissi(pop1m_pasex, Age, digit = 1) # 0.55
#' Noumbissi(pop1m_pasex, Age, digit = 2) # 0.73
#' Noumbissi(pop1m_pasex, Age, digit = 3) # 0.76
#' Noumbissi(pop1m_pasex, Age, digit = 4) # 0.49
#' Noumbissi(pop1m_pasex, Age, digit = 5) # 2.08
#' Noumbissi(pop1m_pasex, Age, digit = 6) # 0.66
#' Noumbissi(pop1m_pasex, Age, digit = 7) # 1.08  7 looks good!
#' Noumbissi(pop1m_pasex, Age, digit = 8) # 0.57
#' Noumbissi(pop1m_pasex, Age, digit = 9) # 0.59

Noumbissi <- function(Value, Age, ageMin = 20, ageMax = 64, digit = 0){
	stopifnot(is_single(Age))
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

#' @description Using the digit-specific modified Whipple's index, this index summarizes
#' all age preference and avoidance effects by taking the sum of the absolute differences between digit-specific Whipple's index and 1 (counting all
#' differences as positive).
#' @param Value numeric. A vector of demographic counts by single age.
#' @param Age numeric. A vector of ages corresponding to the lower integer bound of the counts.
#' @param ageMin integer. The lowest age included in calculations. Default 20.
#' @param ageMax integer. The upper age bound used for calculations. Default 64.
#' 
#' @details  \code{ageMin} and \code{ageMax} are applied to numerator ages, not denominators.
#'  Denominators are always 5-year age groups centered on the digit in question,
#' and these therefore stretch into ages a bit higher or lower than the numerator ages. \code{ageMax} is an inclusive upper bound, treated as interval. 
#' If you want ages 20 to 64, then give \code{ageMin = 20} and \code{ageMax = 64}, not 65. 
#' @return The value of the index.
#' @references 
#' \insertRef{spoorenberg2007quality}{DemoTools}
#' @export 
#' @examples
#' Age <- 0:99
#' Spoorenberg(pop1m_pasex, Age)

Spoorenberg <- function(Value, Age, ageMin = 20, ageMax = 64){
	stopifnot(length(Age) == length(Value))
	stopifnot(is_single(Age))
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
