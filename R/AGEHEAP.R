
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


#' Kannisto's age heaping index
#' 
#' @description This age heaping index is used for particular old ages, such as 90, 95, 100, 105, and so forth.
#' 
#' @param Value numeric. A vector of demographic counts or rates by single age.
#' @param Age integer. A vector of ages corresponding to the lower integer bound of the counts.
#' @param Agei integer. The age on which the index is centered.
#' @param pow either \code{"exp"} (default) or a power such as 2. See details
#' 
#' @details The index looks down two ages and up two ages, so the data must accommodate that range. The denominator is a mean of the counts in the sourrounding 5 single ages. The kind of mean can be controlled with the \code{pow} argument. By default, this takes the antilog of the arithmetic mean of the natural log of the five denominator counts. That will fail if one of the counts is equal to 0. In such cases, another power, such as 2 or 10 or 100 may be used, which is more robust to 0s. The higher the power, the closer the result will resemble the default output. If \code{pow=="exp"} but a 0 is detected among the denominator ages, then \code{pow} is assigned a value of 1000. \code{pow=1} would imply an arithmetic mean in the denominator.
#' @return The value of the index.
#' 
#' @references 
#' \insertRef{kannisto1999assessing}{DemoTools}
#' @export
#' 
#' @examples
#'Age <- 0:99
#' KannistoHeap(pop1m_pasex, Age, 90)
#' KannistoHeap(pop1m_pasex, Age, 95)
#' KannistoHeap(pop1m_pasex, Age, 95, pow = 2) # geometric mean in denom
#' KannistoHeap(pop1m_pasex, Age, 95, pow = 1000) # similar, but robust to 0s
#' KannistoHeap(pop1m_pasex, Age, 95, pow = 1) # arithmetic mean in denom
#' pop1m_pasex[Age==95] / mean(pop1m_pasex[Age >= 93 & Age <= 97])
KannistoHeap <- function(Value, Age, Agei = 90,pow="exp"){
	stopifnot(length(Agei) == 1)
	stopifnot(is_single(Age))
	denomi <- Age %in% ((Agei - 2):(Agei + 2))
	if (any(Value[denomi] == 0)){
		pow <- 1000
	}
	if (pow == "exp"){
		denom  <- exp(mean(log(Value[denomi])))
	} else {
		denom  <- mean((Value[denomi])^(1/pow))^pow
	}
	Value[Age == Agei] / denom
}


#' Calculate Jdanov's old-age heaping index
#' 
#' @description This is a slightly more flexible implementation of Jdanov's formula,
#' with defaults set to match his parameters. The numerator is the sum of 
#' (death counts) in ages 95, 100, and 105. The denominator consists in the sum of the 5-year age groups centered around each of the numerator ages.
#' It probably only makes sense to use this with the default values, however. Used with a single age
#' in the numerator, it is almost the same as \code{Noumbissi()}, except here we pick out particular ages,
#' whereas Noumbissi picks out terminal digits.
#' 
#' @param Value numeric. A vector of demographic counts by single age.
#' @param Age integer. A vector of ages corresponding to the lower integer bound of the counts.
#' @param Agei integer. A vector of ages to put in the numerator, default \code{c(95,100,105)}.
#' 
#' @return The value of the index.
#' @references 
#' \insertRef{jdanov2008beyond}{DemoTools}
#' @export
#' @examples
#'Value <-c(8904, 592, 354, 299, 292, 249, 222, 216, 181, 169, 151, 167, 
#'		170, 196, 249, 290, 425, 574, 671, 724, 675, 754, 738, 695, 597, 
#'		498, 522, 479, 482, 478, 558, 582, 620, 606, 676, 768, 862, 952, 
#'		1078, 1215, 1215, 1357, 1470, 1605, 1723, 1782, 1922, 2066, 2364, 
#'		2561, 2476, 1674, 1664, 1616, 1808, 3080, 3871, 4166, 4374, 4707, 
#'		5324, 5678, 6256, 6382, 6823, 7061, 7344, 8149, 8439, 8308, 8482, 
#'		8413, 8157, 7945, 7503, 7164, 7289, 7016, 6753, 6906, 6797, 6624, 
#'		6416, 5811, 5359, 4824, 4277, 3728, 3136, 2524, 2109, 1657, 1235, 
#'		924, 667, 465, 287, 189, 125, 99, 80, 24, 10, 7, 3, 1, 0, 1, 
#'		1, 0, 0)
#'Age <- 0:110
#'Jdanov(Value, Age, Agei = c(95,100,105))

Jdanov <- function(Value, Age, Agei = seq(95,105,by=5)){
	stopifnot(is_single(Age))
	numi   <- Age %in% Agei
	
	# this way of doing the denom makes it possible to 
	# have duplicate values in the denom, for the case that
	# the numerator ages are closer together than 5. Innocuous otherwise
	denom <- sum(Value[shift.vector(numi, -2)]) +
			sum(Value[shift.vector(numi, -1)]) +
			sum(Value[numi]) +
			sum(Value[shift.vector(numi, 1)]) +
			sum(Value[shift.vector(numi, 2)]) 
	
	500 * sum(Value[numi]) / denom
	
}


#' Induce heaping on terminal digits 0 and 5
#' @description For single age count data, perturb the data in such a way as to induce heaping on ages ending in 0 and 5. 
#' This is a common phenomenon that several age heaping evaluation methods are designed to test for. 
#' In order to estimate how these methods respond to different degrees of heaping in a systematic way, a function such as
#'  this may be useful. The way this works is purely a guess, and no checks are in place, so use with caution.
#' @details We use pascal weights, \code{c(1,4,6,4,1)} scaled to argument \code{p5} and \code{c(1,6,15,20,15,6,1)} 
#' scaled to argument \code{p0} to shift individuals from surrounding ages to ages endings in 5 and 0, respectively.
#'  This assumes that ages closer to digits 0 or 5 are more likely to be declared as 0 or 5, that this bias is symmetrical,
#'   and that such rounding can come from farther away ages for terminal digit 0 than for 5. 
#'   The user can also make heaping scale differently for 0s and 5s, but there is no control over the shape of bias.
#' @export
#' @inheritParams Bachi
#' @param p0 total size of bias to round to 0s. Default 2.
#' @param p5 total size of bias to round to 5s. Defaults to \code{p0}.
#' @param ageMin integer ending in 0 or 5. The lowest age to be heaped upon, default 25.
#' @param ageMax integer ending in 0 or 5. The highest age to be heaped upon, defaults to highest age evenly divisible by 10.
#' @return Value numeric vector perturbed to look like age heaping as happened.
#' @examples
#' # example to show what we're talking about.
#' # pop1m_pasex is already quite heaped:
#' Age <- 0:99
#' plot(Age,pop1m_pasex)
#' # here it is again, smoothed:
#' smoothed <- sprague(
#' 		agesmth(pop1m_pasex, 
#' 				Age, 
#' 				method = "Strong", 
#' 				OAG = FALSE, 
#' 				young.tail = "Arriaga"),
#' 		OAG = FALSE)
#' lines(Age, smoothed)
#' # an OK approximation for testing purposes.
#' points(Age,heapify(smoothed,Age=0:99,1.8,1.1,ageMin=20),pch="x")
#' 

heapify <- function(Value, Age, p0=2, p5=p0, ageMin = 25, ageMax = max(Age[Age %% 5 == 0])){
	# pascal weights
	# x3,x4,x5,x6,x7
	fivepdf <- rescale.vector(c(1,4,6,4,1)) 
	#x7,x8,x9,x0,x1,x2,x3 
	tenpdf  <- rescale.vector(c(1,6,15,20,15,6,1)) * p0
	
	# center ages:
	
	A10 <- Age[Age %% 10 == 0 & Age >= ageMin & Age <= ageMax]
	A5  <- Age[Age %% 5 == 0 & !Age%in%A10 & Age >= ageMin & Age <= ageMax]
	
	# rest could happen in a matrix, but this might be clearer to read:
	ai <- sort(c(A5,A10))
	for (a in ai){
		if (a %% 10 == 0){
			dist <- 3
			pdf.i <- tenpdf * p0
		} else {
			dist <- 2
			pdf.i <- fivepdf * p5
		}
		sendi <- Age >= (a - dist) & Age <= (a + dist)
		Vsend <- Value[sendi] * pdf.i
		Value[sendi] <- Value[sendi] - Vsend
		Value[Age == a] <- Value[Age == a] + sum(Vsend)
	}
	Value
}


#' Detect if heaping is worse on terminal digits 0s than on 5s
#' @description Ages ending in 0 often have higher apparent heaping than ages ending in 5. 
#' In this case, data in 5-year age bins might show a sawtooth pattern. 
#' If heaping ocurrs in roughly the same amount on 0s and 5s, then it may be sufficient to group data into 5-year
#'  age groups and then graduate back to single ages. However, if heaping is worse on 0s, then this procedure tends 
#'  to produce a wavy pattern in count data, with 10-year periodicity. In this case it is recommended to use one of 
#'  the methods of \code{agesmth()} as an intermediate step before graduation. 
#' @details Data is grouped to 5-year age bins. The ratio of each value to the average of its neighboring values 
#' is calculated. If 0s have stronger attraction than 5s then we expect these ratios to be >1 for 0s and <1 for 5s.
#'  Ratios are compared within each 10-year age group in the evaluated age range. If in the evaluated range there are
#'  at most two exceptions to this rule (0s>5s), then the ratio of the mean of these ratios is returned, and it is 
#'  recommended to use a smoother method. Higher values suggest use of a more aggressive method. 
#'  This approach is only slightly different from that of Feeney, as implemented in the \code{zigzag()} functions. 
#'  This is not a general measure of roughness, but rather an indicator of this particular pattern of age attraction. 
#' @export
#' @inheritParams heapify
#' @param ageMin integer evenly divisible by 10. Lower bound of evaluated age range, default 40.
#' @param ageMax integer evently divisibly by 5. Upper bound of evaluated age range, defaults to highest age evenly divisible by 10.
#' @return \code{FALSE} if sawtooth pattern is not detected, numeric otherwise.
#' @references
#' \insertRef{feeney1979}{DemoTools}
#' Feeney, G. 2013 "Removing "Zigzag" from Age Data," http://demographer.com/white-papers/2013-removing-zigzag-from-age-data/

#' @examples
#' Age <- 0:99
#' smoothed <- sprague(
#' 		agesmth(pop1m_pasex, 
#' 				Age, 
#' 				method = "Strong", 
#' 				OAG = FALSE, 
#' 				young.tail = "Arriaga"),
#' 		OAG = FALSE)
#' # not saw-tooth jagged
#' zero_pref_sawtooth(smoothed, Age)
#' # saw-tooth pattern detected in older ages
#' zero_pref_sawtooth(pop1m_pasex, Age)
#' # heaped, but no 0>5 preference
#' h1 <- heapify(smoothed, Age, p0 = 1, p5 = 1)
#' # heaping progressively worse on 0s than on 5s.
#' h2 <- heapify(smoothed, Age, p0 = 1.2, p5 = 1)
#' h3 <- heapify(smoothed, Age, p0 = 1.5, p5 = .8)
#' h4 <- heapify(smoothed, Age, p0 = 2, p5 = .5)
#' 
#' \dontrun{
#' 	plot(Age, smoothed, type='l')
#' 	lines(Age, h1,col="blue")
#' 	lines(Age, h2,col="green")
#' 	lines(Age, h3,col="red")
#' }
#' zero_pref_sawtooth(h1, Age)
#' # 0-preference > 5 pref
#' zero_pref_sawtooth(h2, Age)
#' # increasing values
#' zero_pref_sawtooth(h3, Age)
#' zero_pref_sawtooth(h4, Age,ageMin=35)

zero_pref_sawtooth <- function(Value, Age, ageMin = 40, ageMax = max(Age[Age %% 5 == 0])){
	
	# rather than stopifnot() check, just make it work.
	ageMin <- ageMin - ageMin %% 10
	# group to 5-year ages if not already.
	VH5    <- groupAges(Value, Age, 5)
	A5     <- names2age(VH5)
	# elementwise matched to avg of adjacent values
	adj2   <- avg_adj(VH5)
	# evalauated age range
	ai     <- A5 >= ageMin & A5 <= ageMax
	# matrix of ratios, 0s in row 1, 5s in row 2
	m05    <- suppressWarnings(matrix((VH5 / adj2)[ai], nrow = 2))
	# ensure no recycled values used
	if (sum(ai) %% 2 != 0 | any(is.na(m05))){
		m05    <- m05[, -ncol(m05)]	
	}
	
	# need rather consistent x0 > x5 pattern
	if (sum(diff(sign(log(m05))) == -2) < (ncol(m05) - 2)){
		# i.e. it could be very rough still, but not necessarily
		# a regular sawtooth pattern. Possibly a visual assessment needed,
		# since could be real pattern, or other pattern of roughness that also
		# requires smoothing.
		return(FALSE)
	}
	
	# mean of 0s divided by mean of 5s, that simple.
	1 / ratx(rowMeans(m05, na.rm = TRUE))
}

#' Evaluate roughness of data in 5-year age groups
#' @description For a given age-structured vector of counts, how rough is data after grouping to 5-year age bins?
#'  Data may require smoothing even if there is no detectable sawtooth pattern. It is best to use the value in this 
#'  method together with visual evidence to gauage whether use of \code{agesmth()} is recommended.
#' @details First we group data to 5-year age bins. Then we take first differences (d1) of these within the evaluated
#'  age range. Then we smooth first differences (d1s) using a generic smoother (\code{ogive()}).
#'   Roughness is defined as the mean of the absolute differences between \code{mean(abs(d1 - d1s) / abs(d1s))}. 
#'   Higher values indicate rougher data, and may suggest more aggressive smoothing. Just eyeballing, one could consider
#'    smoothing if the returned value is greater than ca 0.2, and values greater than 0.5 already highly recommend
#'     it (pending visual verification).
#' @export
#' @inheritParams zero_pref_sawtooth
#' @param ageMin integer evenly divisible by 5. Lower bound of evaluated age range, default 20.
#'
#' @examples 
#' Age <- 0:99
#' smoothed <- sprague(
#' 		agesmth(pop1m_pasex, 
#' 				Age, 
#' 				method = "Strong", 
#' 				OAG = FALSE, 
#' 				young.tail = "Arriaga"),
#' 		OAG = FALSE)
#' # not very rough, no need to smooth more
#' five_year_roughness(smoothed, Age)
#'  # quite rough, even after grouping to 5-year ages
#' five_year_roughness(pop1m_pasex, Age)
#' # heaped, but no 0>5 preference
#' h1 <- heapify(smoothed, Age, p0 = 1, p5 = 1)
#' # heaping progressively worse
#' h2 <- heapify(smoothed, Age, p0 = 1.2, p5 = 1.2)
#' h3 <- heapify(smoothed, Age, p0 = 1.5, p5 = 1.5)
#' h4 <- heapify(smoothed, Age, p0 = 2, p5 = 2)
#' h5 <- heapify(smoothed, Age, p0 = 2.5, p5 = 2)
#'\dontrun{
#'	#cols <- RColorBrewer::brewer.pal(7,"Reds")[3:7]
#'  cols <-  c("#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#99000D")
#'  A5 <- seq(0,95,by=5)
#' 	plot(A5, groupAges(smoothed), type='l',xlim=c(20,80),ylim=c(0,3e5))
#'	lines(A5, groupAges(h1),col=cols[1])
#' 	lines(A5, groupAges(h2),col=cols[2])
#' 	lines(A5, groupAges(h3),col=cols[3])
#'	lines(A5, groupAges(h4),col=cols[4])
#'	lines(A5, groupAges(h5),col=cols[5])
#'}
#'five_year_roughness(smoothed, Age)
#'five_year_roughness(h1, Age)
#'five_year_roughness(h2, Age)
#'five_year_roughness(h3, Age)
#'five_year_roughness(h4, Age)
#'five_year_roughness(h5, Age)
#' 

five_year_roughness <- function(Value, Age, ageMin = 20, ageMax = max(Age[Age %% 5 == 0])){
	# rather than stopifnot() check, just make it work.
	ageMin <- ageMin - ageMin %% 5
	# group to 5-year ages if not already.
	VH5    <- groupAges(Value, Age, 5)
	A5     <- names2age(VH5)
	# evalauated age range
	ai     <- A5 >= ageMin & A5 <= ageMax
	
	d1     <- diff(VH5[ai]) 
	# compare with something smooth, loess by default.
	d1s    <- agesmth1(Value = d1, Age = A5[ai][-1], OAG = FALSE)
	
	mean(abs(d1 - d1s) / abs(d1s))
}
