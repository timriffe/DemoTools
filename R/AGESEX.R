
# Author: tim
###############################################################################

#' calculate the PAS age ratio score
#' @description A single ratio is defined as the 100 times twice the size of an age group 
#' to it's neighboring two age groups. Under uniformity these would all be 100. The average 
#' absolute deviation from 100 defines the index. This comes from the PAS spreadsheet called AGESEX

#' @param Value numeric. A vector of demographic counts in 5-year age groups.
#' @param Age an integer vector of ages corresponding to the lower integer bound of the counts.
#' @param ageMin integer. The lowest age included in calcs. Default 0.
#' @param ageMax integer. The upper age bound used for calcs. Default \code{max(Age)}.
#' @param method character. Either \code{"UN"} (default), \code{"Zelnick"}, or \code{"Ramachandran"}
#' @param OAG logical. default \code{TRUE}. Is the top age group open?

#' @details Age groups must be of equal intervals. Five year age groups are assumed.
#'  We also assume that the final age group is open, unless \code{ageMax < max(Age)}. Setting \code{OAG = FALSE} will override this and potentially include \code{max(Age)} in calculations.

#' @return the value of the index
#' @export
#' 
#' @references 
#' \insertRef{accuracyun1952}{DemoTools}
#' \insertRef{ramachandran1967}{DemoTools}
#' @examples 
#' # data from PAS spreadsheet AGESEX.xlsx
#' Males   <- c(4677000,4135000,3825000,3647000,3247000,2802000,2409000,2212000,
#' 		1786000,1505000,1390000,984000,745000,537000,346000,335000)
#' Females <- c(4544000,4042000,3735000,3647000,3309000,2793000,2353000,2112000,
#' 		1691000,1409000,1241000,887000,697000,525000,348000,366000)
#' Age     <- seq(0, 75, by = 5)
#' 
#' ageRatioScore(Males, Age, ageMax = 75)    # 3.9, matches PAS
#' ageRatioScore(Females, Age, ageMax = 75)  # 3.66 matches PAS
#' ageRatioScore(Females, Age, ageMax = 75, method = "Ramachandran") # 1.8
#' ageRatioScore(Females, Age, ageMax = 75, method = "Zelnick")      # 2.4
#' 
#' # de facto unit test:
#' stopifnot(round(ageRatioScore(Males, Age, ageMax = 75) ,3) == 3.907)
#' stopifnot(round(ageRatioScore(Females, Age, ageMax = 75) ,3) == 3.655)

ageRatioScore <- function(Value, Age, ageMin = 0, ageMax = max(Age), method = "UN", OAG = TRUE){
	stopifnot(length(Value) == length(Age))
	method          <- tolower(method)
	stopifnot(method %in% c("un","zelnick","ramachandran"))

	N <- length(Value)
	if (OAG){
		Age        <- Age[-N]
		Value      <- Value[-N]
		if (ageMax > max(Age)){
			ageMax <- max(Age)
		}
	}
	
	# cut down data if necessary
	keep            <- Age >= ageMin & Age <= ageMax
	Value           <- Value[keep]
	Age             <- Age[keep]
	
	# selectors for numerators and denominators
	numi            <- Age > ageMin & Age < max(Age)
	denomleft       <- shift.vector(numi, -1) 
	denommiddle     <- numi
	denomright      <- shift.vector(numi, 1)
	
	# one difference between UN, Zelnick and Ramachandran methods is whether and how
	# much to add the numerator into the denominator, and also how much weight to give
	# to the numerator (2, 3, or 4 times)
	middlemult      <- ifelse(method == "un", 0, ifelse(method == "zelnick", 1, 2))
	nummult         <- ifelse(method == "un", 200, ifelse(method == "zelnick", 300, 400))
	
	numerator       <- nummult * Value[numi]
	denominator     <- Value[denomleft] + Value[denomright] + Value[denommiddle] * middlemult
	# ratio of each age to the age above and below it
	ageratio        <- numerator / denominator
	# absolute deviatios from uniformity
	absres          <- abs(ageratio - 100)
	# average absolute deviation
	sum(absres) / length(absres)
}

#' calculate the PAS sex ratio score
#' @description A single ratio is defined as males per 100 females.
#' Taking the first difference of the ratios over age would give 0s
#' if the ratio were constant. The average absolute difference over
#' age defines the index. This comes from the PAS spreadsheet called AGESEX.

#' @param Males numeric. A vector of population counts for males in 5-year age groups
#' @param Females numeric. A vector of population counts for females in 5-year age groups
#' @param Age integer. A vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin integer. The lowest age included in calcs. Default 0
#' @param ageMax integer. The upper age bound used for calcs. Default \code{max(Age)}
#' @param OAG logical. default \code{TRUE}. Is the top age group open?
#' @details Age groups must be of equal intervals. Five year age groups are assumed.
#'  We also assume that the final age group is open, unless \code{ageMax < max(Age)}. The method
#' argument determines the weighting of numerators and denominators, where the UN method puts
#' twice the numerator over the sum of the adjacent ages classes, Zelnich does thrice the 
#' numerator over the sum of the whole range from the next lowest to the next highest age, and 
#' Ramachandran does four times the numerator over the same sum, but with the central age 
#' double-counted in the numerator. Ramachandran is therefore less judgemental, so to speak.

#' @return the value of the index
#' @references 
#' \insertRef{accuracyun1952}{DemoTools}
#' @export
#' @examples 
#' Males   <- c(4677000,4135000,3825000,3647000,3247000,2802000,2409000,2212000,
#' 		1786000,1505000,1390000,984000,745000,537000,346000,335000)
#' Females <- c(4544000,4042000,3735000,3647000,3309000,2793000,2353000,2112000,
#' 		1691000,1409000,1241000,887000,697000,525000,348000,366000)
#' Age     <- seq(0, 75, by = 5)
#' sexRatioScore(Males, Females, Age)  # 2.2, matches PAS
#' # de facto unit test
#' stopifnot(round(sexRatioScore(Males, Females, Age),3) == 2.249)
#' stopifnot(round(sexRatioScore(Males, Females, Age, ageMax = 70),3) == 2.249)

sexRatioScore <- function(Males, Females, Age, ageMin = 0, ageMax = max(Age), OAG = TRUE){
	stopifnot(length(Males) == length(Age) & length(Males) == length(Females))
	
	# handle open age group
	N <- length(Males)
	if (OAG){
		Age        <- Age[-N]
		Males      <- Males[-N]
		Females    <- Females[-N]
		if (ageMax > max(Age)){
			ageMax <- max(Age)
		}
	}
	# now make 
	keep      <- Age >= ageMin & Age <= ageMax
	Males     <- Males[keep]
	Females   <- Females[keep]
	
	ratio     <- Males / Females
	ratiodiff <- 100 * diff(ratio)
	absdiff   <- abs(ratiodiff)
	sum(absdiff) / length(absdiff)
}


#' calculate an age-sex accuracy index
#' @description This index is a composite consisting in the sum of thrice the sex 
#' ratio index plus the age ratio index for males and females. This function is therefore
#' a wrapper to \code{ageRatioScore()} and \code{sexRatioScore()}.

#' @param Males numeric. A vector of population counts for males in 5-year age groups
#' @param Females numeric. A vector of population counts for females in 5-year age groups
#' @param Age integer. A vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin integer. The lowest age included in calcs. Default 0
#' @param ageMax integer. The upper age bound used for calcs. Default \code{max(Age)}
#' @param method character. Either \code{"UN"} (default), \code{"Zelnick"}, \code{"Ramachandran"}, or \code{"das gupta"}
#' @param adjust logical do we adjust the measure when population size is under one million?
#' @param OAG logical. default \code{TRUE}. Is the top age group open?

#' @details Age groups must be of equal intervals. Five year age groups are assumed.
#' If the final element of \code{Males} and \code{Females} is the open age group,
#' then either make sure \code{ageMax} is lower than it, or leave \code{OAG} as \code{TRUE} so that it is properly removed for calcs. The method argument 
#' is passed to \code{ageRatioScore()}, where it determines weightings of numerators and denominators, 
#' except in the case of Das Gupta, where it's a different method entirely (see \code{ageSexAccuracyDasGupta()}.

#' @return the value of the index
#' @export
#' @references 
#' \insertRef{accuracyun1952}{DemoTools}
#' \insertRef{dasgupta1955}{DemoTools}
#' 
#' @examples
#' Males   <- c(4677000,4135000,3825000,3647000,3247000,2802000,2409000,2212000,
#' 		1786000,1505000,1390000,984000,745000,537000,346000,335000)
#' Females <- c(4544000,4042000,3735000,3647000,3309000,2793000,2353000,2112000,
#' 		1691000,1409000,1241000,887000,697000,525000,348000,366000)
#' Age     <- seq(0, 75, by = 5)
#' ageSexAccuracy(Males, Females, Age)    # 14.3, matches PAS
#' ageSexAccuracy(Males, Females, Age, ageMax = 65)
#' ageSexAccuracy(Males, Females, Age, ageMax = 65, adjust = FALSE)
#' ageSexAccuracy(Males, Females, Age, method = "Zelnick")
#' ageSexAccuracy(Males, Females, Age, method = "Ramachandran")
#' # Das Gupta not a comparable magnitude, FYI.
#' ageSexAccuracy(Males, Females, Age, method = "Das Gupta")
#' 
#' # de facto unit test:
#' stopifnot(round(ageSexAccuracy(Males, Females, Age),3) == 14.308)
ageSexAccuracy <- function(
		Males, 
		Females, 
		Age, 
		ageMin = 0, 
		ageMax = max(Age), 
		method = "UN", 
		adjust =  TRUE, 
		OAG = TRUE){
	method <- tolower(method)
	stopifnot(method %in% c("un","zelnick","ramachandran", "das gupta", "dasgupta"))
	
	# Das Gupta method treated differently,
	# since it prescribes all the details.
	if (method %in% c("das gupta", "dasgupta")){
		ind <- ageSexAccuracyDasGupta(Males = Males,
				Females = Females,
				Age = Age,
				ageMin = ageMin,
				ageMax = ageMax,
				OAG = OAG)
		# early return in this case
		return(ind)
	}
	
	SR <- sexRatioScore(
			Males = Males, 
			Females = Females, 
			Age = Age, 
			ageMin = ageMin, 
			ageMax = ageMax,
			OAG = OAG)
	MA <- ageRatioScore(
			Value = Males,
			Age = Age,
			ageMin = ageMin, 
			ageMax = ageMax,
			method = method,
			OAG = OAG)
	FA <- ageRatioScore(
			Value = Females,
			Age = Age,
			ageMin = ageMin, 
			ageMax = ageMax,
			method = method,
			OAG = OAG)
	# calculate index:
	ind <- 3 * SR + MA + FA
	tot <- sum(Males) + sum(Females)
	if (adjust & tot < 1e6){
		ind <- ind - 3500 / sqrt(tot) + 3.5
	}
	return(ind)
}

#' Calculate Ajit Das Gupta's (1995) age sex accuracy index
#' @description Given population counts in 5-year age groups for males and females, follow Ajit Das Gupta's steps
#'  to calculate a composite index of the quality of the age and sex structure for a given population.
#' 
#' @param Males numeric. A vector of population counts for males by age
#' @param Females numeric. A vector of population counts for females by age
#' @param Age integer. A vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin integer. The lowest age included in calcs. Default 0
#' @param ageMax integer. The upper age bound used for calcs. Default \code{max(Age)}
#' @param OAG logical. default \code{TRUE}. Is the top age group open?
#' 
#' @references 
#' \insertRef{dasgupta1955}{DemoTools}
#' 
#' @details It is assumed that the terminal age group is open, in which case it's ignored.  Set \code{OAG = FALSE} if the top age is indeed a closed interval that you want included in calculations. If \code{ageMax == max(Age)} and \code{OAG} is \code{TRUE}, then \code{ageMax} gets decremented one age class.
#' 
#' @export
#' @examples
#' # data from table for South West Africa (1946) given in reference
#' Males   <- c(2365, 2320, 1859, 1554, 1758, 1534, 1404, 1324,
#'               1118, 872, 795, 745, 743, 574)
#' Females <- c(2244, 2248, 1773, 1594, 1616, 1510, 1478, 1320,
#'               1085, 858, 768, 726, 533, 282)
#' Age     <- seq(0, 65, by = 5)
#' ageSexAccuracyDasGupta(Males, Females, Age)
#' # this method is not on the same scale as the others, so don't directly compare.
#' ageSexAccuracy(Males, Females, Age, method = "das gupta")

ageSexAccuracyDasGupta <- function(Males, Females, Age, ageMin = 0, ageMax = max(Age), OAG = TRUE){
	# handle open age group
	N <- length(Males)
	if (OAG){
		Age        <- Age[-N]
		Males      <- Males[-N]
		Females    <- Females[-N]
		if (ageMax > max(Age)){
			ageMax <- max(Age)
		}
	}
	# cut down data:
	keep      <- Age >= ageMin & Age <= ageMax
	Males     <- Males[keep]
	Females   <- Females[keep]
	# sex ratio pars
	Ri        <- Males / Females
	ratiodiff <- 100 * diff(Ri)
	absdiff   <- abs(ratiodiff)
	
	N         <- length(Males)
	n         <- length(ratiodiff)
	
	d         <- sum(absdiff)
	
	#i         <- mean(d - absdiff)
	
	# rolling sums
	mm        <- Males[-1] + Males[-N]
	mf        <- Females[-1] + Females[-N]
	
	# ratios of first and terminal rolling sums
	ri        <- 100 * mm[1] / mf[1]
	rt        <- 100 * mm[n] / mf[n]
	
	# Das Gupta's sex ratio index
	i         <- (d - abs(ri - rt)) / n
	
	# --------------------------
	# Das Gupta's age ratio index
	
	Rmi        <- 200 * Males[-1] / mm
	Rfi        <- 200 * Females[-1] / mf
	
	# sum of absolute deviations from unity ratio
	Dm         <- sum(abs(Rmi - 100))
	Df         <- sum(abs(Rfi - 100))
	
	# ratio of highest to first rolling sum
	Rhm        <- 200 * max(Males) / mm[1]
	Rhf        <- 200 * max(Females) / mf[1]
	
	# ratio of last to first rolling sum
	Rtm        <- 200 * Males[N] / mm[1]
	Rtf        <- 200 * Females[N] / mf[1]
	
	# age ratio indices
	Im         <- Dm / (Rhm - 100 + Rhm - Rtm)
	If         <- Df / (Rhf - 100 + Rhf - Rtf)
	
	# calculate the composite index
	i * Im * If
}

#' Kannisto's age heaping index
#' 
#' @description This age heaping index is used for particular old ages, such as 90, 95, 100, 105, and so forth.
#' 
#' @param Value a vector of demographic counts or rates by single age
#' @param Age a vector of ages corresponding to the lower integer bound of the counts
#' @param Agei the age on which the index is centered.
#' 
#' @details The index looks down two ages and up two ages, so the data must accomodate that range.
#' @return the index value.
#' 
#' @references 
#' Kannisto (1999) Assessing the Information on Age at Death of Old Persons in National Vital Statistics. 
#' Validation of Exceptional Longevity. [Monographs on Population Aging. Vol. 6]
#' @export
#' 
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
#'Age <- 0:99
#' AHI(Value, Age, 90)
#' AHI(Value, Age, 95)
AHI <- function(Value, Age, Agei = 90){
	denomi <- Age %in% ((Agei - 2):(Agei + 2))
	denom  <- exp(sum(log(Value[denomi])) / 5)
	Value[Age == Agei] / denom
}


#' calculate Jdanov's old-age heaping index
#' 
#' @description This is a slightly more flexible implementation of Jdanov's fomula,
#' with defaults set to match his parameters. The numerator is the sum of 
#' (death counts) in ages 95, 100, and 105. Really the numerator could be an cound and any set of ages.
#' The denominator consists in the sum of the 5-year age groups centered around each of the numerator ages.
#' It probably only makes sense to use this with the default values, however. Used with a single age
#' in the numerator, it's almost the same as \code{Noumbissi()}, except here we pick out particular ages,
#' whereas Noumbissi picks out terminal digits.
#' 
#' @param Value a vector of demographic counts (probably deaths) by single age
#' @param Age a vector of ages corresponding to the lower integer bound of the counts
#' @param Ages a vector of ages to put in the numerator, default \code{c(95,100,105)}.
#' 
#' @return the index value
#' @export
#' 
#' @references 
#' Jdanov (2008) Beyond the Kannisto-Thatcher database on old age mortality - 
#' An assessment of data quality at advanced ages. MPIDR Working Paper WP-2008-013.
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
#'WI(Value, Age, c(95,100,105))

WI <- function(Value, Age, Ages = seq(95,105,by=5)){
	numi   <- Ages %in% Age
	
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