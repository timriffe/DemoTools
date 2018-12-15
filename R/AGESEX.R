
# Author: tim
###############################################################################

#' Calculate the PAS age ratio score
#' @description A single ratio is defined as the 100 times twice the size of an age group 
#' to its neighboring two age groups. Under uniformity these would all be 100. The average 
#' absolute deviation from 100 defines this index. This comes from the PAS spreadsheet called AGESEX

#' @param Value numeric. A vector of demographic counts in 5-year age groups.
#' @param Age numeric. A vector of ages corresponding to the lower integer bound of the counts.
#' @param ageMin integer. The lowest age included in calculations. Default 0.
#' @param ageMax integer. The upper age bound used for calculations. Default \code{max(Age)}.
#' @param method character. Either \code{"UN"} (default), \code{"Zelnick"}, or \code{"Ramachandran"}
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 

#' @details Age groups must be of equal intervals. Five year age groups are assumed.
#'  It is also assumed that the final age group is open, unless \code{ageMax < max(Age)}.
#'  Setting \code{OAG = FALSE} will override this and potentially include \code{max(Age)} in calculations.
#'  The method argument determines the weighting of numerators and denominators, where the UN method puts
#' twice the numerator over the sum of the adjacent ages classes, Zelnich does thrice the
#' numerator over the sum of the whole range from the next lowest to the next highest age, and
#' Ramachandran does four times the numerator over the same sum, but with the central age
#' double-counted in the numerator.

#' @return The value of the index.
#' 
#' @references 
#' \insertRef{accuracyun1952}{DemoTools}
#' \insertRef{ramachandran1967}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#' @export
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

#' Calculate the PAS sex ratio score
#' @description A single ratio is defined as males per 100 females.
#' Taking the first difference of the ratios over age would give 0s
#' if the ratio were constant. The average absolute difference over
#' age defines the index. This comes from the PAS spreadsheet called AGESEX.

#' @param Males numeric.  A vector of demographic counts in 5-year age groups for males.
#' @param Females numeric. A vector of demographic counts in 5-year age groups for females.
#' @param Age numeric. A vector of ages corresponding to the lower integer bound of the counts.
#' @param ageMin integer. The lowest age included in calculations. Default 0.
#' @param ageMax integer. The upper age bound used for calculations. Default \code{max(Age)}.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @details Age groups must be of equal intervals. Five year age groups are assumed.
#'  It is also assumed that the final age group is open, unless \code{ageMax < max(Age)}. 

#' @return The value of the index.
#' @references 
#' \insertRef{accuracyun1952}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#' @export
#' @examples 
#' Males   <- c(4677000,4135000,3825000,3647000,3247000,2802000,2409000,2212000,
#' 		1786000,1505000,1390000,984000,745000,537000,346000,335000)
#' Females <- c(4544000,4042000,3735000,3647000,3309000,2793000,2353000,2112000,
#' 		1691000,1409000,1241000,887000,697000,525000,348000,366000)
#' Age     <- seq(0, 75, by = 5)
#' sexRatioScore(Males, Females, Age)  # 2.2, matches PAS


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


#' Calculate an age-sex accuracy index
#' @description This index is a composite consisting in the sum of thrice the sex 
#' ratio index plus the age ratio index for males and females. This function is therefore
#' a wrapper to \code{ageRatioScore()} and \code{sexRatioScore()}.

#' @param Males numeric.  A vector of demographic counts in 5-year age groups for males.
#' @param Females numeric. A vector of demographic counts in 5-year age groups for females.
#' @param Age numeric. A vector of ages corresponding to the lower integer bound of the counts.
#' @param ageMin integer. The lowest age included in calculations. Default 0.
#' @param ageMax integer. The upper age bound used for calculations. Default \code{max(Age)}.
#' @param method character. Either \code{"UN"} (default), \code{"Zelnick"}, or \code{"Ramachandran"},or \code{"das gupta"}.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @param adjust logical. Whether or not to adjust the measure when population size is under one million. Default \code{TRUE}.

#' @details Age groups must be of equal intervals. Five year age groups are assumed.
#' If the final element of \code{Males} and \code{Females} is the open age group,
#' then either make sure \code{ageMax} is lower than it, or leave \code{OAG} as \code{TRUE} so that it is properly removed for calculations. 
#' The method argument is passed to \code{ageRatioScore()}, where it determines weightings of numerators and denominators, 
#' except in the case of Das Gupta, where it's a different method entirely (see \code{ageSexAccuracyDasGupta()}.

#' @return The value of the index.
#' @references 
#' \insertRef{accuracyun1952}{DemoTools}
#' \insertRef{dasgupta1955}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#' @export
#' @examples
#' Males   <- c(4677000,4135000,3825000,3647000,3247000,2802000,2409000,2212000,
#' 		1786000,1505000,1390000,984000,745000,537000,346000,335000)
#' Females <- c(4544000,4042000,3735000,3647000,3309000,2793000,2353000,2112000,
#' 		1691000,1409000,1241000,887000,697000,525000,348000,366000)
#' Age     <- seq(0, 75, by = 5)
#' ageSexAccuracy(Males, Females, Age)    # 14.3, matches PAS
#' ageSexAccuracy(Males, Females, Age, ageMax = 75)
#' ageSexAccuracy(Males, Females, Age, ageMax = 75, adjust = FALSE)
#' ageSexAccuracy(Males, Females, Age, method = "Zelnick")
#' ageSexAccuracy(Males, Females, Age, method = "Ramachandran")
#' # Das Gupta not a comparable magnitude, FYI.
#' ageSexAccuracy(Males, Females, Age, method = "Das Gupta")
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

#' Calculate Das Gupta's (1995) age sex accuracy index
#' @description Given population counts in 5-year age groups for males and females, follow Das Gupta's steps
#'  to calculate a composite index of the quality of the age and sex structure for a given population.
#' 
#' @param Males numeric.  A vector of demographic counts in 5-year age groups for males.
#' @param Females numeric. A vector of demographic counts in 5-year age groups for females.
#' @param Age numeric. A vector of ages corresponding to the lower integer bound of the counts.
#' @param ageMin integer. The lowest age included in calculations. Default 0.
#' @param ageMax integer. The upper age bound used for calculations. Default \code{max(Age)}.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' 
#' @references 
#' \insertRef{dasgupta1955}{DemoTools}
#' 
#' @details It is assumed that the terminal age group is open, in which case it is ignored.  
#' Set \code{OAG = FALSE} if the top age is indeed a closed interval that you want included in calculations. 
#' If \code{ageMax == max(Age)} and \code{OAG} is \code{TRUE}, then \code{ageMax} gets decremented one age class.
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
