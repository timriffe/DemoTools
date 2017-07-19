
# Author: tim
###############################################################################
#
#Males   <- c(4677000,4135000,3825000,3647000,3247000,2802000,2409000,2212000,
#1786000,1505000,1390000,984000,745000,537000,346000,335000)
#Females <- c(4544000,4042000,3735000,3647000,3309000,2793000,2353000,2112000,
#1691000,1409000,1241000,887000,697000,525000,348000,366000)
#Age     <- seq(0, 75, by = 5)


#' calculate the PASEX age ratio score
#' @description A single ratio is defined as the 100 times twice the size of an age group 
#' to it's neighboring two age groups. Under uniformity these would all be 100. The average 
#' absolute deviation from 100 defines the index. This comes from the PASEX spreadsheet called AGESEX

#' @param Value numeric. A vector of demographic counts by single age.
#' @param Age an integer vector of ages corresponding to the lower integer bound of the counts.
#' @param ageMin integer. The lowest age included in calcs. Default 0.
#' @param ageMax integer. The upper age bound used for calcs. Default \code{max(Age)}.

#' @details Age groups must be of equal intervals. Five year age groups are assumed.
#'  We also assume that the final age group is open, unless \code{ageMax < max(Age)}

#' @return the value of the index
#' @export
#' 
#' @examples 
#' Males   <- c(4677000,4135000,3825000,3647000,3247000,2802000,2409000,2212000,
#' 		1786000,1505000,1390000,984000,745000,537000,346000,335000)
#' Females <- c(4544000,4042000,3735000,3647000,3309000,2793000,2353000,2112000,
#' 		1691000,1409000,1241000,887000,697000,525000,348000,366000)
#' Age     <- seq(0, 75, by = 5)
#' ageRatioScore(Males, Age)    # 3.9, matches pasex
#' ageRatioScore(Females, Age)  # 3.65 matches pasex
ageRatioScore <- function(Value, Age, ageMin = 0, ageMax = max(Age)){
	stopifnot(length(Value) == length(Age))
	
	N           <- length(Value)
	
	# remove open age group (assume final age open)
	if (ageMax == max(Age)){
		Value       <- Value[-N]
		Age         <- Age[-N]
		N           <- N - 1
		ageMax      <- max(Age)
	}
	
	numi        <- Age > ageMin & Age < ageMax
	denomi1     <- shift.vector(numi, -1) 
	denomi2     <- shift.vector(numi, 1)
	
	numerator   <- 200 * Value[numi]
	denominator <- Value[denomi1] + Value[denomi2] 
	# ratio of each age to the age above and below it
	ageratio    <- numerator / denominator
	# absolute deviatios from uniformity
	absres      <- abs(ageratio - 100)
	# average absolute deviation
	sum(absres) / length(absres)
}

#' calculate the PASEX sex ratio score
#' @description A single ratio is defined as males per 100 females.
#' Taking the first difference of the ratios over age would give 0s
#' if the ratio were constant. The average absolute difference over
#' age defines the index. This comes from the PASEX spreadsheet called AGESEX.

#' @param Males numeric. A vector of population counts for males by age
#' @param Females numeric. A vector of population counts for females by age
#' @param Age integer. A vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin integer. The lowest age included in calcs. Default 0
#' @param ageMax integer. The upper age bound used for calcs. Default \code{max(Age)}

#' @details Age groups must be of equal intervals. Five year age groups are assumed.
#'  We also assume that the final age group is open, unless \code{ageMax < max(Age)}

#' @return the value of the index
#' @export
#' @examples 
#' Males   <- c(4677000,4135000,3825000,3647000,3247000,2802000,2409000,2212000,
#' 		1786000,1505000,1390000,984000,745000,537000,346000,335000)
#' Females <- c(4544000,4042000,3735000,3647000,3309000,2793000,2353000,2112000,
#' 		1691000,1409000,1241000,887000,697000,525000,348000,366000)
#' Age     <- seq(0, 75, by = 5)
#' sexRatioScore(Males, Females, Age)  # 2.2, matches pasex

sexRatioScore <- function(Males, Females, Age, ageMin = 0, ageMax = max(Age)){
	stopifnot(length(Males) == length(Age) & length(Males) == length(Females))
	N         <- length(Males)
	
	
	ratio     <- Males / Females
	ratiodiff <- 100 * diff(ratio)
	absdiff   <- abs(ratiodiff)
	absdiff   <- absdiff[-length(absdiff)]
	sum(absdiff) / length(absdiff)
}


#' calculate the PASEX age-sex accuracy index
#' @description This index is a composite consisting in the sum of thrice the sex 
#' ratio index plus the age ratio index for males and females. This function is therefore
#' a wrapper to \code{ageRatioScore()} and \code{sexRatioScore()}.

#' @param Males numeric. A vector of population counts for males by age
#' @param Females numeric. A vector of population counts for females by age
#' @param Age integer. A vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin integer. The lowest age included in calcs. Default 0
#' @param ageMax integer. The upper age bound used for calcs. Default \code{max(Age)}

#' @details Age groups must be of equal intervals. Five year age groups are assumed.
#'  We also assume that the final age group is open, unless \code{ageMax < max(Age)}

#' @return the value of the index
#' @export

#' @examples 
#' Males   <- c(4677000,4135000,3825000,3647000,3247000,2802000,2409000,2212000,
#' 		1786000,1505000,1390000,984000,745000,537000,346000,335000)
#' Females <- c(4544000,4042000,3735000,3647000,3309000,2793000,2353000,2112000,
#' 		1691000,1409000,1241000,887000,697000,525000,348000,366000)
#' Age     <- seq(0, 75, by = 5)
#' ageSexAccuracy(Males, Females, Age)    # 14.3, matches pasex

ageSexAccuracy <- function(Males, Females, Age, ageMin = 0, ageMax = max(Age), adjust = TRUE){
	SR <- sexRatioScore(
			Males = Males, 
			Females = Females, 
			Age = Age, 
			ageMin = ageMin, 
			ageMax = ageMax)
	MA <- ageRatioScore(
			Value = Males,
			Age = Age,
			ageMin = ageMin, 
			ageMax = ageMax)
	FA <- ageRatioScore(
			Value = Females,
			Age = Age,
			ageMin = ageMin, 
			ageMax = ageMax)
	# calculate index:
	index <- 3 * SR + MA + FA
	
	# deflate if population small (formula in PAS)
	if (adjust){
		TOT   <- sum(Males) + sum(Females)
		if (TOT < 1000000){
		   index <- index - 3500 / sqrt(TOT) + 3.5
	    }
	}
	
	return(index)
}

