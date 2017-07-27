
# Author: tim
###############################################################################

#' Spoorenbergs total modified Whipple index

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
