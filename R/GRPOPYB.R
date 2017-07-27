# Author: sean
###############################################################################
#
#Males      <- c(18052, 16523, 15555, 12063, 9654, 9256, 10058, 10035, 9785, 9025,
#7562, 6241, 5362, 4256, 3168, 1956, 1056, 568)
#Females    <- c(17274, 15915, 14455, 11475, 9202, 9226, 10162, 10633, 9962, 9226,
#7966, 6739, 5458, 4575, 3355, 2156, 1121, 636)
#Age        <- seq(0, 85, by = 5)
#CensusDate       <- 1950.25

#' create the historical birth cohorts for a census
#' @description The birth cohorts for a census represent the population born during a set interval
#' of time. The date assigned to the time interval is the midpoint of the interval.
#' This is generalized from the PAS spreadsheet called GRPOP-YB.

#' @param Value numeric vector of demographic population counts.
#' @param Age integer vector of ages corresponding to the lower integer bound of the age range.
#' @param CensusDate the date of the first census. See details for ways to express it.
#' @param cohortSize integer. The length of time (years) surrounding each output birth cohort. Default 5.

#' @details Age groups must be of equal intervals. No specific age structure is assumed for the census. Births
#' are assumed to happen uniformly over the cohorts' intervals. The final age group is assumed to be the same
#' size as all the other age groups. If the cohortSize does not divide evenly into the largest age in the data,
#' any additional (higher) ages needed are set as zero. For example, if cohortSize is 7 and the largest age is 90,
#' then age 91 (necessary for the matrix sum) is set as zero. Dates can be given in three ways 1) a 
#' \code{Date} class object, 2) an unambiguous character string in the format \code{"YYYY-MM-DD"}, or 
#' 3) as a decimal date consisting in the year plus the fraction of the year passed as of the given date. 

#' @return a data frame with a decimal date corresponding to the midpoint of the birth cohort and resulting population counts.
#' @export
#' 
#' @examples 
#' Males      <- c(18052, 16523, 15555, 12063, 9654, 9256, 10058, 10035, 9785, 9025,
#'                 7562, 6241, 5362, 4256, 3168, 1956, 1056, 568)
#' Females    <- c(17274, 15915, 14455, 11475, 9202, 9226, 10162, 10633, 9962, 9226,
#'                 7966, 6739, 5458, 4575, 3355, 2156, 1121, 636)
#' Age        <- seq(0, 85, by = 5)
#' CensusDate       <- 1950.25
#' birthCohorts(Males, Age, CensusDate)
#' birthCohorts(Females, Age, CensusDate)
#' birthCohorts(Females, Age, CensusDate, cohortSize = 10)

birthCohorts <- function(Value, Age, CensusDate, cohortSize = 5){
      
	  # TR: added date handling. dec.date() is in utils.R
	  CensusDate   <- dec.date(CensusDate)
	
      ageMax <- max(Age)             # the lower bound of the largest age group
      N      <- length(Value)        # number of age groups from the census.
      M      <- ageMax / (N - 1)     # length of each age group from the census.
      
      ageGroupBirths   <- Value / M  # vector of the number of births in a single year 
	                                 # for each age group assuming uniformity.
      
      singleAgeGroupBirths  <- rep(ageGroupBirths, each = M)  # vector of the single year 
	                                                          # births for all ages
      
      # Check that the cohort divides into the max age. If not, 
	  # add some zeros to prevent errors when summing across the vector.
      if (length(singleAgeGroupBirths) %% cohortSize != 0){
                  singleAgeGroupBirths <- c(singleAgeGroupBirths, rep(0, each = cohortSize - length(singleAgeGroupBirths)%%cohortSize))
      }
      
      # TR: replaced data.frame code with two vectors
	  CohBirths       <- colSums(matrix(singleAgeGroupBirths, nrow = cohortSize))   
	  yrs             <- rev(seq( 
			                 CensusDate - length(singleAgeGroupBirths) + cohortSize / 2, 
							 CensusDate,
							 by = cohortSize))
      # TR: now create a data.frame from scratch without coerce
	  # is this just reclassifying pops in the census by birth year?
	  # then why call it Births?
	  outputDataFrame <- data.frame(Year = yrs, Births = CohBirths)
      
      return(outputDataFrame)
}
