
# Author: tim
###############################################################################


#' calculate Myer's blended index of age heaping

#' @description Implementation following the PASEX spreadsheet SINGAGE

#' @param Value a vector of demographic counts by single age
#' @param Age a vector of ages corresponding to the lower integer bound of the counts
#' @param ageMin the lowest age included in calcs. Default 10
#' @param ageMax the upper age bound used for calcs. Default 90

#' @details \code{ageMax} is the hard upper bound, treated as interval. If you want ages
#' 20 to 89, then give \code{ageMin = 20} and \code{ageMax = 90}, not 89.

#' @export 

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


#Value <- c(80626,95823,104315,115813,100796,105086,97266,116328,
#75984,89982,95525,56973,78767,65672,53438,85014,
#47600,64363,42195,42262,73221,30080,34391,29072,
#20531,66171,24029,44227,24128,23599,82088,16454,
#22628,17108,12531,57325,17220,28425,16206,17532,
#65976,11593,15828,13541,8133,44696,11165,18543,
#12614,12041,55798,9324,10772,10453,6773,28358,
#9916,13348,8039,7583,42470,5288,5317,6582,
#3361,17949,3650,5873,3279,3336,27368,1965,
#2495,2319,1335,12022,1401,1668,1360,1185,
#9167,424,568,462,282,6206,343,409,333,291,4137,133,169,157,89,2068,68,81,66,57)
#Age <- 0:99

#Myers(Value, Age, 10, 90) * 2 #47.46, replicates SINGAGE males
