
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
#' if it's 0,5. Otherwise just one digit at a time. \code{ageMin} and \code{ageMax} are applied to numerator ages, not denominators.
#'  Denominators are always 5-year age groups centered on the digit in question,
#' and these therefore stretch into ages a bit higher or lower than the numerator ages.
#' @return the index value 
#' @export 

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

#Whipple(Value, Age, 25, 60, digit = c(0,5)) # 2.34,replicates SINGAGE males

# implements formula from GDA_1981_Structures-par-age-et-sexe-en-Afrique_[IREDA]
# p 148
#Whipple(Value, Age, 25, 60, digit = 0)
#Whipple(Value, Age, 25, 60, digit = 5) 

# Whipple types
#Whipple(Value, Age, 25, 60, digit = 3) 
