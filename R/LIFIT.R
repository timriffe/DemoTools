
# Author: tim
###############################################################################

# function design based on Computer programs for demographic analysis 
# Eduardo Arriaga, Patricia Anderson, Larry Heligman.
# Washington : U.S. Dept. of Commerce, Bureau of the Census, International Statistical Programs
# Center, 1976.

#' Mean absolute difference in survival rates.
#' 
#' @description The mean absolute difference between 5-year survival 
#' ratios with optional age trimming. 
#' Based on LIFIT program from the MPCDA code base.
#' 
#' @param lx1 numeric. Vector of the first survival function, any radix size.
#' @param lx2 numeric. Vector of the first survival function, any radix size.
#' @param Age1 integer. A vector of ages of the lower integer bound of the age groups corresponding to \code{lx1}.
#' @param Age2 integer. A vector of ages of the lower integer bound of the age groups corresponding to \code{lx2}.
#' @param ageMin integer. Optional lower age bound for calculations.
#' @param ageMax integer. Optional upper age bound for calculations.
#' @export
#' @details \code{lx1} and \code{lx2} can be of different lengths and grouped differently. 
#' For example, either or both input vectors could be in single or five-year age groups
#' for this function. Results are selected internally for 5-year age groups, and vector lengths are matched.
#' \code{ageMax} is an inclusive upper bound, treated as interval.
#'  If you want the age-group 70-74 to be included in calculations, then give \code{ageMax = 70}, not 75.
#' The input vectors can also be of different lifetable radices. As with the original Fortran code, 
#' the open age group must be included in input vectors, even though it is not included in calcs. 
#' If your final age group is not open, then tack an extra element to it to simulate an open element,
#' and just make sure to match ages.
#' @return The mean absolute difference in survival ratios.
#' @references
#' \insertRef{arriaga1976}{DemoTools}
#' 
#' @examples 
#' lx1 <-c(1, 0.94265, 0.93344, 0.92708, 0.923, 0.92023, 0.91796, 0.91601, 
#'		0.91431, 0.9128, 0.91145, 0.91025, 0.90914, 0.90805, 0.90692, 
#'		0.90566, 0.90419, 0.90251, 0.90063, 0.89858, 0.89642, 0.89417, 
#'		0.8918, 0.88923, 0.88647, 0.88356, 0.88046, 0.87721, 0.87386, 
#'		0.87047, 0.86706, 0.8636, 0.86008, 0.85649, 0.85276, 0.84882, 
#'		0.84465, 0.84027, 0.8357, 0.83102, 0.8264, 0.82193, 0.8175, 0.81296, 
#'		0.80827, 0.80343, 0.79845, 0.79329, 0.78788, 0.78218, 0.77617, 
#'		0.76982, 0.76309, 0.75596, 0.74837, 0.74029, 0.73182, 0.72307, 
#'		0.71396, 0.70434, 0.69402, 0.68295, 0.67102, 0.6581, 0.64442, 
#'		0.63012, 0.61486, 0.59811, 0.57943, 0.55896, 0.53739, 0.51528, 
#'		0.49248, 0.46826, 0.44272, 0.41631, 0.38876, 0.36037, 0.33157, 
#'		0.30217, 0.27334, 0.24599, 0.21948, 0.19412, 0.17029, 0.14826, 
#'		0.12771, 0.10909, 0.09241, 0.07763, 0.0647, 0.05363, 0.04403, 
#'		0.03578, 0.02877, 0.02289, 0.018, 0.01399, 0.01075, 0.00815, 
#'		0.0061, 0.00451, 0.00329, 0.00236, 0.00168, 0.00117, 0.00081, 
#'		0.00055, 0.00036, 0.00024, 0.00015)
#'lx2 <- c(1, 0.99361, 0.99315, 0.99284, 0.99264, 0.99246, 0.9923, 0.99216, 
#'		0.99201, 0.99189, 0.99176, 0.99164, 0.99149, 0.99135, 0.99118, 
#'		0.99096, 0.99069, 0.99033, 0.98993, 0.98946, 0.98902, 0.98857, 
#'		0.98809, 0.98763, 0.98717, 0.98672, 0.98623, 0.98574, 0.98521, 
#'		0.98468, 0.9841, 0.98351, 0.98285, 0.98213, 0.98134, 0.9805, 
#'		0.97959, 0.9786, 0.97748, 0.97628, 0.97495, 0.97353, 0.97197, 
#'		0.97032, 0.96847, 0.96655, 0.96443, 0.96225, 0.95984, 0.95724, 
#'		0.95446, 0.95142, 0.94822, 0.94469, 0.94058, 0.93651, 0.93187, 
#'		0.92678, 0.92112, 0.91508, 0.90846, 0.90128, 0.89336, 0.88466, 
#'		0.87543, 0.86532, 0.85442, 0.84265, 0.82993, 0.81641, 0.80212, 
#'		0.78656, 0.77024, 0.75237, 0.73333, 0.71304, 0.69102, 0.66797, 
#'		0.64371, 0.61749, 0.58998, 0.55953, 0.52836, 0.49505, 0.46051, 
#'		0.42472, 0.38756, 0.34995, 0.31203, 0.27485, 0.23825, 0.20292, 
#'		0.16979, 0.13867, 0.11131, 0.08713, 0.06679, 0.04989, 0.03626, 
#'		0.02561, 0.01755, 0.01165, 0.00749, 0.00466, 0.00281, 0.00163, 
#'		0.00092, 5e-04, 0.00026, 0.00013, 7e-05)
#' ADM(lx1,lx2, 0:110)
#' ADM(lx1,lx2,0:110,ageMax = 80)
ADM <- function(
		lx1, 
		lx2, 
		Age1, 
		Age2 = Age1, 
		ageMin = max(c(min(Age1),min(Age2))), 
		ageMax = min(c(max(Age1), max(Age2)))){
	
    stopifnot(length(lx1) == length(Age1))
    stopifnot(length(lx2) == length(Age2))
	
	# TODO add OAG to control removal
	
	# remove open age group
	n1     <- length(lx1)
	n2     <- length(lx2)
	lx1    <- lx1[-n1]
	lx2    <- lx2[-n2]
	Age1   <- Age1[-n1]
	Age2   <- Age2[-n2]
	
	# ages in common only
	age5   <- AGEN(
			    Age1, 
			    Age2, 
			    N = 5, 
			    consecutive = TRUE, 
			    ageMin = ageMin, 
			    ageMax = ageMax)
	# trim lx to needed ages
	lx1_5  <- lx1[Age1 %in% age5]
	lx2_5  <- lx2[Age2 %in% age5]
	
	# survivor ratios
	Sx1    <- ratx(lx1_5)
	Sx2    <- ratx(lx2_5)
	
	# absolute differences
	AD     <- Sx1 - Sx2
	
	# mean absolute difference
	mean(abs(AD))
}

#' Mean absolute difference in age-ratios of survival rates.
#' 
#' @description Take two survival curves, use them to calculate 5-year survival 
#' ratios, then again take the consecutive age-ratios of these. Finally, take the 
#' mean absolute difference between ratios, with optional age trimming. 
#' Based on LIFIT program from the MPCDA code base.
#' 
#' @param lx1 numeric. Vector of the first survival function, any radix size.
#' @param lx2 numeric. Vector of the first survival function, any radix size.
#' @param Age1 integer. A vector of ages of the lower integer bound of the age groups corresponding to \code{lx1}.
#' @param Age2 integer. A vector of ages of the lower integer bound of the age groups corresponding to \code{lx2}.
#' @param ageMin integer. Optional lower age bound for calculations.
#' @param ageMax integer. Optional upper age bound for calculations.
#' @export
#' @details \code{lx1} and \code{lx2} can be of different lengths and grouped differently. 
#' For example, either or both input vectors could be in single or five-year age groups
#' for this function. Results are selected internally for 5-year age groups, and vector lengths are matched.
#' If you want the age-group 70-74 to be included in calculations, then give \code{ageMax = 70}, not 75.
#' The input vectors can also be of different lifetable radices. As with the original Fortran code, 
#' the open age group must be included in input vectors, even though it is not included in calcs. 
#' If your final age group is not open, then tack an extra element to it to simulate an open element,
#' and just make sure to match ages.
#' @return The mean absolute difference in survival ratios.
#' @references
#' \insertRef{arriaga1976}{DemoTools}
#' 
#' @examples 
#' lx1 <-c(1, 0.94265, 0.93344, 0.92708, 0.923, 0.92023, 0.91796, 0.91601, 
#'		0.91431, 0.9128, 0.91145, 0.91025, 0.90914, 0.90805, 0.90692, 
#'		0.90566, 0.90419, 0.90251, 0.90063, 0.89858, 0.89642, 0.89417, 
#'		0.8918, 0.88923, 0.88647, 0.88356, 0.88046, 0.87721, 0.87386, 
#'		0.87047, 0.86706, 0.8636, 0.86008, 0.85649, 0.85276, 0.84882, 
#'		0.84465, 0.84027, 0.8357, 0.83102, 0.8264, 0.82193, 0.8175, 0.81296, 
#'		0.80827, 0.80343, 0.79845, 0.79329, 0.78788, 0.78218, 0.77617, 
#'		0.76982, 0.76309, 0.75596, 0.74837, 0.74029, 0.73182, 0.72307, 
#'		0.71396, 0.70434, 0.69402, 0.68295, 0.67102, 0.6581, 0.64442, 
#'		0.63012, 0.61486, 0.59811, 0.57943, 0.55896, 0.53739, 0.51528, 
#'		0.49248, 0.46826, 0.44272, 0.41631, 0.38876, 0.36037, 0.33157, 
#'		0.30217, 0.27334, 0.24599, 0.21948, 0.19412, 0.17029, 0.14826, 
#'		0.12771, 0.10909, 0.09241, 0.07763, 0.0647, 0.05363, 0.04403, 
#'		0.03578, 0.02877, 0.02289, 0.018, 0.01399, 0.01075, 0.00815, 
#'		0.0061, 0.00451, 0.00329, 0.00236, 0.00168, 0.00117, 0.00081, 
#'		0.00055, 0.00036, 0.00024, 0.00015)
#'lx2 <- c(1, 0.99361, 0.99315, 0.99284, 0.99264, 0.99246, 0.9923, 0.99216, 
#'		0.99201, 0.99189, 0.99176, 0.99164, 0.99149, 0.99135, 0.99118, 
#'		0.99096, 0.99069, 0.99033, 0.98993, 0.98946, 0.98902, 0.98857, 
#'		0.98809, 0.98763, 0.98717, 0.98672, 0.98623, 0.98574, 0.98521, 
#'		0.98468, 0.9841, 0.98351, 0.98285, 0.98213, 0.98134, 0.9805, 
#'		0.97959, 0.9786, 0.97748, 0.97628, 0.97495, 0.97353, 0.97197, 
#'		0.97032, 0.96847, 0.96655, 0.96443, 0.96225, 0.95984, 0.95724, 
#'		0.95446, 0.95142, 0.94822, 0.94469, 0.94058, 0.93651, 0.93187, 
#'		0.92678, 0.92112, 0.91508, 0.90846, 0.90128, 0.89336, 0.88466, 
#'		0.87543, 0.86532, 0.85442, 0.84265, 0.82993, 0.81641, 0.80212, 
#'		0.78656, 0.77024, 0.75237, 0.73333, 0.71304, 0.69102, 0.66797, 
#'		0.64371, 0.61749, 0.58998, 0.55953, 0.52836, 0.49505, 0.46051, 
#'		0.42472, 0.38756, 0.34995, 0.31203, 0.27485, 0.23825, 0.20292, 
#'		0.16979, 0.13867, 0.11131, 0.08713, 0.06679, 0.04989, 0.03626, 
#'		0.02561, 0.01755, 0.01165, 0.00749, 0.00466, 0.00281, 0.00163, 
#'		0.00092, 5e-04, 0.00026, 0.00013, 7e-05)
#' RDM(lx1, lx2, 0:110)
#' RDM(lx1, lx2, 0:110, ageMax = 80)
RDM <- function(
		lx1, 
		lx2, 
		Age1, 
		Age2 = Age1, 
		ageMin = max(c(min(Age1),min(Age2))), 
		ageMax = min(c(max(Age1), max(Age2)))){
	stopifnot(length(lx1) == length(Age1))
	stopifnot(length(lx2) == length(Age2))
	
	# TODO add OAG to control removal
	
	# remove open age group
	n1     <- length(lx1)
	n2     <- length(lx2)
	lx1    <- lx1[-n1]
	lx2    <- lx2[-n2]
	Age1   <- Age1[-n1]
	Age2   <- Age2[-n2]
	
	# ages in common only
	age5   <- AGEN(
				Age1, 
				Age2, 
				N = 5, 
				consecutive = TRUE, 
				ageMin = ageMin, 
				ageMax = ageMax)
	# select only needed lx
	lx1_5  <- lx1[Age1 %in% age5]
	lx2_5  <- lx2[Age2 %in% age5]
	
	# survivor ratios
	Sx1    <- ratx(lx1_5)
	Sx2    <- ratx(lx2_5)
	
	# ratios of the ratios
	Rx1    <- ratx(Sx1)
	Rx2    <- ratx(Sx2)
	
	# differences of ratio ratios
	RD     <- Rx1 - Rx2
	
	# mean absolute difference
	mean(abs(RD))
}

