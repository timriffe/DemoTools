# Author: tim, based on earlier version by sean
# modified by TR 18 July, 2018
# TODO: see individual TODO notes for functions
# also need to add generic DemoTools wrapper here, will offer tail imputation options
# and call various methods with a method argument.
###############################################################################

#' The Carrier-Farrag method of population count smoothing
#' @description Smooth population counts in 5-year age groups.
#' @details This method does not account for ages < 10 nor for the 10 year age 
#' interval prior to the open age group. These are returned imputed with \code{NA}.
#' Age classes must be cleanly groupable to 5-year age groups.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' # from PASEX AGESMTH
#' MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 
#' 		198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
#' Ages         <- seq(0, 80, by = 5)
#' KKNtest <- c(NA,NA,354871,278502,285508,261429,236513 ,
#' 		204233,162138,126555,90094,65988,54803,41041,NA,NA,NA)
#' 
#' CFmales <- carrier_farrag_smth(MalePop, Ages, TRUE)
#' CFtest <- c(NA,NA,346290,287083,285855,261082,237937,
#' 202809,162973,125720,88730,67352,55187,40657,NA,NA,NA)
#' all(round(CFmales) - CFtest == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, MalePop)
#' lines(as.integer(names(CFmales)),CFmales)
#' }
#' @references
#' Carrier, Norman H., and A. M. Farrag. "The reduction of errors 
#' in census populations for statistically underdeveloped countries." 
#' Population Studies 12.3 (1959): 240-285.

# plz add PASEX citation, and add above to bibtex

carrier_farrag_smth <- function(Value, 
		Age, 
		OAG = TRUE){
	# these values are not used, it's just for lengths, and to make sure we 
	# end on an even 10. Technically we could even provide data in 10-year
	# age groups and it'd still not break.
	Value1     <- splitUniform(Counts = Value, Age = Age, OAG = OAG)
	Value5     <- groupAges(Value1, Age = as.integer(names(Value1)), N = 5)
	N          <- length(Value5)
	Age5       <- as.integer(names(Value5))
	
	# would need to move this up to ensure?
	# or in case of 85+ would we want to keep 80-84, 85+ as-is?
	Value10    <- groupAges(Value, Age = Age, N = 10)
	
	# what OAG is a strange digit? Then take OAG after grouping.
	if (OAG){
		OAGvalue <- Value10[length(Value10)]
		Value10[length(Value10)] <- NA
	}
	
	v10R       <- shift.vector(Value10, 1, fill = NA)
	v10L       <- shift.vector(Value10, -1, fill = NA)
	vevens     <- Value10 / (1 + (v10R / v10L) ^ .25)
	vodds      <- Value10 - vevens
	
	out        <- Value5 * NA	
	# cut back down (depending) and name
	interleaf  <- c(rbind(vodds, vevens))
	n          <- min(c(length(interleaf),N))
	out[1:n]   <- interleaf[1:n]

	out
	# tail behavior will be controlled in top level function.
}

#' Karup-King-Newton method of population count smoothing
#' @description Smooth population counts in 5-year age groups.
#' @details This method does not account for ages < 10 nor for the 10 year age 
#' interval prior to the open age group. These are returned imputed with \code{NA}.
#' Age classes must be cleanly groupable to 5-year age groups.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' # from PASEX AGESMTH
#' MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 
#' 		198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
#' Ages         <- seq(0, 80, by = 5)
#' KKNtest <- c(NA,NA,354871,278502,285508,261429,236513 ,
#' 		204233,162138,126555,90094,65988,54803,41041,NA,NA,NA)
#' 
#' KKNmales <- kkn_smth(MalePop, Ages, TRUE)
#' all(round(KKNmales) - KKNtest == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, MalePop)
#' lines(as.integer(names(KKNmales)),KKNmales)
#' }
#' @references
#' Carrier, Norman H., and A. M. Farrag. "The reduction of errors 
#' in census populations for statistically underdeveloped countries." 
#' Population Studies 12.3 (1959): 240-285.
# plz add PASEX citation, and add above to bibtex

kkn_smth <- function(Value, 
		Age, 
		OAG = TRUE){
	
	# these values are not used, it's just for lengths, and to make sure we 
	# end on an even 10. Technically we could even provide data in 10-year
	# age groups and it'd still not break.
	Value1     <- splitUniform(Counts = Value, Age = Age, OAG = OAG)
	Value5     <- groupAges(Value1, Age = as.integer(names(Value1)), N = 5)
	N          <- length(Value5)
	Age5       <- as.integer(names(Value5))
	
	# would need to move this up to ensure?
	# or in case of 85+ would we want to keep 80-84, 85+ as-is?
	Value10    <- groupAges(Value, Age = Age, N = 10)
	
	# what OAG is a strange digit? Then take OAG after grouping.
	if (OAG){
		OAGvalue <- Value10[length(Value10)]
		Value10[length(Value10)] <- NA
	}
	
	v10R       <- shift.vector(Value10, 1, fill = NA)
	v10L       <- shift.vector(Value10, -1, fill = NA)
	
	# this is the KNN operation
	vodds      <- Value10 / 2 + (v10R - v10L) / 16
	# constrained in 10-year age groups
	vevens     <- Value10 - vodds
	# stagger odd even 5s
	interleaf  <- c(rbind(vodds, vevens))
	# produce results vector
	out        <- Value5 * NA 
	n          <- min(c(length(interleaf),N))
    out[1:n]   <- interleaf[1:n]

	out
}

#' E. Arriaga's method of population count smoothing
#' @description Smooth population counts in 5-year age groups.
#' @details The open age group is aggregated down to be evenly divisible by 10. 
#' This method accounts for the youngest and oldest age groups. Age classes must be cleanly groupable to 5-year age groups.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' # from PASEX AGESMTH
#' MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 
#' 		198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
#' Ages         <- seq(0, 80, by = 5)
#' AMales       <- arriaga_smth(Value = MalePop, Age = Ages, OAG = TRUE)
#' # PAS spreadsheet result:
#' Atest        <- c(662761, 495126, 345744, 287629, 285919, 261018, 237469, 203277, 
#' 161733, 126960, 88586, 67496, 54587, 41257, 28790, 17189, 34729) 
#' all(round(AMales) - Atest == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, MalePop)
#' lines(as.integer(names(AMales)),AMales)
#' }
#' @references 
#' Arriaga, Eduardo, 1968. New Life Tables for Latin American Populations in the
#' Nineteenth and Twentieth Centuries, Population Monograph Series, no. 3,
#' appendix 3, page 295. University of California, Berkeley.

# plz add ref to bibtex, also cite PASEX and 
# DemGen_1994_Method_Arriaga....pdf page 40
arriaga_smth <- function(Value, 
		Age, 
		OAG = TRUE){
	# these values are not used, it's just for lengths, and to make sure we 
	# end on an even 10. Technically we could even provide data in 10-year
	# age groups and it'd still not break.
	Value1     <- splitUniform(Counts = Value, Age = Age, OAG = OAG)
	Value5     <- groupAges(Value1, Age = as.integer(names(Value1)), N = 5)
	N          <- length(Value5)
	Age5       <- as.integer(names(Value5))
	
	# would need to move this up to ensure?
	# or in case of 85+ would we want to keep 80-84, 85+ as-is?
	Value10    <- groupAges(Value, Age = Age, N = 10)
	
	# what OAG is a strange digit? Then take OAG after grouping.
	if (OAG){
		OAGvalue <- Value5[length(Value5)]
		Value10[length(Value10)] <- NA
		Value5[length(Value5)]   <- NA
	}
	
	# again get staggered vectors
	Value10LL   <- shift.vector(Value10, -2, fill = NA)
	Value10L    <- shift.vector(Value10, -1, fill = NA)
	Value10R    <- shift.vector(Value10, 1, fill = NA)
	Value10RR   <- shift.vector(Value10, 2, fill = NA)
	
	# alternating calc, with differences at tails
	vevens      <- (-Value10R + 11 * Value10 + 2 * Value10L) / 24
	# tails different
	vevens[1]   <- (8 * Value10[1] + 5 * Value10L[1] - Value10LL[1]) / 24
	lastind     <- which(is.na(vevens))[1]
	vevens[lastind] <-  Value10[lastind] - (8 * Value10[lastind] + 5 * Value10R[lastind] - Value10RR[lastind]) / 24
	# odds are complement
	vodds       <- Value10 - vevens
	
	# prepare output
	interleaf  <- c(rbind(vodds, vevens))
	# produce results vector
	out        <- Value5 * NA 
	n          <- min(c(length(interleaf),N))
	out[1:n]   <- interleaf[1:n]
	
	# if OA ends in 5, then we can save penultimate value too
	na.i       <- is.na(out)
	out[na.i]  <- Value5[na.i]
	if (OAG){
		out[N] <- OAGvalue
	}
	
	out
}

#' The old United Nations method of population count smoothing
#' @description Smooth population counts in 5-year age groups.
#' @details The open age group is aggregated down to be evenly divisible by 10. 
#' This method accounts for the youngest and oldest age groups. Age classes must be cleanly groupable to 5-year age groups.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' # from PASEX AGESMTH
#' MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 
#' 		198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
#' Ages         <- seq(0, 80, by = 5)
#' un_test <- c(NA,NA,364491,279123,268724,272228,243638,200923,162752,126304,
#' 		91662,67432,54677,38833,NA,NA,NA)
#' un_result <- united_nations_smth(MalePop, Ages, TRUE)
#' all(round(un_result) - un_test == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, MalePop)
#' lines(as.integer(names(un_result)),un_result)
#' }
#' @references
#' Carrier, Norman H., and A. M. Farrag. "The reduction of errors 
#' in census populations for statistically underdeveloped countries." 
#' Population Studies 12.3 (1959): 240-285.
# plz add Demgen 1994 and PASEX. I can't download above citation while 
# on vacation. plz check if this method is discussed there as well.
united_nations_smth <- function(Value, 
		Age, 
		OAG = TRUE){
	N <- length(Value)
	if (OAG){
		Value[N] <- NA
	}
	
	Value5     <- groupAges(Value, Age = Age, N = 5)
	N          <- length(Value5)
	Age5       <- as.integer(names(Value5))
	
	# get staggered vectors
	Value5LL   <- shift.vector(Value5, -2, fill = NA)
	Value5L    <- shift.vector(Value5, -1, fill = NA)
	Value5R    <- shift.vector(Value5, 1, fill = NA)
	Value5RR   <- shift.vector(Value5, 2, fill = NA)
	
	# this is the funny KNN operation
	# B11 is central
	out  <-  (-Value5RR + 4 * (Value5L + Value5R) + 10 * Value5 - Value5LL) / 16
	
	# cut back down (depending) and name
	names(out) <- Age5
	out
}

#' A strong method of population count smoothing
#' @description Smooth population counts in 5-year age groups.
#' @details The open age group is aggregated down to be evenly divisible by 10. 
#' This method accounts for the youngest and oldest age groups. Age classes must be cleanly 
#' groupable to 5-year age groups. All age classes are returned, but the strongest adjustment
#' ocurrs in ages bounded by \code{minA} and \code{maxA}. To be clear \code{maxA} refers to 
#' the lower bound of the highest age class, inclusive. So, if you want a ceiling of 70 (default), specify 65.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @param minA integer. The lowest age included included in intermediate adjustment. Default 10.
#' @param maxA integer. The highest age class included in intermediate adjustment. Default 65.
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' # from PASEX AGESMTH
#' MalePop      <- c(642367, 515520, 357831, 275542, 268336, 278601, 242515, 
#' 		198231, 165937, 122756, 96775, 59307, 63467, 32377, 29796, 16183, 34729)
#' Ages         <- seq(0, 80, by = 5)
#' strongtest <- c(646617,511270,386889,317345,273736,240058,218645,188297, 
#' 		153931, 124347,93254,71858,53594,39721,27887,18092,34729 ) 
#' strong_result <- strong_smth(MalePop,Ages,TRUE)
#' # differences due to intermediate rounding in spreadsheet (bad practice IMO)
#' all(abs(strong_result - strongtest) < 1, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, MalePop)
#' lines(as.integer(names(strong_result)),strong_result)
#' }
#' @references
#' Carrier, Norman H., and A. M. Farrag. "The reduction of errors 
#' in census populations for statistically underdeveloped countries." 
#' Population Studies 12.3 (1959): 240-285.
# plz add Demgen 1994 and PASEX. I can't download above citation while 
# on vacation. plz check if this method is discussed there as well.

strong_smth <- function(Value, 
		Age, 
		OAG = TRUE,
		minA = 10,
		maxA = 65){
	
	# these values are not used, it's just for lengths, and to make sure we 
	# end on an even 10. Technically we could even provide data in 10-year
	# age groups and it'd still not break.
	Value1     <- splitUniform(Counts = Value, Age = Age, OAG = OAG)
	Value5     <- groupAges(Value1, Age = as.integer(names(Value1)), N = 5)
	N          <- length(Value5)
	Age5       <- as.integer(names(Value5))
	
	# would need to move this up to ensure?
	# or in case of 85+ would we want to keep 80-84, 85+ as-is?
	Value10    <- groupAges(Value, Age = Age, N = 10)
	
	# what OAG is a strange digit? Then take OAG after grouping.
	if (OAG){
		OAGvalue <- Value5[length(Value5)]
		Value10[length(Value10)] <- NA
		Value5[length(Value5)]   <- NA
	}
	Age5       <- as.integer(names(Value5))
	Age10      <- as.integer(names(Value10))
	
	# subtotal
	indsub     <- Age10 >= minA & Age10 <= maxA
	SubTot     <- sum(Value10[indsub])
	
#	# get staggered vectors
	Value10L    <- shift.vector(Value10, -1, fill = NA)
	Value10R    <- shift.vector(Value10, 1, fill = NA)
	# this is the funny KNN operation
	# B11 is central
	Value10Pert <- (Value10 * 2 + Value10L + Value10R) / 4
	Value10Pert[is.na(Value10Pert)] <- Value10[is.na(Value10Pert)]
	
	# rescale ages between min and max to sum to original
	Value10Adj  <- Value10Pert
	Value10Adj[indsub] <- Value10Adj[indsub] * SubTot / sum(Value10Adj[indsub])
	
	# again get staggered vectors
	V10adjLL    <- shift.vector(Value10Adj, -2, fill = NA)
	V10adjL     <- shift.vector(Value10Adj, -1, fill = NA)
	V10adjR     <- shift.vector(Value10Adj, 1, fill = NA)
	V10adjRR    <- shift.vector(Value10Adj, 2, fill = NA)
	
	# alternating calc, with differences at tails
	vevens     <- (-V10adjR + 11 * Value10Adj + 2 * V10adjL) / 24
	# tails different
	vevens[1]  <- (8 * Value10Adj[1] + 5 * V10adjL[1] - V10adjLL[1]) / 24
	lastind     <- which(is.na(vevens))[1]
	vevens[lastind] <-  Value10Adj[lastind] - (8 * Value10Adj[lastind] + 5 * V10adjR[lastind] - V10adjRR[lastind]) / 24
	# odds are complement
	vodds      <- Value10Adj - vevens
	
	# prepare output
    # prepare output
	interleaf  <- c(rbind(vodds, vevens))
    # produce results vector
	out        <- Value5 * NA 
	n          <- min(c(length(interleaf),N))
	out[1:n]   <- interleaf[1:n]
	
    # if OA ends in 5, then we can save penultimate value too
	na.i       <- is.na(out)
	out[na.i]  <- Value5[na.i]
	if (OAG){
		out[N] <- OAGvalue
	}
	
	# what if OAis e.g. 85?
	out
}




