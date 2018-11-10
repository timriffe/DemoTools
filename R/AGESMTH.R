# Author: tim, based on earlier version by sean
# modified by TR 18 July, 2018
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
#' Ages         <- seq(0, 80, by = 5)
#' KKNtest <- c(NA,NA,354871,278502,285508,261429,236513 ,
#' 		204233,162138,126555,90094,65988,54803,41041,NA,NA,NA)
#' 
#' CFmales <- carrier_farrag_smth(pop5m_pasex, Ages, TRUE)
#' CFtest <- c(NA,NA,346290,287083,285855,261082,237937,
#' 202809,162973,125720,88730,67352,55187,40657,NA,NA,NA)
#' all(round(CFmales) - CFtest == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, pop5m_pasex)
#' lines(as.integer(names(CFmales)),CFmales)
#' }
#' @references
#' \insertRef{carrier1959reduction}{DemoTools}

# plz add PASEX citation, and add above to bibtex

carrier_farrag_smth <- function(Value, 
		Age, 
		OAG = TRUE){
	# these values are not used, it's just for lengths, and to make sure we 
	# end on an even 10. Technically we could even provide data in 10-year
	# age groups and it'd still not break.
	Value1     <- splitUniform(Value = Value, Age = Age, OAG = OAG)
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
#' Ages         <- seq(0, 80, by = 5)
#' KKNtest <- c(NA,NA,354871,278502,285508,261429,236513 ,
#' 		204233,162138,126555,90094,65988,54803,41041,NA,NA,NA)
#' 
#' KKNmales <- kkn_smth(pop5m_pasex, Ages, TRUE)
#' all(round(KKNmales) - KKNtest == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, pop5m_pasex)
#' lines(as.integer(names(KKNmales)),KKNmales)
#' }
#' @references
#' \insertRef{carrier1959reduction}{DemoTools}
# plz add PASEX citation, and add above to bibtex

kkn_smth <- function(Value, 
		Age, 
		OAG = TRUE){
	
	# these values are not used, it's just for lengths, and to make sure we 
	# end on an even 10. Technically we could even provide data in 10-year
	# age groups and it'd still not break.
	Value1     <- splitUniform(Value = Value, Age = Age, OAG = OAG)
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
#' Ages         <- seq(0, 80, by = 5)
#' AMales       <- arriaga_smth(Value = pop5m_pasex, Age = Ages, OAG = TRUE)
#' # PAS spreadsheet result:
#' Atest        <- c(662761, 495126, 345744, 287629, 285919, 261018, 237469, 203277, 
#' 161733, 126960, 88586, 67496, 54587, 41257, 28790, 17189, 34729) 
#' all(round(AMales) - Atest == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, pop5m_pasex)
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
	Value1     <- splitUniform(Value = Value, Age = Age, OAG = OAG)
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
#' Ages         <- seq(0, 80, by = 5)
#' un_test <- c(NA,NA,364491,279123,268724,272228,243638,200923,162752,126304,
#' 		91662,67432,54677,38833,NA,NA,NA)
#' un_result <- united_nations_smth(pop5m_pasex, Ages, TRUE)
#' all(round(un_result) - un_test == 0, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, pop5m_pasex)
#' lines(as.integer(names(un_result)),un_result)
#' }
#' @references
#' \insertRef{carrier1959reduction}{DemoTools}
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
#' occurs in ages bounded by \code{ageMin} and \code{ageMax}. To be clear \code{ageMax} refers to 
#' the lower bound of the highest age class, inclusive. So, if you want a ceiling of 70 (default), specify 65.
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age numeric vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @param ageMin integer. The lowest age included included in intermediate adjustment. Default 10.
#' @param ageMax integer. The highest age class included in intermediate adjustment. Default 65.
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' @examples
#' Ages         <- seq(0, 80, by = 5)
#' strongtest <- c(646617,511270,386889,317345,273736,240058,218645,188297, 
#' 		153931, 124347,93254,71858,53594,39721,27887,18092,34729 ) 
#' strong_result <- strong_smth(pop5m_pasex,Ages,TRUE)
#' # differences due to intermediate rounding in spreadsheet (bad practice IMO)
#' all(abs(strong_result - strongtest) < 1, na.rm = TRUE)
#' \dontrun{
#' plot(Ages, pop5m_pasex)
#' lines(as.integer(names(strong_result)),strong_result)
#' }
#' @references
#' \insertRef{carrier1959reduction}{DemoTools}
# plz add Demgen 1994 and PASEX. I can't download above citation while 
# on vacation. plz check if this method is discussed there as well.

strong_smth <- function(Value, 
		Age, 
		OAG = TRUE,
		ageMin = 10,
		ageMax = 65){
	
	# these values are not used, it's just for lengths, and to make sure we 
	# end on an even 10. Technically we could even provide data in 10-year
	# age groups and it'd still not break.
	Value1     <- splitUniform(Value = Value, Age = Age, OAG = OAG)
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
	indsub     <- Age10 >= ageMin & Age10 <= ageMax
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


#' G. Feeney's method of smoothing counts in 5-year age groups.
#' @description If age heaping is much worse on 0's than on 5's then even counts in 5-year age bins can preserve a sawtooth pattern. Most graduation techniques translate the zig-zag/sawtooth pattern to a wave pattern. It is not typically desired. This method redistributes counts 'from' every second 5-year age group in a specified range 'to' the adjacent age groups. How much to redistribute depends on a detection of roughness in the 5-year binned data, which follows the formulas recommended by Feeney. 
#' @details This function calls \code{zigzag()}, but prepares data in a way consistent with other methods called by \code{agesmth()}. It is probably preferable to call \code{zigzag()} from the top level, or else call this method from \code{agesmth()} for more control over tail imputations.
#' @param Value numeric vector of (presumably) counts in 5-year age groups. 
#' @param Age integer vector of age group lower bounds.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @param ageMin integer. Lower age bound to adjust values.
#' @param ageMax integer. Upper age bound to adjust values.
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export 
#' @references
#' Feeney, G. 2013 "Removing "Zigzag" from Age Data," http://demographer.com/white-papers/2013-removing-zigzag-from-age-data/
#' @examples
#' Value <- c(13331,4151,1746,1585,3859,8354,11146,12076,
#' 		12216,12016,12473,11513,12899,11413,12710,11516,11408,6733,4031,2069)
#' Age <- c(0,1,seq(5,90,by=5))
#' # defaults
#' zz <- zigzag_smth(Value, Age, OAG = TRUE, ageMin = 40, ageMax = 80)
#' \dontrun{
#' plot(Age, Value)
#' lines(as.integer(names(zz)),zz)
#' }
zigzag_smth <- function(Value, 
		Age, 
		OAG = TRUE, 
		ageMin = 40, 
		ageMax = max(Age) - max(Age)%%10 - 5){
	
	# insist on 5-year age groups
	Value <- groupAges(Value, Age = Age, N = 5)
	Age   <- as.integer(names(Value))
	
	Smoothed <- zigzag(
			Value = Value,
			Age = Age,
			OAG = OAG,
			ageMin = ageMin,
			ageMax = ageMax)
	# put NAs in for unadjusted ages,
	Smoothed[Smoothed == Value] <- NA
	
	Smoothed
}

#' Smooth in 5-year age groups using a moving average
#' @description Smooth data in 5-year age groups. 
#' @details This function calls \code{zigzag()}, but prepares data in a way consistent with other methods called by \code{agesmth()}. It is probably preferable to call \code{zigzag()} from the top level, or else call this method from \code{agesmth()} for more control over tail imputations.
#' @param Value numeric vector of (presumably) counts in 5-year age groups. 
#' @param Age integer vector of age group lower bounds.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @param n integer. The width of the moving average. Default 3 intervals (x-5 to x+9).
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @details This function calls \code{mav()}, which itself relies on the more general \code{ma()}. We lose the lowest and highest ages with this method, unless \code{n=1}, in which case data is returned in the original 5-year age groups.
#' @examples
#' Value <- c(13331,4151,1746,1585,3859,8354,11146,12076,
#' 		12216,12016,12473,11513,12899,11413,12710,11516,11408,6733,4031,2069)
#' Age <- c(0,1,seq(5,90,by=5))
#' # defaults
#' ns   <- sapply(1:5,mav_smth,Value=Value,Age=Age,OAG=TRUE)
#' cols <- paste0(gray(seq(.8,0,length=5)),"A0")
#' lwds <- seq(4,1,length=5)
#' \dontrun{
#' plot(Age, Value,pch=16)
#' matplot(as.integer(rownames(ns)),ns,type='l',
#' 		col = cols,
#' 		lty = 1,
#' 		add = TRUE,
#' 		lwd = lwds)
#' legend("topright", col = cols, lty = 1, lwd = lwds, legend = paste("n =",1:5))
#' }
#' @export

mav_smth <- function(Value, 
		Age, 
		OAG = TRUE,
		n = 3){
	Value <- groupAges(Value, Age = Age, N = 5)
	Age   <- as.integer(names(Value))
	
	Smoothed <- mav(
			Value = Value,
			Age = Age,
			OAG = OAG,
			n = n)
	
	Smoothed
}

#' Smooth populations in 5-year age groups using various methods
#' 
#' @description Smooth population counts in 5-year age groups using the Carrier-Farrag, 
#' Karup-King-Newton, Arriaga, United Nations, Stong, or Zigzag methods. Allows for imputation 
#' of values in the youngest and oldest age groups for the Carrier-Farrag, Karup-King-Newton, 
#' and United Nations methods.

#' @details The Carrier-Farrag, Karup-King-Newton, and Arriaga methods do not modify the totals 
#' in each 10-year age group, whereas the United Nations, Strong, Zigzag, and moving average (MAV) methods do. The age intervals 
#' of input data could be any integer structure (single, abridged, 5-year, other), but 
#' output is always in 5-year age groups. All methods except the United Nations and MAV methods
#' operate based on 10-year age group totals, excluding the open age group. 
#' 
#' The Carrier-Farrag, Karup-King-Newton, and United Nations methods do not produce estimates 
#' for the first and final 10-year age groups. By default, these are imputed with the original 5-year age group totals, but
#' you can also specify to impute with \code{NA}, or the results of the Arriaga or
#' Strong methods. If the terminal digit of the open age group is 5, then the terminal 10-year 
#' age group shifts down, so imputations may affect more ages in this case. Imputation can follow 
#' different methods for young and old ages. 
#' 
#' Method names are simplified using \code{simplify.text} and checked against a set of plausible matches 
#' designed to give some flexibility in case you're not sure 
#' 
#' In accordance with the description of these methods in Arriaga (1994), it is advised to 
#' compare the results from a variety of methods. 
#' 
#' @param Value numeric vector of counts in single, abridged, or 5-year age groups.
#' @param Age integer vector of ages corresponding to the lower integer bound of the counts.
#' @param method character string. Options include \code{"Carrier-Farrag"},\code{"Arriaga"},\code{"Karup-King-Newton"},\code{"United Nations"}, \code{"Strong"}, and \code{"Zigzag"}. See details. Default \code{"Carrier-Farrag"}.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @param ageMin integer. The lowest age included included in intermediate adjustment. Default 10. Only relevant for Strong method.
#' @param ageMax integer. The highest age class included in intermediate adjustment. Default 65. Only relevant for Strong method.
#' @param n integer. The width of the moving average. Default 3 intervals (x-5 to x+9). Only relevant for moving average method.
#' @param young.tail \code{NA} or character. Method to use for ages 0-9. See details. Default \code{"original"}.
#' @param old.tail \code{NA} or character. Method to use for the final age groups. See details. Default \code{"original"}.
#' 
#' @return numeric vector of smoothed counts in 5-year age groups.
#' @export
#' 
#' @examples
#' Ages         <- seq(0, 80, by = 5)
#'
#' # names a bit flexible:
#' cf <- agesmth(Value = pop5m_pasex, 
#'		Age = Ages, 
#'		method = "Carrier-Farrag", 
#'		OAG = TRUE)
#'# old.tail be default same as young.tail
#'# "cf" also works
#'
#'# no need to specify tails for Arriaga or Strong
#'arr <- agesmth(Value = pop5m_pasex, 
#'		Age = Ages, 
#'		method = "Arriaga", 
#'		OAG = TRUE)
#'strong <- agesmth(Value = pop5m_pasex, 
#'		Age = Ages, 
#'		method = "Strong", 
#'		OAG = TRUE)
#'# other methods:
#'un <- agesmth(Value = pop5m_pasex, 
#'		Age = Ages, 
#'		method = "United Nations", 
#'		OAG = TRUE)
#'kkn <- agesmth(Value = pop5m_pasex, 
#'		Age = Ages, 
#'		method = "Karup-King-Newton", 
#'		OAG = TRUE)
#' # zigzag, not plotted.
#' zz <- agesmth(pop5m_pasex,Ages,OAG=TRUE,method="Zigzag",ageMin = 30, ageMax = 70)
#' # mav, not plotted.
#' ma3 <- agesmth(pop5m_pasex,Ages,OAG=TRUE,method="MAV",n=3)
#' 
#'\dontrun{
#'	plot(Ages,pop5m_pasex,pch=16)
#'	lines(Ages, cf)
#'	lines(Ages, arr, col = "red")
#'	lines(Ages, strong, col = "#FF000080", lwd = 3)
#'	lines(Ages, kkn, col = "blue")
#'	lines(Ages, un, col = "magenta")
#'	legend("topright",
#'			pch=c(16,NA,NA,NA,NA,NA),
#'			lty = c(NA,1,1,1,1,1),
#'			lwd = c(NA,1,1,3,1,1),
#'			col = c("black","black","red","#FF000080","blue","magenta"),
#'			legend = c("orig 5","Carrier-Farrag","Arriaga","Strong","KKN","UN"))
#'}
#' # an extreme case:
#'  Age <- 0:99
#'  
#'  V5 <- groupAges(pop1m_pasex, Age=Age)
#'  Age5 <- as.integer(names(V5))
#'  cf2 <- agesmth(Value = pop1m_pasex, 
#'		  Age = Age, 
#'		  method = "Carrier-Farrag", 
#'		  OAG = TRUE)
#'  st2 <- agesmth(Value = pop1m_pasex, 
#'		  Age = Age, 
#'		  method = "Strong", 
#'		  OAG = TRUE)
#'  \dontrun{
#'  plot(Age,pop1m_pasex,pch=16)
#'  lines(Age,splitUniform(V5,Age=Age5,OAG=FALSE), lty=2, lwd = 2)
#'  lines(Age,splitUniform(cf2,Age=Age5,OAG=FALSE),col="blue")
#'  lines(Age,splitUniform(st2,Age=Age5,OAG=FALSE),col="red")
#'  legend("topright",
#'		  pch=c(16,NA,NA,NA),
#'		  lty=c(NA,2,1,1),
#'		  col=c("black","black","blue","red"),
#'		  lwd=c(NA,2,1,1),
#'		  legend=c("orig single","orig 5","Carrier-Farrag","Strong"))
#'}
#'  
#'# it might make sense to do this level of smoothing as intermediate step
#'# in Sprague-like situation. Compare:
#'spr1 <- sprague(Value, Age=Age,OAG=FALSE)
#'spr2 <- sprague(cf2, Age=Age5,OAG=FALSE)
#'spr3 <- sprague(st2, Age=Age5,OAG=FALSE)
#'\dontrun{
#'plot(Age,Value,pch=16, main = "Smoothing as pre-step to graduation")
#'lines(Age,spr1,lty=2)
#'lines(Age,spr2,col="blue")
#'lines(Age,spr3,col="red")
#'legend("topright",
#'		pch=c(16,NA,NA,NA),
#'		lty=c(NA,2,1,1),
#'		col=c("black","black","blue","red"),
#'		lwd=c(NA,2,1,1),
#'		legend=c("orig single","orig->Sprague","C-F->Sprague","Strong->Sprague"))
#'}
#' 
agesmth <- function(Value, 
		Age, 
		method = c("Carrier-Farrag","KKN","Arriaga","United Nations","Strong","Zigzag","MAV")[1], 
		OAG = TRUE, 
		ageMin = 10, ageMax = 65,
		n = 3,
		young.tail = c("Original","Arriaga","Strong",NA)[1], 
		old.tail = young.tail){
	method     <- simplify.text(method)
	young.tail <- simplify.text(young.tail)
	old.tail   <- simplify.text(old.tail)
	
	if (missing(Age)){
		Age <- as.integer(names(Value))
	}
	stopifnot(length(Age) == length(Value))
	# carrierfarrag or cf
	if (method %in% c("cf", "carrierfarrag")){
		out <- carrier_farrag_smth(Value = Value, Age = Age, OAG = OAG)
	}
	
	# stong
	if (method == "strong"){
		out <- strong_smth(Value = Value, Age = Age, OAG = OAG, ageMin = ageMin, ageMax = ageMax)
	}
	
	# un or unitednations
	if (method %in% c("un", "unitednations")){
		out <- united_nations_smth(Value = Value, Age = Age, OAG = OAG)
	}
	
	# arriaga
	if (method  == "arriaga"){
		out <- arriaga_smth(Value = Value, Age = Age, OAG = OAG)
	}
	
	# kkn kkingnewton karupkingnewton
	if (method %in% c("kkn", "kkingnewton", "karupkingnewton")){
		out <- kkn_smth(Value = Value, Age = Age, OAG = OAG)
	}
	# TR: new Feeney method added July 31, 2018
	if (method %in% c("feeney","zigzag")){
		# however, need to make it so NAs returned in unaffected ages?
		# or make the user call it in various runs and graft together.
		out <- zigzag_smth(Value = Value, Age = Age, OAG = OAG, ageMin = ageMin, ageMax = ageMax)
	}
	# TR: MAV added Aug 7
	if (method %in% c("mav","ma","movingaverage")){
		# however, need to make it so NAs returned in unaffected ages?
		# or make the user call it in various runs and graft together.
		out <- mav_smth(Value = Value, Age = Age, OAG = OAG, n = n)
	}
	# -------------------------------
	# clean tails
	nas     <- is.na(out)
	if (any(nas) & (!is.na(old.tail) | !is.na(young.tail))){
		nrle         <- rle(as.integer(nas))
		original     <- groupAges(Value, Age = Age, N = 5)
		arriaga      <- arriaga_smth(Value, Age = Age, OAG = OAG)
		strong       <- strong_smth(Value, Age = Age, OAG = OAG)
		# are the final entries NAs?
		if (nrle$values[length(nrle$values)] == 1 & !is.na(old.tail)){
			nrle$values[1] <- 0
			old.ind        <- as.logical(rep(nrle$values, times = nrle$lengths))
			# do we want original values?
			if (old.tail %in% c("o","orig","original")){
				stopifnot(length(original) == length(out))
				out[old.ind] <- original[old.ind]
			}
			# or arriaga?
			if (old.tail == "arriaga"){
				stopifnot(length(arriaga) == length(out))
				out[old.ind] <- arriaga[old.ind]
			}
			# or strong?
			if (old.tail == "strong"){
				stopifnot(length(strong) == length(out))
				out[old.ind] <- strong[old.ind]
			}
			
		}
		nrle             <- rle(as.integer(nas))
		# take care of young tail
		if (nrle$values[1] == 1 & !is.na(young.tail)){
			nrle$values[length(nrle$values)] <- 0
			young.ind        <- as.logical(rep(nrle$values, times = nrle$lengths))
			
			if (young.tail %in% c("o","orig","original")){
				stopifnot(length(original) == length(out))
				out[young.ind] <- original[young.ind]
			}
			# or arriaga?
			if (young.tail == "arriaga"){
				stopifnot(length(arriaga) == length(out))
				out[young.ind] <- arriaga[young.ind]
			}
			# or strong?
			if (young.tail == "strong"){
				stopifnot(length(strong) == length(out))
				out[young.ind] <- strong[young.ind]
			}
		}
	} # end tail operations
	
	out
}




