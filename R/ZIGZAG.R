# Author: Tim Riffe, based on earlier version from Juan Galeano
###############################################################################

#' Objective function to minimize Feeney's zigzag method residual
#' @description This function is auxiliary to \code{zigzag()}, see \code{?zigzag} for a description.
#' @details This function is not intended to be used at the top level, but just in case, make sure that \code{ageMax = ageMin + 10 * length(p)}. Age groups \code{ >= ageMin} AND \code{ <= ageMin} must be in 5-year age groups. This function does no checks.
#' @param Value numeric vector of (presumably) counts in 5-year age groups. 
#' @param Age integer vector of age group lower bounds.
#' @param ageMin integer. Lower age bound to adjust values.
#' @param ageMax integer. Upper age bound to adjust values.
#' @param p numeric vector of adjustment parameters.
#' @return positive residual to minimize.
#' @export
#' @importFrom stats optim
#' @references
#' Feeney, G. 2013 "Removing "Zigzag" from Age Data," http://demographer.com/white-papers/2013-removing-zigzag-from-age-data/
#' @examples
#'Age <- c(0,1,seq(5,90,by=5))
#'zigzag_min(dth5_zigzag,Age,ageMin=40,ageMax = 80,p=rep(.05,4))
#'# it's used like this in zigzag()
#'(p <- optim(
#'			rep(.05,4), 
#'			zigzag_min, 
#'			Value = dth5_zigzag, 
#'			Age = Age, 
#'			ageMin = 40, 
#'			ageMax = 80)$par)
#' Smoothed <- zigzag_p(dth5_zigzag,Age,40,80,p)
#' # de facto unit test:
#' # check result using results frozen in Feeney spreadsheet
#' # after fixing probable cell range 'error'
#' p.feeney <- c(0.0235802695087692,0.0724286618207911,
#' 		      0.0242327829742267,0.0883411499065237)
#' ans      <- 106.1147411629
#' stopifnot(abs(zigzag_min(dth5_zigzag, Age, 40,80,p.feeney) - ans) < 1e-6)
zigzag_min <- function(Value, Age, ageMin = 40, ageMax = 80, p){
	
	# first, we need an odd number of age groups to smooth over.
	# we enforce this
	
	
	Smoothed <- zigzag_p(
			Value = Value, 
			Age = Age, 
			ageMin = ageMin, 
			ageMax = ageMax, 
			p = p)
	
	# multiplier
	id <- abs(Smoothed - Value) > 1e-6
	
	# adjacent avg
	avg2  <- avg_adj(Smoothed) * id
	
	# difference squared
	diffsq <- (Smoothed - avg2)^2 * id
	
	# return weighted sum sq dev
	sum(diffsq * 1 / Smoothed, na.rm=TRUE)
	
}

#' Smooth population counts using Feeney's zigzag method and smoothing parameters.
#' @description This function is auxiliary to \code{zigzag()}, see \code{?zigzag} for a description.
#' @details This function is not intended to be used at the top level, but just in case, make sure that \code{ageMax = ageMin + 10 * length(p)}. Age groups \code{ >= ageMin} AND \code{ <= ageMin} must be in 5-year age groups. This function does no checks.
#' @param Value numeric vector of (presumably) counts in 5-year age groups. 
#' @param Age integer vector of age group lower bounds.
#' @param ageMin integer. Lower age bound to adjust values.
#' @param ageMax integer. Upper age bound to adjust values.
#' @param p numeric vector of adjustment parameters.
#' @return numeric vector of smoothed population counts.
#' @export
#' @references
#' Feeney, G. 2013 "Removing "Zigzag" from Age Data," http://demographer.com/white-papers/2013-removing-zigzag-from-age-data/
#' @examples
#'Age <- c(0,1,seq(5,90,by=5))
#'zigzag_min(dth5_zigzag,Age,ageMin=40,ageMax = 80,p=rep(.05,4))
#'# it's used like this in zigzag()
#'(p <- optim(
#'			rep(.05,4), 
#'			zigzag_min, 
#'			Value = dth5_zigzag, 
#'			Age = Age, 
#'			ageMin = 40, 
#'			ageMax = 80)$par)
#' Smoothed <- zigzag_p(dth5_zigzag,Age,40,80,p)

zigzag_p <- function(
		Value, 
		Age, 
		ageMin = 40, 
		ageMax = 80, 
		p = rep(0.05, floor((ageMax - ageMin) / 10))){
	
	# 1) p used differently
	p_ages      <- as.integer(
			seq(ageMin + 5, ageMin + length(p) * 10 - 5, length = length(p))
	)
	mid5        <- Age %in% p_ages
	down5       <- Age %in% (p_ages - 5)
	up5         <- Age %in% (p_ages + 5)
	
	# 2) this can be a multiplier later
	id          <- rep(0,length(Age))
	id[Age >= (min(p_ages) - 5) & Age <= (max(p_ages) + 5)] <- 1
	
	# 2) average of value above and below
	avg1        <- avg_adj(Value) * id
	
	# 3) take difference
	diff1       <- (Value - avg1) * id
	
	# 5) get Np
	Np          <- Value[mid5] * p # should be 4 numbers in example
	
	# 6) use Np to smooth
	Smoothed                <- Value
	Smoothed[mid5]          <- Value[mid5] - Np
	Smoothed[down5]         <- Smoothed[down5] + Np / 2
	Smoothed[up5]           <- Smoothed[up5] + Np / 2
	
	# return smoothed
	names(Smoothed) <- Age
	Smoothed
}

#' G. Feeney's method of removing the zigzag from counts in 5-year age groups.
#' @description If age heaping is much worse on 0's than on 5's then even counts in 5-year age bins can preserve a sawtooth pattern. Most graduation techniques translate the zig-zag/sawtooth pattern to a wave pattern. It is not typically desired. This method redistributes counts 'from' every second 5-year age group in a specified range 'to' the adjacent age groups. How much to redistribute depends on a detection of roughness in the 5-year binned data, which follows the formulas recommended by Feeney. 
#' @details Determining the degree to redistribute population counts is an optimization problem, so this function has two auxiliary functions, \code{p_zigzag()}, which redistributes counts according to a set of specified proportions \code{p}, and \code{zigzag_min()} which is the function minimized to determine the optimal vector of \code{p}. If data is not in 5-year age groups then it is grouped as necessary (unless abridged, in which case grouping is preserved). Only ages \code{>= ageMin} and \code{<= ageMax} will be adjusted. \code{ageMax} is inclusive and interpreted as the lower bound of the highest age group to include. The number of 5-year age groups adjusted must be odd, and \code{ageMax} may be reduced internally without warning in order to force this condition. The open age group is excluded from adjustment. This function also has a wrapper \code{zigzag_smth()}, called by \code{agesmth()}.
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
#' Age <- c(0,1,seq(5,90,by=5))
#' # defaults
#' zz1 <- zigzag(dth5_zigzag, Age, OAG = TRUE, ageMin = 40, ageMax = 80)
#' # shift age range up by 5
#' zz2 <- zigzag(dth5_zigzag, Age, OAG = TRUE, ageMin = 45, ageMax = 85)
#' \dontrun{
#' plot(Age, dth5_zigzag,pch = 16)
#' lines(Age,zz1,col="red",lty=1)
#' lines(Age,zz2,col="blue",lty=1)
#' # even smoother:
#' lines(Age,(zz2+zz1)/2,col="green",lty=1)
#' legend("bottom",pch=c(16,NA,NA,NA),
#' 		lty=c(NA,1,1,1),
#' 		col = c("black","red","blue","green"),
#' 		legend = c("Original","40-80","45-85","avg"))
#' }

zigzag <- function(
		Value, 
		Age, 
		OAG = TRUE, 
		ageMin = 40, 
		ageMax = max(Age) - max(Age)%%10 - 5){
	stopifnot(length(Age) == length(Value))
	
	Age        <- as.integer(Age)
	
	if (!is_abridged(Age)){
		Value5     <- groupAges(Value, Age = Age, N = 5)
		Age5       <- as.integer(names(Value5))
	} else {
		Value5     <- Value
		Age5       <- Age
	}
	
	N          <- length(Value5)
	NN         <- N
	# make sure we have a suitable span of age.
	Ndecades <- floor((ageMax - ageMin)/10) * 10
	if ((ageMax - ageMin) - Ndecades != 5){
		ageMax   <- ageMin + Ndecades + 5
	}
	# enforce 5-year age groups. Data could be given as such,
	# in abridged ages, or in single ages, but will be returned
	# in 5-year age groups.
	if (OAG){
		OAGvalue  <- Value5[N]
		Value5[N] <- NA
		NN        <- N - 1
	}
	
	# make sure ageMax isn't too high. Likely doesn't come into play.
	if (ageMax > max(Age5[1:NN])){
		hm     <- ceiling((ageMax - max(Age5[1:NN]))/10) * 10
		ageMax <- ageMax - hm
	}
	Ndecades <- floor((ageMax - ageMin)/10)
	# get starting values for p. Number of values to optimize
	# depends on age range in question.
	p        <- rep(0.05, Ndecades)
	# get optimal p
	p        <- optim(p, zigzag_min, Value = Value5, Age = Age5, ageMin = ageMin, ageMax = ageMax)$par
	# and resulting smoothed pops
	Smoothed <- zigzag_p(Value = Value5, Age = Age5, ageMin = ageMin, ageMax = ageMax, p = p)
	
	# put OAG back.
	if (OAG){
		Smoothed[N] <- OAGvalue
	}
	names(Smoothed) <- Age5
	
	Smoothed
}





