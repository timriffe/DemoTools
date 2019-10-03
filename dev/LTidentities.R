
# Author: tim
# a bunch of identities that assume abridged ages, but will work for single ages
# if AgeInt is a bunch of 1s
###############################################################################


#' derive nMx from nqx and nax
#' @description This is the common ax-qx-mx identity as found in any demography text.
#' 
#' @param nqx numeric vector. Age specific death probabilities.
#' @param ax numeric vector. Average time spent in interval by those dying in interval
#' @param AgeInt integer width of each age interval
#' 
#' @return nMx vector of age specific death rates derived via identity.
#' @export
qxax2mx <- function(nqx, nax, AgeInt = inferAgeIntAbr(vec = nqx)) {
	nqx / (AgeInt - (AgeInt - nax) * nqx)      
}

#' derive nqx from nMx and nax
#' @description This is the common ax-qx-mx identity as found in any demography text.
#' 
#' @param nMx numeric vector. Age specific death rates.
#' @param nax numeric vector. Average time spent in interval by those dying in interval
#' @param AgeInt integer width of each age interval
#' 
#' @return nqx vector of age specific death probabilities derived via identity.
#' @export
mx2qx <- function(nMx, nax, AgeInt = inferAgeIntAbr(vec = nMx)) {
	(AgeInt * nMx) / (1 + (AgeInt - nax) * nMx)
}

#' derive nax from nqx and nMx
#' @description This is the common ax-qx-mx identity as found in any demography text.
#' 
#' @param nMx numeric vector. Age specific death rates.
#' @param nax numeric vector. Average time spent in interval by those dying in interval
#' @param AgeInt integer width of each age interval
#' 
#' @return nax numeric vector of average time spent in interval by those dying in interval
#' @export
qxmx2ax <- function(nqx, nMx, AgeInt){
	1 / nMx - AgeInt / nqx + AgeInt
}

#' derive nqx from nMx and nax
#' @description This is the common ax-qx-mx identity as found in any demography text. 
#' This is a more full-service wrapper of \code{mx2qx()}, with closeoput options and optional age 0 
#' treatment.
#' 
#' @param nMx numeric vector. Age specific death rates.
#' @param nax numeric vector. Average time spent in interval by those dying in interval
#' @param AgeInt integer width of each age interval
#' @param closeout logical (default \code{TRUE}). Set to 1 if TRUE. Identity otherwise
#' @param IMR numeric. Optional q0 to impute, in case available separately.
#' 
#' @return nqx vector of age specific death probabilities derived via identity.
#' @export
mxax2qx <- function(nMx, nax, AgeInt, closeout = TRUE, IMR){
	
	qx <- mx2qx(nMx, nax, AgeInt)
	if (closeout){
		qx[length(qx)] <- 1
		if (length(nMx) == 1){
			warning("only a single nMx was given, and it was treated as age omega, with qx = 1.\nSpecify closeout = FALSE otherwise")
		}
	}
	if (!missing(IMR) & !is.na(IMR)){
		qx[1] <- IMR
	}
	qx
}

#' derive lifetable survivorship from death probabilities
#' @description This lifetable identity is the same no matter what kind of lifetable
#' we're making. You can find it in any demography textbook.
#' @details set \code{radix = 1} for the probability of surviving until age x. The vector returned is
#' the same length as \code{nqx}, thereby throwing out the final value of qx, which is probably set to 1 anyway.
#' 
#' @param nqx numeric vector. Age-specific death probabilities in any age classes.
#' @param radix numeric (default 100000). The lifetable starting population.
#' @return lx numeric vector of lifetable survivorship
#' @export
qx2lx <- function(nqx, radix = 1e5){
	radix * cumprod(c(1, 1 - nqx[-length(nqx)]))
}
#' derive lifetable deaths from survivorship
#' @description This lifetable identity is the same no matter what kind of lifetable
#' we're making. You can find it in any demography textbook.
#' @details The vector returned is the same length as \code{lx} and it sums to the lifetable radix. 
#' If the radix is one then this is the discrete deaths distribution.
#' 
#' @param lx numeric vector. Age-specific lifetable survivorship
#' @return ndx numeric vector of lifetable deaths.
#' @export
lx2dx <- function(lx){
	diff(-c(lx,0))
}

#' derive lifetable exposure from lx, ndx and nax.
#' @description This is a common approximation of lifetable exposure: 
#' All persons surviving to the end of the interval time the interval width, plus all those that died 
#' in the interval multiplied by their average time spent in the interval.
#' @details There is no checking of equal vector lengths at this time.
#' @param lx numeric vector of lifetable survivorship
#' @param ndx numeric vector of lifetable deaths, summing to radix of \code{lx}.
#' @param nax numeric vector of average time spent in the age interval of those dying in the interval.
#' @param AgeInt integer vector of age class widths
#' @return nLx numeric vector of lifetable exposure.
#' @export
lxdxax2Lx <- function(lx, ndx, nax, AgeInt){
	N                   <- length(lx)
	nLx                 <- rep(0, N)
	nLx[1:(N - 1)]      <- AgeInt[1:(N - 1)] * lx[2:N] + nax[1:(N - 1)] * ndx[1:(N - 1)]
	nLx[N]		        <- lx[N] * nax[N]
	nLx
}

#' derive lifetable total person years left to live from exposure
#' @description A lifetable identity. Tx is interpreted as the total years
#' left to live above age x in the life table stationary population.
#' @details No \code{NA} or other error checking here. This is taken as the numerator in the classic lifetable method
#'  of calculation life expectancy.
#' @param Lx numeric vector of lifetable exposure.
#' @return Tx total years left to live above age x
#' @export
Lx2Tx <- function(Lx){
	rev(cumsum(rev(Lx)))
}


#' calculate survivor ratios
#' @description An extra lifetable column for use in projections, which require uniform time steps both both age and period. Intervals are either single age (\code{N=1}) or five-year ages (\code{N=5}). Input vectors are assumed to come from either single or standard abridged ages.
#' @details This function does not account for \code{nLx} having been pre-binned into uniform 5-year age widths, which will throw an error. Just leave them in abridged ages instead. Note that in the case of abridged ages, the interpretation for the first and second value don't follow the original abridged age intervals: the first value in the probability of surviving from birth into ages 0-4 in the first five years, and the second value is the probability of surviving from 0-4 to 5-9. This represents a slight misalignment with the rest of the lifetable, user beware.
#' @inheritParams LTabr 
#' @param N integer, the age width for survivor ratios, either 5 or 1. Default 5.
#' @export
Lxlx2Sx <- function(nLx, lx, AgeInt, N = c(5,1)){
  n  <- length(nLx)
  stopifnot(length(lx) == n)
  # either we're in 1 or 5 year age groups
  N  <- match.arg(N)
  ## compute Sx (missing from the LTbr computation
  Sx <- rep(NA, n)
  # first age group is survival from births to the second age group		
  if (N == 5){
    # double check because assuming abridged nLx is given...
    stopifnot(length(AgeInt) == n)
    ageintcompare <- inferAgeIntAbr(vec=nLx)
    stopifnot(all(ageintcompare[-n] == AgeInt[-n]))
    # birth until 0-4
    Sx[1]      <- (nLx[1] + nLx[2]) / ((AgeInt[1] + AgeInt[2]) * lx[1])
    # second age group is survival age 0-4 to age 5-9
    Sx[2]      <- nLx[3] / (nLx[1] + nLx[2])
    # middle age groups 
    mind       <- 3:(n - 2)
    Sx[mind]   <- nLx[mind + 1] / nLx[mind]
  }
  if (N == 1){
    LLXX       <- c(lx[1], nLx)
    mind       <- 1:(n - 2)
    Sx[mind]   <-  LLXX[mind + 1] / LLXX[mind]
  }
  
  # penultimate age group
  Sx[n - 1]    <- nLx[n] / (nLx[n - 1] + nLx[n])	
  # closeout
  Sx[n]        <- 0.0
  
  Sx
}