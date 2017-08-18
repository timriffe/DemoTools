
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


