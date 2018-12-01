
# Author: tim
# a bunch of identities that assume abridged ages, but will work for single ages
# if AgeInt is a bunch of 1s
###############################################################################


#' Derive nMx from nqx and nax.
#' @description This is the standard identity to derive nMx from nax and nqx.
#' 
#' @param nqx numeric. Vector of age specific death probabilities.
#' @param nax numeric. Vector of average time spent in interval by those dying in interval.
#' @param AgeInt integer. Vector of age class widths.
#' @references 
#' \insertRef{preston2000demography}{DemoTools}
#' @return nMx vector of age specific death rates derived via identity.
#' @export
qxax2mx <- function(nqx, nax, AgeInt = inferAgeIntAbr(vec = nqx)) {
	nqx / (AgeInt - (AgeInt - nax) * nqx)      
}

#' Derive nqx from nMx and nax.
#' @description This is the standard identity to derive nqx from nax and nMx.
#' @details qx values calculated as greater than 1 are imputed with 1.
#' @param nMx numeric. Vector of age-specific death rates.
#' @param nax numeric. Vector of average time spent in interval by those dying in interval.
#' @param AgeInt integer. Vector of age class widths.
#' @references 
#' \insertRef{preston2000demography}{DemoTools}
#' @return nqx vector of age specific death probabilities derived via identity.
#' @export
mx2qx <- function(nMx, nax, AgeInt = inferAgeIntAbr(vec = nMx)) {
	qx <- (AgeInt * nMx) / (1 + (AgeInt - nax) * nMx)
	if (any(qx > 1)){
		cat("at least 1 q(x) > 1, imputed as 1")
		qx[qx > 1] <- 1
	}
	qx
}

#' Derive nax from nqx and nMx.
#' @description This is the standard identity to derive nax from nqx and nMx.
#' 
#' @param nqx numeric. Vector of age specific death probabilities.
#' @param nMx numeric. Vector of age-specific death rates.
#' @param AgeInt integer. Vector of age class widths.
#' @references 
#' \insertRef{preston2000demography}{DemoTools}
#' @return nax numeric vector of average time spent in interval by those dying in interval via identity.
#' @export
qxmx2ax <- function(nqx, nMx, AgeInt){
	1 / nMx - AgeInt / nqx + AgeInt
}

#' Derive nqx from nMx and nax.
#' @description This is the standard identity to derive nqx from nax and nMx.
#' This is a more full-service wrapper of \code{mx2qx()}, with closeout options and optional age 0 
#' treatment.
#' @details qx values calculated as greater than 1 are imputed with 1.
#' 
#' @param nMx numeric. Vector of age-specific death rates.
#' @param nax numeric. Vector of average time spent in interval by those dying in interval.
#' @param AgeInt integer. Vector of age class widths.
#' @param closeout logical. Set to 1 if TRUE. Identity otherwise. Default \code{TRUE}.
#' @param IMR numeric. Optional q0 to impute, in case available separately.
#' @references 
#' \insertRef{preston2000demography}{DemoTools}
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
	if (any(qx > 1)){
		cat("at least 1 q(x) > 1, imputed as 1")
		qx[qx > 1] <- 1
	}
	qx
}

#' Derive lifetable survivorship (lx) from death probabilities.
#' @description This lifetable identity is the same no matter what kind of lifetable is required.
#'  You can find it in any demography textbook.
#' @details set \code{radix = 1} for the probability of surviving until age x. The vector returned is
#' the same length as \code{nqx}, thereby throwing out the final value of qx, which is usually set to 1.
#' 
#' @param nqx numeric. Vector of age specific death probabilities.
#' @param radix numeric. The lifetable starting population. Default 100000.
#' @references 
#' \insertRef{preston2000demography}{DemoTools}
#' @return lx vector of lifetable survivorship.
#' @export
qx2lx <- function(nqx, radix = 1e5){
	radix * cumprod(c(1, 1 - nqx[-length(nqx)]))
}

#' Derive lifetable deaths from survivorship.
#' @description This lifetable identity is the same no matter what kind of lifetable is required.
#'  You can find it in any demography textbook.
#' @details The vector returned is the same length as \code{lx} and it sums to the lifetable radix. 
#' If the radix is one then this is the discrete deaths distribution.
#' 
#' @param lx numeric.  Vector of age-specific lifetable survivorship.
#' @references 
#' \insertRef{preston2000demography}{DemoTools}
#' @return ndx vector of lifetable deaths.
#' @export
lx2dx <- function(lx){
	diff(-c(lx,0))
}

#' Derive lifetable exposure from lx, ndx and nax.
#' @description This is a common approximation of lifetable exposure: 
#' All persons surviving to the end of the interval time the interval width, plus all those that died 
#' in the interval multiplied by their average time spent in the interval.
#' @details There is no checking of equal vector lengths at this time.
#' @param lx numeric.  Vector of age-specific lifetable survivorship.
#' @param ndx numeric. Vector of lifetable deaths, summing to radix of \code{lx}.
#' @param nax numeric. Vector of average time spent in interval by those dying in interval.
#' @param AgeInt integer. Vector of age class widths.
#' @references 
#' \insertRef{preston2000demography}{DemoTools}
#' @return nLx numeric vector of lifetable exposure.
#' @export
lxdxax2Lx <- function(lx, ndx, nax, AgeInt){
	N                   <- length(lx)
	nLx                 <- rep(0, N)
	nLx[1:(N - 1)]      <- AgeInt[1:(N - 1)] * lx[2:N] + nax[1:(N - 1)] * ndx[1:(N - 1)]
	nLx[N]		        <- lx[N] * nax[N]
	nLx
}

#' Derive lifetable total person years left to live from exposure.
#' @description A lifetable identity. Tx is interpreted as the total years
#' left to live above age x in the life table stationary population.
#' @details No \code{NA} or other error checking here. This is taken as the numerator in the classic lifetable method
#'  of calculation of life expectancy.
#' @param Lx numeric. Vector of lifetable exposure.
#' @references 
#' \insertRef{preston2000demography}{DemoTools}
#' @return Tx total years left to live above age x.
#' @export
Lx2Tx <- function(Lx){
	rev(cumsum(rev(Lx)))
}


