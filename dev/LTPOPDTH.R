
# Author: tim
###############################################################################


# 
source("/home/tim/git/DemoTools/dev/nAx.R")

#' derive nMx from nqx and nax
#' @description This is the common ax-qx-mx identity as found in any demography text.
#' 
#' @param nqx numeric vector. Age specific death probabilities.
#' @param ax numeric vector. Average time spent in interval by those dying in interval
#' @param AgeInt integer width of each age interval
#' 
#' @return nMx vector of age specific death rates derived via identity.
#' @export
qx2mx <- function(nqx, nax, AgeInt = inferAgeIntAbr(vec = nqx)) {
	nqx / (AgeInt - (AgeInt - nax) * nqx)      }

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


qx2lx <- function(nqx, radix=1e5){
	
	radix * cumprod(c(1, 1 - nqx[-length(nqx)]))
}

lx2dx <- function(lx){
	diff(-c(lx,0))
}
# dx <- ndx; ax <- nAx
lxdxax2Lx <- function(lx,dx,ax,AgeInt){
	N                   <- length(lx)
	Lx                  <- rep(0,N)
	Lx[1:(N - 1)]       <- AgeInt[1:(N - 1)] * lx[2:N] + ax[1:(N - 1)] * dx[1:(N - 1)]
	Lx[N]		        <- lx[N] * ax[N]
	Lx
}
Lx2Tx <- function(Lx){
	rev(cumsum(rev(Lx)))
}




#' @examples
#' # trial code from PAS LTPOPDTH, North, Males, IMR = .1
#' Exposures <- c(100958,466275,624134,559559,446736,370653,301862,249409,
#' 		247473,223014,172260,149338,127242,105715,79614,53660,
#' 		31021,16805,8000,4000,2000,1000)
#' 
#' Deaths <- c(8674,1592,618,411,755,1098,1100,1357,
#' 		1335,3257,2200,4023,2167,4578,2956,4212,
#' 		2887,2351,1500,900,500,300)
#' # lower age bounds
#' Age    <- c(0, 1, seq(5, 100, by = 5))
#' AgeInt <- c(diff(Age), NA)
#' 
#' 
#' 
#' PASLT <- LTabr(Deaths = Deaths, 
#' 		Exposures = Exposures, 
#' 		Age = Age,
#' 		AgeInt =AgeInt,
#' 		axmethod = "PAS",
#' 		IMR = .1,
#' 		region = "N",
#' 		Sex = "M")
#' # de facto unit test. The unsmoothed output from PAS
#' # spreadsheet (allow for rounding error in last decimal)
#' excheck <- c(56.31,61.53,58.35,53.63, 
#' 		48.81,44.21,39.83,35.52,
#' 		31.43,27.22,24.09,20.52,
#' 		18.12,14.51,12.42,9.45,7.85, 
#' 		6.09,4.95,4.28,3.85,3.33)
#' stopifnot(max(abs(round(PASLT$ex,2) - excheck)) <= .01)
#' # examples based on UN 1982 (p. 34)
#' Mx <- c(.23669,.04672,.00982,.00511,.00697,.01036,.01169,
#' 		.01332,.01528,.01757,.02092,.02517,.03225,.04241,.06056,
#' 		.08574,.11840,.16226,.23745)
#' excheckUN <-  c(35.000,42.901,47.190,44.438,
#' 		40.523,36.868,33.691,30.567,27.500,24.485,21.504,18.599,
#' 		15.758,13.080,10.584,8.466,6.729,5.312,4.211)
#' AgeInt <- inferAgeIntAbr(vec = Mx)
#' 
#' 
#' UNLT1 <- LTabr(nMx = Mx,
#' 		Age = c(0,1,seq(5,85,by=5)),
#' 		AgeInt = AgeInt,
#' 		axmethod = "UN",
#' 		Sex = "M", 
#' 		kludge = FALSE)
#' UNLT2 <- LTabr(nMx = Mx,
#' 		Age = c(0,1,seq(5,85,by=5)),
#' 		AgeInt = AgeInt,
#' 		axmethod = "UN",
#' 		Sex = "M", 
#' 		kludge = TRUE)
#' # TODO: possibly there is more intermediate rounding happening in the FOR code?
#' round(UNLT2$ex,2) - round(excheckUN,2)
#' UNLT2$ex - UNLT1$ex
LTabr <- function(
		Deaths, 
		Exposures, 
		nMx,
		Age,
		AgeInt = inferAgeIntAbr(Age = Age), 
		axmethod = "pas", 
		Sex = "m", 
		region = "w",
		SRB = 1.05, 
		IMR,
		kludge = FALSE){
	
	if (missing(Deaths)){
		Deaths    <- NULL
	}
	if (missing(Exposures)){
		Exposures <- NULL
	}
	
	axmethod      <- tolower(axmethod)
	Sex           <- tolower(Sex)
	region        <- tolower(region)
	
	if (missing(nMx)){
		nMx       <- Deaths / Exposures
    }
	if (missing(IMR)){
		IMR       <- ifelse(axmethod == "pas", nMx[1],NA)
	}
	
	# take care of ax first, two ways presently
	if (axmethod == "pas"){
	    nAx       <- axPAS(
				      nMx, 
				      AgeInt = AgeInt, 
				      IMR = IMR, 
				      Sex = Sex, 
				      region = region)
      }
	if (axmethod == "un"){
		# UN method just CD west for now, so no region arg
		nAx       <- axUN(
				      nMx = nMx, 
				      AgeInt = AgeInt, 
				      Sex = Sex,
					  kludge = kludge)
	}
    # these are the same always. Only thing different is ax method.
	nqx <- mxax2qx(nMx = nMx, nax = nAx, AgeInt = AgeInt, closeout = TRUE, IMR = IMR)
	# TODO: differences might begin here
	lx  <- qx2lx(nqx)
	ndx <- lx2dx(lx)
	nLx <- lxdxax2Lx(lx = lx, dx = ndx, ax = nAx, AgeInt = AgeInt)
	Tx  <- Lx2Tx(nLx)
	ex  <- Tx / lx
	
	# output is an unrounded, unsmoothed lifetable
	out <- data.frame(
			   Age = Age,
			   AgeInt = AgeInt,
			   nMx = nMx,
			   nAx = nAx,
			   nqx = nqx,
			   lx = lx,
			   ndx = ndx,
			   nLx = nLx,
			   Tx = Tx,
			   ex = ex)
	return(out)
}






# example data from United Nations Model Life Tables -
# Males Latin American Pattern

Mx <- c(.23669,.04672,.00982,.00511,.00697,.01036,.01169,
.01332,.01528,.01757,.02092,.02517,.03225,.04241,.06056,
.08574,.11840,.16226,.23745)
excheckUN <-  c(35.000,42.901,47.190,44.438,
40.523,36.868,33.691,30.567,27.500,24.485,21.504,18.599,
15.758,13.080,10.584,8.466,6.729,5.312,4.211)
AgeInt <- inferAgeIntAbr(vec = Mx)


UNLT1 <- LTabr(nMx = Mx,
		Age = c(0,1,seq(5,85,by=5)),
		AgeInt = AgeInt,
		axmethod = "UN",
		Sex = "M", 
		kludge = FALSE)
UNLT2 <- LTabr(nMx = Mx,
		Age = c(0,1,seq(5,85,by=5)),
		AgeInt = AgeInt,
		axmethod = "UN",
		Sex = "M", 
		kludge = TRUE)
# TODO: possibly there is more intermediate rounding happening in the FOR code?
round(UNLT2$ex,2) - round(excheckUN,2)
UNLT2$ex - UNLT1$ex

Mx2 <-  c(.16326,.01499,.00303,.00235,
.00401,.00583,.00685,.00805,
.00922,.01059,.01290,.01654,
.02315,.03290,.04950,.07084,
.10214,.14112,.21389)
UNLT2 <- LTabr(nMx = Mx2,
		Age = c(0,1,seq(5,85,by=5)),
		AgeInt =inferAgeIntAbr(vec = Mx),
		axmethod = "UN",
		Sex = "F")
qxmx2ax(0.16326000, 0.16326, 1)
qx2mx(0.16326, .35, 1)
mx2qx(0.16326, .35, 1)

# a MORTPACK unit test:
#MPnMx <- c(0.12846,0.02477,0.00603,0.0034,
#0.00417,0.00513,0.00581,0.00645,0.00725,
#0.00813,0.00913,0.01199,0.01647,
#0.0256,0.04047,0.06624,0.10638,0.19611)
#
#UNLT2 <- LTabr(
#		nMx = MPnMx,
#		Age = c(0,1,seq(5,80,by=5)),
#		AgeInt = inferAgeIntAbr(vec = MPnMx),
#		axmethod = "UN",
#		Sex = "F")