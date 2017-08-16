
# Author: tim
###############################################################################

# TR: lots of approaches to a(0). Sometimes called 'separation factors'. These ought to be 
# collected, organized, and standardized here.

# these are the PAS versions:

#' Coale-Demeny a0 as function of m0, region, and sex
#' 
#' @description Coale-Demeny a0 from Manual X Table 164. This is just a rule of thumb. 
#' In this and some other older texts, a(0) is known as a 'separation factor'.
#' 
#' @details If sex is given as both, \code{"b"}, then female 
#' values are taken, per the PAS spreadsheet. This function is not vectorized. 
#' Formulas for North, South, and West are identical- only East is different. If \code{IMR} is not given,
#' then \code{M0} is used in its stead.
#' 
#' @param M0 numeric. Event exposure infant mortality rate
#' @param IMR numeric. Optional. q0, the death probability in first year of life, in case available separately.
#' @param Sex character. \code{"m"}, \code{"f"} or \code{"b"} for male, female, or both.
#' @param region character. \code{"n"}, \code{"e"}, \code{"s"} or \code{"w"} for North, East, South, or West.
#' 
#' @return a0 the average age at death in the first year of life.
#' @export
#' @references \insertRef{united1983manual}{DemoTools}
#' @examples
#' m0 <- seq(.001,.2,by=.001)
#' \dontrun{
#' plot(m0, sapply(m0, geta0CD, Sex = "m", region = "e"), ylab = "a0", 
#' 		type = 'l', ylim = c(0,.36), lty = 2, col = "blue")
#' lines(m0,sapply(m0, geta0CD, Sex = "m", region = "w"), col = "blue")
#' lines(m0,sapply(m0, geta0CD, Sex = "f", region = "e"), lty = 2, col = "red")
#' lines(m0,sapply(m0, geta0CD, Sex = "f", region = "w"), col = "red")
#' text(.15, geta0CD(.15,Sex = "m", region = "e"),"males E",font=2)
#' text(.15, geta0CD(.15,Sex = "m", region = "w"),"males N,W,S",font=2)
#' text(.15, geta0CD(.15,Sex = "f", region = "e"),"females E",font=2)
#' text(.15, geta0CD(.15,Sex = "f", region = "w"),"females N,W,S",font=2)
#' 
#' # compare with the Preston approximation
#' # constants identical after m0 = .107
#' m0 <- seq(.001,.107,by=.001)
#' a0CDm0 <- sapply(m0, geta0CD, Sex = "m", region = "w")
#' a0CDpr <- 0.045 + 2.684 * m0
#' plot(m0, a0CDm0, type = 'l', lty = 2, col = "red")
#' lines(m0, a0CDpr)
#' plot(m0, (a0CDm0 - a0CDpr) * 365, main = "difference (days)", ylab = "days")
#'}
# this is called a separation factor in the spreadsheet?
# separate estimate of IMR optional
# Markus: see if the model lifetables in refs/UN_1982... follow this rule of thumb for age 0.
# You may need to try giving q0 as IMR, or else M0 as M0 to the function, not sure.
geta0CD <- function(M0, IMR = NA, Sex = "M", region = "W"){
	# sex can be "m", "f", or "b"
	# region can be "n","e","s","w",or
	Sex    <- tolower(Sex)
	region <- tolower(region)
	
	Age0Const <- matrix(c(0.33, 	0.35, 	0.3500, 
					0.33, 	0.35, 	0.3500, 
					0.29, 	0.31, 	0.3100,
					0.33, 	0.35, 	0.3500 
			), ncol = 3, byrow = TRUE,
			dimnames = list(c("w", "n", "e", "s"), c("m","f","b")))
	
	Intercept <- matrix(c(0.0425, 	0.05, 	0.0500, 
					0.0425, 	0.05, 	0.0500, 
					0.0025, 	0.01, 	0.0100, 
					0.0425, 	0.05, 	0.0500
			), ncol = 3, byrow = TRUE,
			dimnames = list(c("w", "n", "e", "s"), c("m","f","b")))
	Slope        <- c(2.875, 	3.000, 	3.0000 )
	names(Slope) <- c("m","f","b")
	
	Alpha <- Intercept[region, Sex]
	Beta  <- Slope[Sex]
	# IMR optional here, use approximation
	if (missing(IMR) | is.na(IMR)){
		# formula from PAS LTPOPDTH
		a   <- M0 * Beta
		b   <- 1 + M0 * (1 - Alpha)
		IMR <- (b - sqrt(b^2 - 4 * a * M0)) / (2 * a)		
	}
	ifelse(IMR > .1, Age0Const[region, Sex], {Alpha + Beta * IMR})
}
# separate estimate of IMR optional
# TR: I think it's funny that a1-4 doesn't depend at all on m1-4

#' Coale-Demeny 4a1 as function of m0, region, and sex
#' 
#' @description Coale-Demeny 4a1. This is just a rule of thumb. 
#' In this and some other older texts, 4a1 is known as a 'separation factor'. These coefficients
#' were pulled from the PAS spreadsheets \code{LTPOPDTH.XLS} and not located in the original
#' Manual X.
#' 
#' @details If sex is given as both, \code{"b"}, then female 
#' values are taken, per the PAS spreadsheet. This function is not vectorized. 
#' If \code{IMR} is not given, then \code{M0} is used in its stead.
#' 
#' @param M0 numeric. Event exposure infant mortality rate
#' @param IMR numeric. Optional. q0, the death probability in first year of life, in case available separately.
#' @param Sex character. \code{"m"}, \code{"f"} or \code{"b"} for male, female, or both.
#' @param region character. \code{"n"}, \code{"e"}, \code{"s"} or \code{"w"} for North, East, South, or West.
#' 
#' @return a0 the average age at death in the first year of life.
#' @export
#' @references \insertRef{united1983manual}{DemoTools}
#' @examples
#' m0 <- seq(.001,.2,by=.001)
#' \dontrun{
#' plot(m0, sapply(m0, geta1_4CD, Sex = "m", region = "e"), ylab = "4a1", 
#' 		type = 'l', ylim = c(1,2), lty = 2, col = "blue")
#' lines(m0,sapply(m0, geta1_4CD, Sex = "m", region = "w"), col = "blue")
#' lines(m0,sapply(m0, geta1_4CD, Sex = "m", region = "n"), col = "blue", lty = "8383",lwd=2)
#' lines(m0,sapply(m0, geta1_4CD, Sex = "m", region = "s"), col = "blue", lty = "6464",lwd=2)
#' lines(m0,sapply(m0, geta1_4CD, Sex = "f", region = "e"), lty = 2, col = "red")
#' lines(m0,sapply(m0, geta1_4CD, Sex = "f", region = "w"), col = "red")
#' lines(m0,sapply(m0, geta1_4CD, Sex = "f", region = "n"), col = "red", lty = "8383",lwd=2)
#' lines(m0,sapply(m0, geta1_4CD, Sex = "f", region = "s"), col = "red", lty = "6464",lwd=2)
#' 
#' text(.05, geta1_4CD(.05,Sex = "m", region = "e"),"males E",font=2,pos=4)
#' text(.05, geta1_4CD(.05,Sex = "m", region = "w"),"males W",font=2,pos=4)
#' text(.05, geta1_4CD(.05,Sex = "m", region = "s"),"males S",font=2,pos=4)
#' text(.05, geta1_4CD(.05,Sex = "m", region = "n"),"males N",font=2,pos=4)
#' 
#' text(0, geta1_4CD(.01,Sex = "f", region = "e"),"females E",font=2,pos=4)
#' text(0, geta1_4CD(.01,Sex = "f", region = "w"),"females W",font=2,pos=4)
#' text(0, geta1_4CD(.01,Sex = "f", region = "s"),"females S",font=2,pos=4)
#' text(0, geta1_4CD(.01,Sex = "f", region = "n"),"females N",font=2,pos=4)
#' 
#' }
# Markus: see if the model lifetables in refs/UN_1982... follow this rule of thumb for age 1-4.
# You may need to try giving q0 as IMR, or else M0 as M0 to the function, not sure.
geta1_4CD <- function(M0, IMR = NA, Sex = "M", region = "W"){
	# sex can be "m", "f", or "b"
	# region can be "n","e","s","w",or
	Sex    <- tolower(Sex)
	region <- tolower(region)
	Age1_4Const <- matrix(c(1.352, 	1.361, 	1.3610, 
					1.558, 	1.570, 	1.5700, 
					1.313, 	1.324, 	1.3240, 
					1.240, 	1.239, 	1.2390
			), ncol = 3, byrow = TRUE,
			dimnames = list(c("w", "n", "e", "s"), c("m","f","b")))
	
	Intercept <- matrix(c(1.653, 	1.524, 	1.5240, 
					1.859, 	1.733, 	1.7330, 
					1.614, 	1.487, 	1.4870, 
					1.541, 	1.402, 	1.4020 
			), ncol = 3, byrow = TRUE,
			dimnames = list(c("w", "n", "e", "s"), c("m","f","b")))
	Slope        <- c(3.013, 	1.627, 	1.6270 )
	names(Slope) <- c("m","f","b")
	if (missing(IMR) | is.na(IMR)){
		a0 <- geta0CD(M0, IMR = NA, Sex = Sex, region = region)
		IMR <- mxax2qx(nMx = M0, nax = a0, AgeInt = 1, closeout = FALSE, IMR = NA)
	}
	ifelse(IMR > .1, Age1_4Const[region, Sex], Intercept[region, Sex] - Slope[Sex] * IMR)
}


#' PAS a(x) rule of thumb
#' 
#' @description ax is calculated following the Coale-Demeny rules for ages 0 and 1-4, and assumes interval midpoints in higher ages. 
#' This is just a rule of thumb. This procedure is as found in the PAS spreadsheet \code{LTPOPDTH.XLS}.
#' 
#' @details If sex is given as both, \code{"b"}, then female 
#' values are taken for a0 and 4a1, per the PAS spreadsheet. If IMR is not given, the M0 is used in its 
#' stead for ages < 5. This function is not vectorized. ax closeout assumes constant mortality hazard in the open age group.
#' 
#' @param M0 numeric. Event exposure infant mortality rate
#' @param IMR numeric. Optional. q0, the death probability in first year of life, in case available separately.
#' @param Sex character. \code{"m"}, \code{"f"} or \code{"b"} for male, female, or both.
#' @param region character. \code{"n"}, \code{"e"}, \code{"s"} or \code{"w"} for North, East, South, or West.
#' @param OAG logical (default \code{TRUE}). Is the last element of \code{nMx} the open age group?
#' 
#' @return nAx average contribution to exposure of those dying in the interval.
#' @export
#' @references \insertRef{united1983manual}{DemoTools}
# Markus: see if the model lifetables in refs/UN_1982... follow this rule of thumb for age 0.
# You may need to try giving q0 as IMR, or else M0 as M0 to the function, not sure.
axPAS <- function(nMx, AgeInt, IMR = NA, Sex = "M", region = "W", OAG = TRUE){
	# sex can be "m", "f", or "b"
	# region can be "n","e","s","w",or
	Sex    <- tolower(Sex)
	region <- tolower(region)
	
	N      <- length(nMx)
	ax     <- AgeInt / 2
	
	ax[1]  <- geta0CD(M0 = nMx[1], IMR = IMR, Sex = Sex, region = region)
	ax[2]  <- geta1_4CD(M0 = nMx[1], IMR = IMR, Sex = Sex, region = region)
	ax[N]  <- ifelse(OAG, 1 / nMx[N], ax[N])
	ax
}

#' UN version of the Greville formula for a(x) from M(x)
#' 
#' @description The UN ax formula uses Coale-Demeny for ages 0, and 1-4, values of 2.5 
#' for ages 5-9 and 10-14, and the Greville formula thereafter. In the original sources 
#' these are referred to as separation factors.
#' 
#' @details ax for age 0 and age group 1-4 are based on Coale-Demeny q0-based lookup tables.
#' We use an approximation to get from m0 to q0 for the sake of generating a0 and 4a1. 
#' The final ax value is closed out as if the final Mx value were constant thereafter, 
#' which is a common lifetable closeout choice. Age groups must be standard abridged, 
#' and we do not check age groups here! 
#' 
#' @param nMx numeric vector. Mortality rates in standard abridged age groups.
#' @param nqx numeric vector. Age specific death probabilities in standard abridged age groups.
#' @param nlx numeric vector. Lifetable survivorship in standard abridged age groups.
#' @param IMR numeric infant death probability. Optional if using \code{nMx}, not required otherwise.
#' @param Sex character. \code{"m"}, \code{"f"} or \code{"b"} for male, female, or both.
#' @param region character. \code{"n"}, \code{"e"}, \code{"s"} or \code{"w"} for North, East, South, or West.
#' @param mod
#' 
#' @return ax numeric vector of a(x) the average time spent in the interval of those that die in the interval.
#' @export
#' 
#' @references
#' \insertRef{greville1977short}{DemoTools}
#' \insertRef{un1982model}{DemoTools}
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertref{mortpak1988}{DemoTools}

ax.greville.mortpak <- function(nMx, nqx, lx, IMR = NA, Sex = "m", region = "W", mod = TRUE){
	Sex     <- tolower(Sex)
	region  <- tolower(region)
	DBL_MIN <- .Machine$double.xmin
	stopifnot(!missing(nMx) | !missing(nqx) | !missing(lx))
	# sort out arguments:
	if (missing(nqx) & !missing(lx)){
		nqx <- lx2dx(lx) / lx
	}
	# now we have either qx or mx
	
	if (missing(nqx) & !missing(nMx)){
		a0     <- geta0CD(M0 = nMx[1], IMR = IMR, Sex = Sex, region = region)
		a1_4   <- geta1_4CD(M0 = nMx[1], IMR = IMR, Sex = Sex, region = region)
		qind   <- FALSE
	}
	if (missing(nMx) & !missing(nqx)){
		a0     <- geta0CD(M0 = NA, IMR = nqx[1], Sex = Sex, region = region)
		a1_4   <- geta1_4CD(M0 = NA, IMR = nqx[1], Sex = Sex, region = region)
		# just call it nMx for rest of calcs:
		nMx    <- nqx
		qind   <- TRUE
	}
	
	# some setup 
	N      <- length(nMx)
	AgeInt <- inferAgeIntAbr(vec=nMx)
	
	# use ages rather than index positions to select.
	Age    <- AgefromAgeInt(AgeInt)
	# default midpoints to overwrite
	ax     <- AgeInt / 2
	
	for (i in 3:(N - 1)) {
		## Mortpak LIFETB for age 5-9 and 10-14
		# ax[j] = 2.5 
		## for ages 15-19 onward
		## AK <- log(QxMx[j+1]/QxMx[j-1])/10
		## ax[j] <- 2.5 - (25.0/12.0) * (QxMx[j] - AK)
		
		## improved Greville formula for adolescent ages 5-9 and 10-14
		## Let the three successive lengths be n1, n2 and n3, the formula for 5a5 is:
		## ax[i] = 2.5 - (25 / 12) * (mx[i] - log(mx[i + 1] / mx[i-1])/(n1/2+n2+n3/2))
		## for age 5-9, coefficient should be 1/9.5, because age group 1-4 has only 4 ages (not 5), while the other 5-year age group are 1/10
		## ax[i] = 2.5 - (25 / 12) * (mx[i] - (1/9.5)* log(mx[i + 1] / mx[i-1]))
		## Age 20-25, ..., 95-99
		## Greville (based on Mortpak LIFETB) for other ages, new implementation
		
		# front term
		Af     <- (AgeInt[i] / 2) - (AgeInt[i]^2 / 12)
		# back term
		Ab     <- 1 / (AgeInt[i - 1] / 2 + AgeInt[i] + AgeInt[i + 1] / 2)
		# subtract
	    if (i < (N - 1)){
		# N-1 uses K from N-2...
		   K      <- rlog(max(nMx[i + 1] / max(nMx[i - 1], DBL_MIN))) 
	    }
		# main formula
		ax[i]  <- Af * (nMx[i] - K) * Ab

		## add constraint at older ages (in Mortpak and bayesPop)
		## 0.97 = 1-5*exp(-5)/(1-exp(-5)), for constant mu=1, Kannisto assumption  
		## (Mortpak used 1 instead of 0.97).
		if (Age[i] > 35 && ax[i] < 0.97) {
			ax[i] <- 0.97
		}
		
		## Extra condition based on Mortpak LIFETB for age 65 onward
		# TR: why .8a[x-1] only for qx?
        tmp   <- 0.8 * ax[i - 1]
		ax[i] <- ifelse(qind & Age[i] >= 65 & ax[i] < tmp, tmp, ax[i])

	}

	if (!mod){
		ax[3:4] <- 2.5
	}
    # closeout
	aomega         <- max(c(1 / nMx[N], .97))
	aomega         <- ifelse(qind, max(aomega, .8 * ax[N - 1]), aomega)
	ax[N]          <- aomega
	# 
	ax
}

#' UN a(x) estimates from either M(x), q(x), or both
#' 
#' @description The UN ax formula uses Coale-Demeny for ages 0, and 1-4, values of 2.5 
#' for ages 5-9 and 10-14, and the Greville formula for higher ages. In the original sources 
#' these are referred to as separation factors.
#' 
#' @details ax for age 0 and age group 1-4 are based on Coale-Demeny q0-based lookup tables. If the main 
#' input is \code{nMx}, and if \code{IMR} is not given, we first approximate q0 for the CD formula 
#' before applying the formula.
#' The final ax value is closed out as the reciprocal of the final m(x), which is standard lifetable practice, but . For nMx inputs 
#' this method is rather direct, but for qx or lx inputs it is iterative. Age groups must be standard abridged, 
#' and we do not check age groups here!
#' 
#' @param nMx numeric vector. Mortality rates in standard abridged age groups.
#' @param nqx numeric vector. Age specific death probabilities in standard abridged age groups.
#' @param nlx numeric vector. Lifetable survivorship in standard abridged age groups.
#' @param IMR numeric. Infant death probability, (q0). Optional.
#' @param Sex character. \code{"m"}, \code{"f"} or \code{"b"} for male, female, or both.
#' @param region character. \code{"n"}, \code{"e"}, \code{"s"} or \code{"w"} for North, East, South, or West.
#' @param tol the tolerance for the qx-based iterative method (default \code{.Machine$double.eps}).
#' @param maxit the maximum numbder of iterations for the qx-based iterative method (default 1000).
#' @param mod logical (default \code{TRUE}). Use Gerland's modification for ages 5-14?
#' 
#' @return ax numeric vector of a(x) the average time spent in the interval of those that die in the interval.
#' @export
#' 
#' @references
#' \insertRef{greville1977short}{DemoTools}
#' \insertRef{un1982model}{DemoTools}
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertref{mortpak1988}{DemoTools}
#' @examples 
#' # example data from UN 1982 Model Life Tables for Developing Countries.
#' # first Latin American model table for males (p. 34).
#' Mx <- c(.23669,.04672,.00982,.00511,.00697,.01036,.01169,
#' 		.01332,.01528,.01757,.02092,.02517,.03225,.04241,.06056,
#' 		.08574,.11840,.16226,.23745)
#' ax <- c(0.330,1.352,2.500,2.500,2.633,2.586,2.528,2.528,
#' 		2.526,2.529,2.531,2.538,2.542,2.543,2.520,2.461,2.386,2.295,4.211)
#' 
#' AgeInt     <- inferAgeIntAbr(vec = Mx)
#' 
#' nAx1       <- axUN(nMx = Mx, 
#' 		              AgeInt = AgeInt, 
#' 		              Sex = "m",
#' 		              mod = FALSE)
#' nAx2       <- axUN(nMx = Mx, 
#' 		              AgeInt = AgeInt, 
#' 		              Sex = "m",
#' 		              mod = TRUE)
#' # this is acceptable...
#' round(nAx2,3) - ax # only different in ages 5-9 and 10-14
#' # default unit test...
#' stopifnot(all(round(nAx1,3) - ax == 0)) # spot on
axUN <- function(
		nMx,
		nqx, 
		lx, 
		IMR = NA, 
		AgeInt, 
		Sex = "m", 
		region = "W", 
		tol = .Machine$double.eps, 
		maxit = 1e3, 
		kludge = FALSE){
	stopifnot(!missing(nqx) | !missing(nMx))
	smsq    <- 99999
	Sex     <- tolower(Sex)
	region  <- tolower(region)
	
	
	if (missing(nqx) & !missing(lx)){
		nqx <- lx2dx(lx) / lx
	}
	# now we have either nqx or nMx
	
	if (missing(nqx) & !missing(nMx)){
# UN (1982) p 31 
# http://www.un.org/esa/population/publications/Model_Life_Tables/Model_Life_Tables.htm
#		For ages 15 and over, the expression for nQx is derived
#		from Greville" as ,nax, = 2.5 - (25/12) (nmx) - k), where
#		k = 1/10 log(nmx+5/nmx-5). For ages 5 and 10, nQx = 2.5
#		and for ages under 5, nQx values from the Coale and
#		Demeny West region relationships are used." 
		
		axi <- ax.greville.mortpak(nMx = nMx, IMR = IMR, Sex = Sex, region = region, mod = mod)
	}
	if (!missing(nqx) & missing(nMx)){
# UN (1982) p 31 
# http://www.un.org/esa/population/publications/Model_Life_Tables/Model_Life_Tables.htm
#		With nqx as input, the procedure is identical, except
#		that an iterative procedure is used to find the nmx and nqx
#		values consistent with the given nqx and with the Greville
#		expression. 
		
		axi <- ax.greville.mortpak(nqx = nqx, Sex = Sex, region = region, mod = mod)
		
		mxi <- mx.from.qx(nqx, ax = axi)
		for (i in 1:maxit) {
			mxi   <- qx2mx(nqx, axi)
			axi   <- ax.greville.mortpak(nMx = mxi, IMR = nqx[1], Sex = Sex, region = region, mod = mod)
			qxnew <- mx2qx(mxi, axi)
			smsq  <- sum((qxnew - nqx)^2)
			if (smsq < tol){
				break
			}
		}
		# no need for approximate a0 and 4a1 values:
		
	}
	# if both given, then we have ax via identity:
	if (!missing(nqx) & !missing(nMx)){
		axi <- qxmx2ax(nqx = nqx, nmx = nMx, AgeInt = inferAgeIntAbr(vec = nMx))
	}
	
	# if mx, qx, or both are given, then by now we have ax
	axi
}



