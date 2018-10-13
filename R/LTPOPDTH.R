# TODO: introducing new OAG argument, but not being handled well by UN ax method.
# when FALSE it appears to still be doing an abacus extrapolation which goes haywire
# with patrick's qx test case.
# Author: tim
###############################################################################

#' Calculate an abridged-age lifetable.
#' @description Given vectors for Deaths and Exposures, or Mx, or qx, or lx, calculate a full abridged lifetable.
#' 
#' @details The main variations here are in the treatment of nAx. In all cases, the lifetable is extended and closed out 
#' using one from a selection of methods implemented in the \code{MortalityLaws} package. 
#' For this, a desired open age must be specified, defaulting to the present open age group, but which can not exceed 100 in the present implementation.
#' The ax method beyond and including the original open age group in the data follows the Mortpak (UN) method only.
#' 
#' @param Deaths numeric. Vector of death counts in abridged age classes.
#' @param Exposures numeric. Vector of population exposures in abridged age classes.
#' @param nMx numeric. Vector of mortality rates in abridged age classes.
#' @param nqx numeric. Vector of conditional death probabilities in abridged age classes.
#' @param lx numeric. Vector of lifetable survivorship at abridged ages.
#' @param Age integer. A vector of ages of the lower integer bound of the age classes.
#' @param AgeInt integer. Vector of age class widths. Ddefault \code{inferAgeIntAbr(Age = Age)}.
#' @param radix numeric. Lifetable radix, \ifelse{html}{\out{l<sub>0}}{\eqn{l_0}}. Default 100000.
#' @param axmethod character. Either \code{"pas"} or \code{"un"}. 
#' @param Sex character. Either male \code{"m"}, female \code{"f"}, or both \code{"b"}.   
#' @param region character. North, East, South, or West: code{"n"}, code{"e"}, code{"s"}, code{"w"}. Default code{"w"}.
#' @param IMR numeric. Infant mortality rate \ifelse{html}{\out{q<sub>0}}{\eqn{q_0}}, in case available and \code{nqx} is not specified. Default \code{NA}.
#' @param mod logical. If \code{"un"} specified for \code{axmethod}, whether or not to use Patrick Gerland's modification for ages 5-14. Default \code{TRUE}.
#' @param OAnew integer. Desired open age group (5-year ages only). Default \code{max(Age)}. If higher then rates are extrapolated.
#' @param extrapLaw character. If extrapolating, which parametric mortality law should be invoked? Options include  \code{"Kannisto","Makeham", "Gompertz", "Gamma-Gompertz", "Beard",
#'                    "Makeham-Beard", "Quadratic"}. Default \code{"Kannisto"}. See details.
#' @param extrapFrom integer. Age from which to impute extrapolated mortality.
#' @param extrapFit integer vector. Ages to include in model fitting. Defaults to all ages \code{>=60}.
#' @export 
#' @return Lifetable in data.frame with columns
#' \itemize{
#'   \item{Age}{integer. Lower bound of abridged age class},
#'   \item{AgeInt}{integer. Age class widths.}
#'   \item{nMx}{numeric. Age-specific central death rates.} 
#'   \item{nAx}{numeric. Average time spent in interval by those deceased in interval. } 
#'   \item{nqx}{numeric. Age-specific conditional death probabilities.} 
#'   \item{lx}{numeric. Lifetable survivorship} 
#'   \item{ndx}{numeric. Lifetable deaths distribution.} 
#'   \item{nLx}{numeric. Lifetable exposure.} 
#'   \item{Tx}{numeric. Lifetable total years left to live above age x.} 
#'   \item{ex}{numeric. Age-specific remaining life expectancy.}
#' }
#' @references 
#' \insertRef{greville1977short}{DemoTools}
#' \insertRef{un1982model}{DemoTools}
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertRef{mortpak1988}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#' 
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
#' PASLT <- LTabr(Deaths = Deaths, 
#' 		Exposures = Exposures, 
#' 		Age = Age,
#' 		AgeInt =AgeInt,
#' 		axmethod = "pas",
#' 		IMR = .1,
#' 		region = "n",
#' 		Sex = "m")
#' 
#' # examples based on UN 1982 (p. 34)
#' Mx <- c(.23669,.04672,.00982,.00511,.00697,.01036,.01169,
#' 		.01332,.01528,.01757,.02092,.02517,.03225,.04241,.06056,
#' 		.08574,.11840,.16226,.23745)
#' excheckUN <-  c(35.000,42.901,47.190,44.438,
#' 		40.523,36.868,33.691,30.567,27.500,24.485,21.504,18.599,
#' 		15.758,13.080,10.584,8.466,6.729,5.312,4.211)
#' AgeInt <- inferAgeIntAbr(vec = Mx)
#' 
#' # generate two variants: with and without PG's variants
#' # for ages 5-14
#' UNLT1 <- LTabr(nMx = Mx,
#' 		Age = c(0,1,seq(5,85,by=5)),
#' 		AgeInt = AgeInt,
#' 		axmethod = "un",
#' 		Sex = "m", 
#' 		mod = FALSE)
#' UNLT2 <- LTabr(nMx = Mx,
#' 		Age = c(0,1,seq(5,85,by=5)),
#' 		AgeInt = AgeInt,
#' 		axmethod = "un",
#' 		Sex = "m", 
#' 		mod = TRUE)
#' 		
#' #  \dontrun{
#' # 	 plot(UNLT2$ex - UNLT1$ex)
#' #  }
#' 
#' # a Mortpak unit test:
#' # data from  p. 82 United Nations (1988) Mortpak - ...
#' MPnMx <- c(0.12846,0.02477,0.00603,0.0034,
#' 		0.00417,0.00513,0.00581,0.00645,0.00725,
#' 		0.00813,0.00913,0.01199,0.01647,
#' 		0.0256,0.04047,0.06624,0.10638,0.19611)
#' 
#' MPexcheck <- c(49.997,55.675,57.245,53.921,
#' 		49.803,45.799,41.922,38.084,34.249,
#' 		30.420,26.578,22.701,18.945,
#' 		15.349,12.095,9.240,6.903,5.099)
#' 
#' # First with lifetable extention to 100
#' MP_UNLT100 <- LTabr(
#' 		nMx = MPnMx,
#' 		Age = c(0,1,seq(5,80,by=5)),
#' 		AgeInt = inferAgeIntAbr(vec = MPnMx),
#' 		axmethod = "un",
#' 		Sex = "f",
#' 		mod = FALSE,
#' 		OAnew = 100)
#' 
#' # lifetable to original open age group
#' MP_UNLT80 <- LTabr(
#' 		nMx = MPnMx,
#' 		Age = c(0,1,seq(5,80,by=5)),
#' 		AgeInt = inferAgeIntAbr(vec = MPnMx),
#' 		axmethod = "un",
#' 		Sex = "f",
#' 		mod = FALSE,
#' 		OAnew = 80)
#' 
#' # same, but truncated at 60
#' MP_UNLT60 <- LTabr(
#' 		nMx = MPnMx,
#' 		Age = c(0,1,seq(5,80,by=5)),
#' 		AgeInt = inferAgeIntAbr(vec = MPnMx),
#' 		axmethod = "un",
#' 		Sex = "f",
#' 		mod = FALSE,
#' 		OAnew = 60)


LTabr <- function(
		Deaths, 
		Exposures, 
		nMx,
		nqx,
		lx,
		Age,
		AgeInt = inferAgeIntAbr(Age = Age), 
		radix = 1e5,
		axmethod = "pas", 
		Sex = "m", 
		region = "w",
		IMR = NA,
		mod = TRUE,
		OAG = TRUE,
		OAnew = max(Age),
		extrapLaw = c("Kannisto", "Makeham", "Gompertz", "Gamma-Gompertz", "Beard",
                   "Makeham-Beard", "Quadratic")[1],
		extrapFrom = max(Age),
		extrapFit = Age[Age >= 60]){
	# this is a hard rule for now. May be relaxed if the 
	# ABACUS code is relaxed. For now mostly unmodified.
	#stopifnot(OAnew <= 100)
	# now overwriting raw nMx is allowed by lowering this
	# arbitrary lower bound to accept the fitted model. Really
	# this functionality is intended for extrapolation and not
	# model overwriting of rates.
	stopifnot(extrapFrom <= max(Age) & (extrapFrom > 50))
	# need to make it possible to start w (D,E), M, q or l...
	
	# 1) if lx given but not qx:
	if (missing(nqx) & !missing(lx)){
		nqx <- lx2dx(lx) / lx
	}
	# 2) if still no nqx then make sure we have or can get nMx
	if (missing(nqx) & missing(nMx)){
		stopifnot((!missing(Deaths)) & (!missing(Exposures)))
		nMx       <- Deaths / Exposures
	}
	
	# TR: what is this for?
#	if (missing(Deaths)){
#		Deaths    <- NULL
#	}
#	if (missing(Exposures)){
#		Exposures <- NULL
#	}
	
	axmethod      <- tolower(axmethod)
	Sex           <- tolower(Sex)
	region        <- tolower(region)
	extrapLaw     <- tolower(extrapLaw)
	# take care of ax first, two ways presently
	if (missing(nMx)){
		nAx <- mxorqx2ax(
				nqx = nqx, 
				axmethod = axmethod,
				AgeInt = AgeInt,
				Sex = Sex,
				region = region, 
				OAG = OAG, 
				mod = mod,
				IMR = IMR)
	} else {
		nAx <- mxorqx2ax(
				nMx = nMx, 
				axmethod = axmethod,
				AgeInt = AgeInt,
				Sex = Sex,
				region = region, 
				OAG = OAG, 
				mod = mod,
				IMR = IMR)
	}

#	# as of here we have nAx either way. And we have either mx or qx.
	
	if (missing(nqx)){
		nqx <- mxax2qx(nMx = nMx, nax = nAx, AgeInt = AgeInt, closeout = TRUE, IMR = IMR)
	}
	if (missing(nMx)){
		nMx <- qxax2mx(nqx = nqx, nax = nAx, AgeInt = AgeInt)
	}
	# now we have all three, [mx,ax,qx] guaranteed.
	
	# TR: begin new chunk 28 Feb, 2018 to allow extention 
	# and at a minimum give a better ex and ax estimate 
	# in open age group. But only if OA <= 100.
# TR: Oct 11, 2018: consider ability to overwrite mx using the 
# extrapolation law invoked. Ask PG if desired. Deprecate
# abacus code, too hard to maintain. Variety of laws now avail
# no guarantee of monotonicity in old age however. Consider extrapolating
# to OAnew + 10 and then truncating lifetable, which would give a better
# closeout ex and ax, but wouldn't affect the other columns.

# TR: 13 Oct 2018. NOTE switch to always extrapolate to 120 no matter what,
# then truncate to OAnew in all cases. This will ensure more robust closeouts
# and an e(x) that doesn't depend on OAnew.
	OA     <- max(Age)
	if (OA < OAnew){
		#x_fit  <- Age[Age >= 60]
		x_extr <- seq(extrapFrom, OAnew + 10, by = 5)
		Mxnew <- extra_mortality(
				x = Age, 
				mx = nMx, 
				x_fit = extrapFit,
				x_extr = x_extr,
				law = extrapLaw)

		nMxext <- Mxnew$values
		Age2   <- names2age(nMxext)
		# for now, no overwriting of earlier mx
	    # sum(Age <= extrapFrom)
		nMxext[Age2 <= extrapFrom] <- nMx[Age <= extrapFrom]
		nMx    <- nMxext
		# i.e. open age may or may not be treated as such as of here
	    # if extrapolating from now on OAG is true.
		AgeInt <- age2int(Age2,OAG=TRUE,OAvalue=max(AgeInt))
		# redo ax and qx for extended ages
		nAx <- mxorqx2ax(
				nMx = nMx, 
				axmethod = axmethod,
				AgeInt = AgeInt,
				Sex = Sex,
				region = region, 
				OAG = TRUE, 
				mod = mod)
		nqx <- mxax2qx(
				nMx = nMx, 
				nax = nAx, 
				AgeInt = AgeInt, 
				closeout = TRUE, 
				IMR = IMR)
		Age <- Age2
		
	}

	# TR: the lifetable is the shortest part of this code!
	lx  <- qx2lx(nqx, radix = radix)
	ndx <- lx2dx(lx)
	nLx <- lxdxax2Lx(lx = lx, ndx = ndx, nax = nAx, AgeInt = AgeInt)
	Tx  <- Lx2Tx(nLx)
	ex  <- Tx / lx
	
	# TR: now cut down due to closeout method (added 11 Oct 2018)
	ind <- Age <= OAnew
	Age <- Age[ind]
	AgeInt <- AgeInt[ind]
	nAx <- nAx[ind]
	nMx <- nMx[ind]
	nqx <- nqx[ind]
	nqx[length(nqx)] <- 1
	lx  <- lx[ind]
	ndx <- lx2dx(lx)
	nLx <- nLx[ind]
	Tx  <- Tx[ind]
	ex  <- ex[ind]
	
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


