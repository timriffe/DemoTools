# Author: tim
###############################################################################

#' Calculate an abridged-age lifetable.
#' @description Given vectors for Deaths and Exposures, or Mx, or qx, or lx,
#' calculate a full abridged lifetable.
#'
#' @details The main variations here are in the treatment of \code{nAx}.
#' In all cases, the lifetable is extended and closed out using one from a
#' selection of mortality age extrapolation methods implemented in the
#'  \code{MortalityLaws} package rather than the common practice of taking
#'  the inverse of the final \code{nMx} value (if it's an open interval).
#'  For this, a desired open age must be specified, defaulting to the present
#'  open age group, but which can not exceed 110 in the present implementation.
#'  By default, the extrapolation model is fit to ages 60 and higher, but you
#'  can control this using the \code{extrapFit} argument (give the vector of ages,
#'  which must be a subset of \code{Age}). By default extrapolated values are used
#'  starting with the input open age, but you can lower this age using the
#'  \code{extrapFrom} argument. The \code{Sx} output column (survivor ratios)
#'  is aligned with the other columns in all 5-year age groups, but note the
#'  first two values have a slightly different age-interval interpretation:
#'  In Age 0, the interpretation is survival from birth until interval 0-4.
#'  In Age 1, it is survival from 0-4 into 5-9. Therafter the age groups align.
#'  This column is required for population projections.
#'
#' @param Deaths numeric. Vector of death counts in abridged age classes.
#' @param Exposures numeric. Vector of population exposures in abridged age classes.
#' @param nMx numeric. Vector of mortality rates in abridged age classes.
#' @param nqx numeric. Vector of conditional death probabilities in abridged age classes.
#' @param lx numeric. Vector of lifetable survivorship at abridged ages.
#' @param AgeInt integer. Vector of age class widths. Default \code{inferAgeIntAbr(Age = Age)}.
#' @param radix numeric. Lifetable radix, \ifelse{html}{\out{l<sub>0}}{\eqn{l_0}}. Default 100000.
#' @param axmethod character. Either \code{"pas"} or \code{"un"}.
#' @param a0rule character. Either \code{"ak"} (default) or \code{"cd"}.
#' @param Sex character. Either male \code{"m"}, female \code{"f"}, or both \code{"b"}.
#' @param region character. North, East, South, or West: code{"n"}, code{"e"}, code{"s"}, code{"w"}. Default code{"w"}.
#' @param IMR numeric. Infant mortality rate \ifelse{html}{\out{q<sub>0}}{\eqn{q_0}}, in case available and \code{nqx} is not specified. Default \code{NA}.
#' @param mod logical. If \code{"un"} specified for \code{axmethod}, whether or not to use Nan Li's modification for ages 5-14. Default \code{TRUE}.
#' @param SRB the sex ratio at birth (boys / girls), detault 1.05
#' @param OAnew integer. Desired open age group (5-year ages only). Default \code{max(Age)}. If higher then rates are extrapolated.
#' @param OAG logical. Whether or not the last element of \code{nMx} (or \code{nqx} or \code{lx}) is an open age group. Default \code{TRUE}.
#' @param extrapLaw character. If extrapolating, which parametric mortality law should be invoked? Options include
#'   \code{"Kannisto", "Kannisto_Makeham", "Makeham", "Gompertz", "GGompertz", "Beard",	"Beard_Makeham", "Quadratic"}. Default \code{"Kannisto"}. See details.
#' @inheritParams lt_a_closeout
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
#'   \item{Sx}{numeric. Survivor ratios in uniform 5-year age groups.}
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
#' PASLT <- lt_abridged(Deaths = Deaths,
#' 		Exposures = Exposures,
#' 		Age = Age,
#' 		AgeInt =AgeInt,
#' 		axmethod = "pas",
#' 		IMR = .1,
#' 		region = "n",
#' 		Sex = "m",
#'		a0rule ="cd")
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
#' UNLT1 <- lt_abridged(nMx = Mx,
#' 		Age = c(0,1,seq(5,85,by=5)),
#' 		AgeInt = AgeInt,
#' 		axmethod = "un",
#' 		Sex = "m",
#' 		mod = FALSE)
#' UNLT2 <- lt_abridged(nMx = Mx,
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
#'  MPnMx <- c(0.12846,0.02477,0.00603,0.0034,
#'  		0.00417,0.00513,0.00581,0.00645,0.00725,
#'  		0.00813,0.00913,0.01199,0.01647,
#'  		0.0256,0.04047,0.06624,0.10638,0.19611)
#'  Age <- c(0,1,seq(5,80,by=5))
#'  AgeInt <- age2int(Age,OAvalue = 5)
#'  MPexcheck <- c(49.997,55.675,57.245,53.921,
#'  		49.803,45.799,41.922,38.084,34.249,
#'  		30.420,26.578,22.701,18.945,
#'  		15.349,12.095,9.240,6.903,5.099)
#'
#'  # First with lifetable extention to 100
#'  MP_UNLT100 <- lt_abridged(
#'  		nMx = MPnMx,
#'  		Age = Age,
#'  		AgeInt = AgeInt,
#'  		axmethod = "un",
#'  		Sex = "f",
#'  		mod = FALSE,
#'  		OAnew = 100,
#'  		a0rule ="cd")
#' #'
#' #' # lifetable to original open age group
#'  MP_UNLT80 <- lt_abridged(
#'  		nMx = MPnMx,
#'  		Age = Age,
#'  		AgeInt = AgeInt,
#'  		axmethod = "un",
#'  		Sex = "f",
#'  		mod = FALSE,
#'  		OAnew = 80,
#'  		a0rule ="cd")
#'
#' # same, but truncated at 60
#' MP_UNLT60 <- lt_abridged(
#' 		nMx = MPnMx,
#' 		Age = Age,
#' 		AgeInt = AgeInt,
#' 		axmethod = "un",
#' 		Sex = "f",
#' 		mod = FALSE,
#' 		OAnew = 60,
#'  		a0rule ="cd")

lt_abridged <- function(Deaths = NULL,
                  Exposures = NULL,
                  nMx = NULL,
                  nqx = NULL,
                  lx = NULL,
                  Age,
                  AgeInt = age2int(Age = Age, OAvalue = 5),
                  radix = 1e5,
                  # redesign backwards from the ax args. Possibly add some options here?
                  axmethod = "pas", 
                  a0rule = "ak", 
                  Sex = "m", 
                  IMR = NA, 
                  region = "w", 
                  mod = TRUE,
                  SRB = 1.05,
                  OAG = TRUE,
                  OAnew = max(Age),
                  extrapLaw = "kannisto",
                  extrapFrom = max(Age),
                  extrapFit = Age[Age >= 60 & ifelse(OAG, Age < max(Age), TRUE)],
                  ...) {
  axmethod <- match.arg(axmethod, choices = c("pas","un"))
  Sex      <- match.arg(Sex, choices = c("m","f","b"))
  a0rule   <- match.arg(a0rule, choices = c("ak","cd"))
  extrapLaw      <- match.arg(extrapLaw, choices = c("kannisto",
                                                     "kannisto_makeham",
                                                     "makeham",
                                                     "gompertz",
                                                     "ggompertz",
                                                     "beard",
                                                     "beard_makeham",
                                                     "quadratic"
  ))
  region   <- match.arg(region, choices = c("w","n","s","e"))
  # ages must be abridged.
  stopifnot(is_abridged(Age))

  # now overwriting raw nMx is allowed by lowering this
  # arbitrary lower bound to accept the fitted model. Really
  # this functionality is intended for extrapolation and not
  # model overwriting of rates.
  stopifnot(extrapFrom <= max(Age))
  # need to make it possible to start w (D,E), M, q or l...

  # TR: make sure IMR propagates
  imr_flag <- !is.na(IMR)
  # and use more consistent flags
  mxflag   <- !is.null(nMx)
  qxflag   <- !is.null(nqx)
  # 1) if lx given but not qx:
  if ((!qxflag) & (!is.null(lx))) {
    nqx          <- lt_id_l_d(lx) / lx
    nqx[1]       <- ifelse(imr_flag, IMR, nqx[1])
    qxflag       <- TRUE
  }
  # 2) if still no nqx then make sure we have or can get nMx
  if ((!qxflag) & (!mxflag)) {
    stopifnot((!is.null(Deaths)) & (!is.null(Exposures)))
    nMx          <- Deaths / Exposures
    mxflag       <- TRUE
  }

  # axmethod       <- tolower(axmethod)
  # Sex            <- tolower(Sex)
  # region         <- tolower(region)
  # extrapLaw      <- tolower(extrapLaw)

  # take care of ax first, two ways presently
  if (!mxflag) {
    nqx[1]       <- ifelse(imr_flag, IMR, nqx[1])

    # TR: expedient hack
    nqx[nqx > 1] <- 1

    # TR: note q0 likely different from IMR. In this case use IMR
    nAx          <- lt_id_morq_a(
                      nqx = nqx,
                      axmethod = axmethod,
                      Age = Age,
                      AgeInt = AgeInt,
                      a0rule = a0rule,
                      Sex = Sex,
                      region = region,
                      OAG = OAG,
                      mod = mod,
                      IMR = IMR,
                      SRB = SRB,
                      extrapLaw = extrapLaw,
                      extrapFrom = extrapFrom,
                      extrapFit = extrapFit
                      )
  } else {
    nAx          <- lt_id_morq_a(
                      nMx = nMx,
                      axmethod = axmethod,
                      Age = Age,
                      AgeInt = AgeInt,
                      a0rule = a0rule,
                      Sex = Sex,
                      region = region,
                      OAG = OAG,
                      mod = mod,
                      IMR = IMR,
                      SRB = SRB,
                      extrapLaw = extrapLaw,
                      extrapFrom = extrapFrom,
                      extrapFit = extrapFit)
  }
  # TR, these nAx ought to turn out to be the same...

  #	# as of here we have nAx either way. And we have either mx or qx.
  # TR: 31-12-2019 comment out, code chunk appears pointless, since
  # qx redefined
  # if (missing(nqx)) {
  #   nqx          <- lt_id_ma_q(
  #                     nMx = nMx,
  #                     nax = nAx,
  #                     AgeInt = AgeInt,
  #                     closeout = TRUE,
  #                     IMR = IMR)
  # }
  # TR: in the case that nMx is missing, then we must have nqx by now
  if (!mxflag){
    # TR note, if OAG, then final nMx is NA here!
    nMx          <- lt_id_qa_m(nqx = nqx,
                               nax = nAx,
                               AgeInt = AgeInt)
  }
  # now we have all three, [mx,ax,qx] guaranteed.

  OA             <- max(Age)
  # TR: save for later, in case OAG preserved
  # TR: 5-1-2020, having doubts re this, can also
  # back out momega from extrapolation, and this would
  # handle closeout agreement in case of nqx inputs
  if (OAG & OAnew == OA & mxflag) {
    momega       <- nMx[length(nMx)]
  }
  # --------------------------------
  # begin extrapolation:
  # TR: 13 Oct 2018. always extrapolate to 130 no matter what,
  # then truncate to OAnew in all cases. This will ensure more robust closeouts
  # and an e(x) that doesn't depend on OAnew. 130 is used similarly by HMD.
  x_extr         <- seq(extrapFrom, 130, by = 5)

  Mxnew          <- lt_rule_m_extrapolate(
                      x = Age,
                      mx = nMx,
                      x_fit = extrapFit,
                      x_extr = x_extr,
                      law = extrapLaw,
                      ...)

  nMxext         <- Mxnew$values
  Age2           <- names2age(nMxext)

  keepi          <- Age2 < extrapFrom
  nMxext[keepi]  <- nMx[Age < extrapFrom]
  nMx            <- nMxext
  Age            <- Age2
  AgeInt         <- age2int(
                      Age,
                      OAG = TRUE,
                      OAvalue = max(AgeInt, na.rm = TRUE))
  # redo ax and qx for extended ages
  nAx            <- lt_id_morq_a(
                      nMx = nMx,
                      axmethod = axmethod,
                      Age = Age,
                      AgeInt = AgeInt,
                      a0rule = a0rule,
                      Sex = Sex,
                      region = region,
                      OAG = TRUE,
                      mod = mod,
                      IMR = IMR,
                      SRB = SRB)

  nqx            <- lt_id_ma_q(
                      nMx = nMx,
                      nax = nAx,
                      AgeInt = AgeInt,
                      closeout = TRUE,
                      IMR = IMR)

  # one more step: make sure M0 agrees with
  # q0 and a0, which may have been determined by IMR parameter
  nMx[1] <- lt_id_qa_m(nqx = nqx[1], nax = nAx[1], AgeInt = 1)


  # end extrapolation
  # ---------------------------------

  # TR: the lifetable is the shortest part of this code!
  lx            <- lt_id_q_l(nqx, radix = radix)
  ndx           <- lt_id_l_d(lx)
  nLx           <- lt_id_lda_L(lx = lx,
                               ndx = ndx,
                               nax = nAx,
                               AgeInt = AgeInt)
  Tx            <- lt_id_L_T(nLx)
  ex            <- Tx / lx

  # TR: now cut down due to closeout method (added 11 Oct 2018)
  ind           <- Age <= OAnew
  Age           <- Age[ind]
  AgeInt        <- AgeInt[ind]
  nAx           <- nAx[ind]
  nMx           <- nMx[ind]
  nqx           <- nqx[ind]

  lx            <- lx[ind]
  # recalc dx from chopped lx
  ndx           <- lt_id_l_d(lx)
  nLx           <- nLx[ind]
  Tx            <- Tx[ind]
  ex            <- ex[ind]

  # some closeout considerations
  N          <- length(nqx)
  nqx[N]     <- 1
  nLx[N]     <- Tx[N]
  nAx[N]     <- ex[N]
  AgeInt[N]  <- NA
  # TR: https://github.com/timriffe/DemoTools/issues/83
  # (Added 3 Oct 2019)
  if (OAG) {
    if (OAnew == OA & mxflag) {
      nMx[N] <- momega
    } else {
      # Otherwise inner coherence
      nMx[N] <-  lx[N] / Tx[N]
    }
  } else {
    nMx[N]   <-  lx[N] / Tx[N]
  }

  Sx <- lt_id_Ll_S(nLx, lx, AgeInt, N = 5)
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
    Sx = Sx,
    Tx = Tx,
    ex = ex
  )
  return(out)
}
