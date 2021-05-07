# A single age lifetable that depends only on mx, first draft, not tested yet.

#' calculate a single age lifetable
#' @description Fuller description forthcoming
#' @details Similar to \code{lt_abridged()} details, forthcoming
#' @inheritParams lt_abridged
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
#' @export
lt_single_mx <- function(nMx,
                        Age = 1:length(nMx) - 1,
                        radix = 1e5,
                        a0rule = "ak", 
                        Sex = "m",
                        region = "w",
                        IMR = NA,
                        mod = TRUE,
                        SRB = 1.05,
                        OAG = TRUE,
                        OAnew = max(Age),
                        extrapLaw = NULL,
                        extrapFrom = max(Age),
                        extrapFit = NULL,
                        ...) {

  stopifnot(extrapFrom <= max(Age))
  Sex      <- match.arg(Sex, choices = c("m","f","b"))
  a0rule   <- match.arg(a0rule, choices = c("ak","cd"))
  if (!is.null(extrapLaw)){
    extrapLaw      <- tolower(extrapLaw)
    extrapLaw      <- match.arg(extrapLaw, choices = c("kannisto",
                                                       "kannisto_makeham",
                                                       "makeham",
                                                       "gompertz",
                                                       "ggompertz",
                                                       "beard",
                                                       "beard_makeham",
                                                       "quadratic"
    ))
  } else {
    extrapLaw <- ifelse(max(Age)>=90, "kannisto","makeham")
  }
  
  region   <- match.arg(region, choices = c("w","n","s","e"))
  
  if (is.null(extrapFit)){
    maxAclosed <- ifelse(OAG, Age[which.max(Age)-1],max(Age))
    if (maxAclosed < 85){
      extrapFit  <- Age[Age >= (maxAclosed - 20) & Age <= maxAclosed]
    } else {
      extrapFit  <- Age[Age >= 60 & Age <= maxAclosed]
    }
  } else {
    stopifnot(all(extrapFit %in% Age))
  }
  
  # setup Open Age handling
  OA            <- max(Age)
  # TR: save for later, in case OAG preserved
  if (OAG & OAnew == OA) {
    momega <- nMx[length(nMx)]
  }
  # --------------------------
  # Now all vectors may end up being longer
  x_extr <- seq(extrapFrom, 130, by = 1)
  Mxnew  <- lt_rule_m_extrapolate(
              x = Age,
              mx = nMx,
              x_fit = extrapFit,
              x_extr = x_extr,
              law = extrapLaw,
              ...)

  nMxext        <- Mxnew$values
  Age2          <- names2age(nMxext)

  keepi         <- Age2 < extrapFrom
  nMxext[keepi] <- nMx[Age < extrapFrom]

  # overwrite some variables:
  nMx           <- nMxext
  Age           <- Age2
  N             <- length(Age)
  AgeInt        <- rep(1, N)

  # get ax:
  nAx           <- rep(.5, N)
  if (Age[1] == 0){
    nAx[1]        <- lt_rule_1a0(
                       rule = a0rule,
                       M0 = nMx[1],
                       IMR = IMR,
                       Sex = Sex,
                       region = region,
                       SRB = SRB)
  }
  # get qx (if pathological qx > 1, ax replaced, assumed constant hazard)
  qx            <- lt_id_ma_q(
                     nMx = nMx,
                     nax = nAx,
                     AgeInt = AgeInt,
                     closeout = TRUE)

  lx            <- lt_id_q_l(qx, radix = radix)
  ndx           <- lt_id_l_d(lx)
  nLx           <- lt_id_lda_L(
                     lx = lx,
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
  qx            <- qx[ind]

  lx            <- lx[ind]

  # recalc dx from chopped lx
  ndx           <- lt_id_l_d(lx)
  nLx           <- nLx[ind]
  Tx            <- Tx[ind]
  ex            <- ex[ind]

  # some closeout considerations
  N             <- length(qx)
  qx[N]         <- 1
  nLx[N]        <- Tx[N]
  nAx[N]        <- ex[N]
  AgeInt[N]     <- NA

  # Survival ratios computed only after  nLx is closed out
  Sx            <- lt_id_Ll_S(nLx, lx, AgeInt = AgeInt, N = 1)

  if (OAG) {
    if (OAnew == OA) {
      # keep first estimate if it was open and if we didn't extend
      nMx[N]    <- momega
    } else {
      # Otherwise inner coherence
      nMx[N]    <-  lx[N] / Tx[N]
    }
  } else {
    nMx[N]      <-  lx[N] / Tx[N]
  }

  # output is an unrounded, unsmoothed lifetable
  out <- data.frame(
           Age = Age,
           AgeInt = AgeInt,
           nMx = nMx,
           nAx = nAx,
           nqx = qx,
           lx = lx,
           ndx = ndx,
           nLx = nLx,
           Sx = Sx,
           Tx = Tx,
           ex = ex
  )
  return(out)
}
