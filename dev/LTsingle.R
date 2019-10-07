
# lt_single <- function(Deaths, 
#                       Exposures, 
#                       nMx,
#                       nqx,
#                       lx,
#                       Age,
#                       radix = 1e5,
#                       Sex = "m", 
#                       region = "w",
#                       IMR = NA,
#                       mod = TRUE,
#                       OAG = TRUE,
#                       OAnew = max(Age),
#                       extrapLaw = c("Kannisto", "Kannisto_Makeham", "Makeham","Gompertz", "GGompertz", "Beard",
#                                     "Beard_Makeham", "Quadratic")[1],
#                       extrapFrom = max(Age),
#                       extrapFit = Age[Age >= 60],
#                       ...){
#   
#   
# }

#' calculate a single age lifetable
#' @description Fuller description forthcoming
#' @details Similar to \code{LTabr()} details, forthcoming
#' @inheritParams LTabr
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
lt_single_simple <- function(nMx,
                             Age = 1:length(nMx) - 1,
                             radix = 1e5,
                             Sex = "m", 
                             region = "w",
                             IMR = NA,
                             mod = TRUE,
                             OAG = TRUE,
                             OAnew = max(Age),
                             extrapLaw = c("Kannisto", "Kannisto_Makeham", "Makeham","Gompertz", "GGompertz", "Beard",
                                           "Beard_Makeham", "Quadratic")[1],
                             extrapFrom = max(Age),
                             extrapFit = Age[Age >= 60],
                             ...){
  stopifnot(extrapFrom <= max(Age))
  
  Sex           <- tolower(Sex)
  region        <- tolower(region)
  extrapLaw     <- tolower(extrapLaw)
  
  # setup Open Age handling
  OA            <- max(Age)
  # TR: save for later, in case OAG preserved
  if (OAG & OAnew == OA){
    momega <- nMx[length(nMx)]
  }
  # --------------------------
  # Now all vectors may end up being longer
  x_extr <- seq(extrapFrom, 130, by = 1)
  Mxnew  <- extra_mortality(
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
  nAx[1]        <- geta0CD(M0 = nMx[1], IMR = IMR, Sex = Sex, region = region)
  
  # get qx:
  qx            <- mxax2qx_Backstop(nMx = nMx, nax = nAx, AgeInt = AgeInt)
 
  lx            <- qx2lx(nqx, radix = radix)
  ndx           <- lx2dx(lx)
  nLx           <- lxdxax2Lx(lx = lx, ndx = ndx, nax = nAx, AgeInt = AgeInt)
  Tx            <- Lx2Tx(nLx)
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
  ndx           <- lx2dx(lx)
  nLx           <- nLx[ind]
  Tx            <- Tx[ind]
  ex            <- ex[ind]
  
  Sx            <- Lxlx2Sx(nLx, lx, AgeInt = AgeInt, N = 1) 

  # some closeout considerations
  N             <- length(nqx)
  nqx[N]        <- 1
  nLx[N]        <- Tx[N]
  nAx[N]        <- ex[N]
  AgeInt[N]     <- NA
  
  if (OAG){
    if (OAnew == OA){
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
           nqx = nqx,
           lx = lx,
           ndx = ndx,
           nLx = nLx,
           Sx = Sx,
           Tx = Tx,
           ex = ex)
  return(out)
}


