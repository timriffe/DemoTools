# A single age lifetable that depends only on qx, first draft, not tested yet.

#' calculate a single age lifetable from qx
#' @description Computes lifetable columns from single age qx by first computing 1ax, then computing
#' 1mx from 1qx and 1ax, and finally passing the 1mx to the lt_single_mx() function
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
#'   \item{Sx}{numeric. Survivor ratios in uniform single-year age groups.}
#'   \item{Tx}{numeric. Lifetable total years left to live above age x.}
#'   \item{ex}{numeric. Age-specific remaining life expectancy.}
#' }
#' @importFrom utils head
#' @export
lt_single_qx <- function(nqx,
                         Age = 1:length(nqx) - 1,
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
  
  # Remove open age group 1qx=1 if it is included in the input vector

  if (OAG == TRUE | nqx[length(nqx)] >= 1.0) {
    Age <- head(Age, -1)
    nqx <- head(nqx, -1)
    extrapFrom <- ifelse(extrapFrom > max(Age), max(Age), extrapFrom)
    extrapFit <- extrapFit[extrapFit <= max(Age)]
  }
  N             <- length(Age)
  AgeInt        <- rep(1, N)
  
  # compute ax:
  nAx           <- rep(.5, N)
  if (Age[1] == 0){
    nAx[1]        <- lt_rule_1a0(rule = a0rule,
                                 q0 = nqx[1],
                                 IMR = IMR,
                                 Sex = Sex,
                                 region = region,
                                 SRB = SRB)
  }
  # compute 1mx from 1qx and 1ax
  nMx          <- lt_id_qa_m(nqx = nqx,
                             nax = nAx,
                             AgeInt = AgeInt)
  
  # pass 1mx to lt_single_mx function
  out           <- lt_single_mx(nMx = nMx,
                                Age = Age,
                                radix = radix,
                                a0rule = "ak", 
                                Sex = Sex,
                                region = region,
                                IMR = IMR,
                                mod = mod,
                                SRB = SRB,
                                OAG = FALSE,
                                OAnew = OAnew,
                                extrapLaw = extrapLaw,
                                extrapFrom = extrapFrom,
                                extrapFit = extrapFit)
  
  return(out)
}
