# TODO this needs to be speed profiled. Why is pclm() slow? Is it just my machine?

# A life table by single year of age obtained by graduating the abridged lt using ungroup package

#' create a life table by single year of age by graduating an abridged life table
#' @description Computes single year of age life table by graduating the mortality schedule of an abridged life table, using the `ungroup::pclm()` to ungroup binned count data. Returns complete single-age lifetable.
#' @details Similar to `lt_abridged()` details, forthcoming. 
#' @inheritParams lt_abridged
#' @param ... optional arguments passed to `pclm()`. For example, if you pass an expicit `lambda` parameter via the `control` argument, you can speed up estimation
#' @return Single-year lifetable in data.frame with columns
#' \itemize{
#'   \item{Age}{integer. Lower bound of single year age class},
#'   \item{AgeInt}{integer. Age class widths.}
#'   \item{nMx}{numeric. Age-specific central death rates.}
#'   \item{nAx}{numeric. Average time spent in interval by those deceased in interval. }
#'   \item{nqx}{numeric. Age-specific conditional death probabilities.}
#'   \item{lx}{numeric. Lifetable survivorship}
#'   \item{ndx}{numeric. Lifetable deaths distribution.}
#'   \item{nLx}{numeric. Lifetable exposure.}
#'   \item{Sx}{numeric. Survivor ratios.}
#'   \item{Tx}{numeric. Lifetable total years left to live above age x.}
#'   \item{ex}{numeric. Age-specific remaining life expectancy.}
#' }
#' 
#' @export
#' @importFrom ungroup pclm
#' @examples
#'  Mx <- c(.23669,.04672,.00982,.00511,.00697,.01036,.01169,
#'          .01332,.01528,.01757,.02092,.02517,.03225,.04241,.06056,
#'          .08574,.11840,.16226,.23745)
#'  Age = c(0,1,seq(5,85,by=5))
#'  AgeInt <- inferAgeIntAbr(vec = Mx)
#'  LTabr <- lt_abridged(nMx = Mx,
#'                       Age = Age, 
#'                       axmethod = "un",
#'                       Sex = "m",
#'                       mod = TRUE)
#'  
#'  LT1 <- lt_abridged2single(nMx = Mx,
#'                     Age = Age, 
#'                     axmethod = "un",
#'                     Sex = "m",
#'                     mod = TRUE)
#' LTabr$ex[1]
#' LT1$ex[1]
#' \dontrun{
#' plot(Age, LTabr$nMx,type = 's', log = 'y')
#' lines(LT1$Age, LT1$nMx)
#' 
#' plot(Age, LTabr$lx,type='S')
#' lines(LT1$Age, LT1$lx)
#' }
lt_abridged2single <- function(
  Deaths = NULL, 
  Exposures = NULL, 
  nMx = NULL, 
  nqx = NULL, 
  lx = NULL,
  Age,
  radix = 1e5,
  axmethod = "un",
  a0rule = "ak", 
  Sex = "m",
  region = "w",
  IMR = NA,
  mod = TRUE,
  SRB = 1.05,
  OAG = TRUE,
  OAnew = max(Age),
  extrapLaw = "kannisto",
  extrapFrom = max(Age),
  extrapFit = Age[Age >= 60 & Age < max(Age)],
  ...) {
  
  stopifnot(is_abridged(Age))
  NN <- length(Age)
  #stopifnot(length(nMx) == NN)
  
  # first extend the abridged life table to OAG = 130 with a big radix so that we don't lose info later when rounding ndx and nLx to integers
  lt_abr <- lt_abridged(Deaths = Deaths, 
                        Exposures = Exposures, 
                        nMx = nMx, 
                        nqx = nqx, 
                        lx = lx, 
                        Age = Age,
                        Sex = Sex, 
                        radix = 1e8, 
                        axmethod = axmethod, 
                        a0rule = a0rule,
                        region = region, 
                        IMR = IMR, 
                        mod = mod, 
                        SRB = SRB, 
                        OAG = OAG, 
                        OAnew = 130,
                        extrapLaw = extrapLaw, 
                        extrapFrom = extrapFrom, 
                        extrapFit = extrapFit)
  
  # use pclm to ungroup to single year of age from 1 to 129
  # need to round ndx and nLx since pclm doesn't perform with values bw 0 and 1
  ind <- lt_abr$Age >= 1 & lt_abr$Age <= 125
  M <- pclm(x      = lt_abr$Age[ind],
            y      = round(lt_abr$ndx[ind]),
            nlast  = 5,
            offset = round(lt_abr$nLx[ind]),
            ...)
  
  # splice original 1M0 with fitted 1Mx and momega from extended abridged LT
  oai <- lt_abr$Age == 130
  M <- c(nMx[1], M$fitted, lt_abr$nMx[oai])
  
  # redefine Age and extrapFit for single year ages
  Age = 1:length(M) - 1
  extrapFit = Age[Age >= min(extrapFit) & Age <= max(Age)] 
  
  # compute life table columns from single year mx
  LT <- lt_single_mx(nMx = M, 
                     Age = Age, 
                     radix = radix,
                     a0rule = a0rule, 
                     Sex = Sex,
                     region = region,
                     IMR = IMR,
                     mod = mod,
                     SRB = SRB,
                     OAG = TRUE,
                     OAnew = OAnew,
                     extrapLaw = extrapLaw,
                     extrapFrom = min(extrapFrom, 110), # always refit from 110 even if extrapFrom > 110
                     extrapFit = extrapFit) 
  
  return(LT)
  
}
