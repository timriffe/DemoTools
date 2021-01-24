# A life table by single year of age obtained by graduating the abridged lt using ungroup package

#' create a life table by single year of age by graduating an abridged life table
#' @description Computes single year of age life table by graduating the mortality schedule of
#' an abridged life table, using the pclm() function of the ungroup package to ungroup binned count
#' data (ndx of the abridged life table) offset by nLx life table exposures. Returns fitted single year nMx.
#' @details Similar to \code{lt_single_mx()} details, forthcoming
#' @inheritParams lt_abridged
#' @param nMx numeric. Vector of mortality rates at abridged ages.
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
#' 
lt_abridged2single_simple <- function(Age,
                               nMx,
                               radix = 1e5,
                               axmethod = "un",
                               a0rule = "ak", 
                               Sex = "m",
                               region = "w",
                               IMR = NA,
                               mod = TRUE,
                               SRB = 1.05,
                               OAG = TRUE,
                               OAnew = OAnew,
                               extrapLaw = "kannisto",
                               extrapFrom = max(Age),
                               extrapFit = Age[Age >= 60 & Age < max(Age)]) {
  
  stopifnot(is_abridged(Age))
  NN <- length(Age)
  stopifnot(length(nMx) == NN)
  
  # first extend the abridged life table to OAG = 130 with a big radix so that we don't lose info later when rounding ndx and nLx to integers
  lt_abr <- lt_abridged(Age = Age, nMx = nMx, Sex = Sex, radix = 1e8, axmethod = axmethod, a0rule = a0rule,
                        region = region, IMR = IMR, mod = mod, SRB = SRB, OAG = OAG, OAnew = 130,
                        extrapLaw = extrapLaw, extrapFrom = extrapFrom, extrapFit = extrapFit)
  
  # use pclm to ungroup to single year of age from 1 to 129
  # need to round ndx and nLx since pclm doesn't perform with values bw 0 and 1
  M <- pclm(x      = lt_abr$Age[2:27],
            y      = round(lt_abr$ndx[2:27]),
            nlast  = 5,
            offset = round(lt_abr$nLx[2:27]))
  
  # splice original 1M0 with fitted 1Mx and momega from extended abridged LT
  M <- c(nMx[1], fitted(M), lt_abr$ndx[28]/lt_abr$nLx[28])
  
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

