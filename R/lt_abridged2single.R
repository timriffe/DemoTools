# A life table by single year of age obtained by graduating the abridged lt using ungroup package

#' create a life table by single year of age by graduating an abridged life table
#' @description Computes single year of age life table by graduating the mortality schedule of
#' an abridged life table, using the `ungroup::pclm()` to ungroup binned count
#' data (`ndx` of the abridged life table) offset by `nLx` life table exposures. Returns fitted single year `nMx`.
#' @details Similar to `lt_single_mx()` details, forthcoming. All steps are constrained to the original `ndx` and `nLx` values.
#' @inheritParams lt_abridged
#' @param ndx numeric. Vector of lifetable deaths at abridged ages. 
#' @param nLx numeric. Vector of lifetable exposure at abridged ages.
#' @param ... optional arguments passed to `pclm()`
#' @details `ndx` and `nLx` could refer to death counts and exposures, respectively.
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
#'
#'
#'  Mx <- c(.23669,.04672,.00982,.00511,.00697,.01036,.01169,
#'          .01332,.01528,.01757,.02092,.02517,.03225,.04241,.06056,
#'          .08574,.11840,.16226,.23745)
#'  Age = c(0,1,seq(5,85,by=5))
#'  AgeInt <- inferAgeIntAbr(vec = Mx)
#'  #
#'  LTabr  <- lt_abridged(nMx = Mx,
#'                       Age = Age,
#'                       AgeInt = AgeInt,
#'                       axmethod = "un",
#'                       Sex = "m",
#'                       mod = TRUE)
#'  ndx <- LTabr$ndx
#'  nLx <- LTabr$nLx
#'  
#'  LT1 <- lt_abridged2single(ndx, 
#'                     nLx, 
#'                     Age, 
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
lt_abridged2single <- function(ndx,
                               nLx,
                               Age,
                               AgeInt = NULL,
                               a0rule = "ak", 
                               Sex = "m",
                               region = "w",
                               IMR = NA,
                               mod = TRUE,
                               SRB = 1.05,
                               OAG = TRUE,
                               OAnew = max(Age),
                               extrapLaw = "kannisto",
                               extrapFrom = max(Age)-1,
                               extrapFit = Age[Age >= 60 & Age < max(Age)],
                               ...) {
  # SH: decimals < 1 don't perform so well in PCLM
  BigRadix  <- 5e6
  radixKeep <- sum(ndx)
  ndx       <- ndx * BigRadix / radixKeep
  nLx       <- nLx * BigRadix / radixKeep
  

  
  
  stopifnot(is_abridged(Age))
  
  NN    <- length(Age)
  stopifnot(length(nLx) == NN & length(ndx) == NN)

  if (is.null(AgeInt)){
    AgeInt <- age2int(Age, OAvalue = 5)
  }
  # graduation will be better if we first extend OAG:
  LTA <- lt_abridged(Deaths = ndx,
                     Exposures = nLx,
                     Age = Age,
                     AgeInt = AgeInt,
                     radix = 5e6,
                     axmethod = "un", 
                     a0rule = a0rule, 
                     Sex = Sex, 
                     IMR = IMR, 
                     region = region, 
                     mod = mod, 
                     SRB = SRB, 
                     OAG = OAG, 
                     OAnew = 130, # we cut down later
                     extrapFit = extrapFit)
  # use pclm to ungroup to single year of age from 1 to maxage_closed+5
  
  ndxE    <- LTA$ndx
  nLxE    <- LTA$nLx
  AgeE    <- LTA$Age
  AgeIntE <- LTA$AgeInt
  AgeIntE[is.na(AgeIntE)] <- 5
  NE      <- length(AgeIntE)
  # TR:
  # adding in a pre-split of nLx, not because it's really needed,
  # but so that we can do an intermediate rescale step to make
  # sure our single-age mx is on the right scale. Hard to describe,
  # but if we *don't* do this life expectancy is greatly exaggerated.
  L1.1 <- pclm(x = AgeE[-NE],
               y = nLxE[-NE],
               nlast = 5)$fitted

  age1    <- 0:length(L1.1 )
  ageint1 <- diff(age1)

  # we can ensure scale of nLx because it's a count
  L1.2 <- rescaleAgeGroups(Value1 = L1.1,
                           AgeInt1 = ageint1,
                           Value2 = nLxE[-NE],
                           AgeInt2 = AgeIntE[-NE],
                           splitfun = graduate_mono)

  # TR: idea is to remove step-artifacts of scaling
  L1.3 <- graduate(Value = L1.2, 
           Age = names2age(L1.2),
           AgeInt = rep(1,length(L1.2)),
           method = "sprague",
           OAG = FALSE,
           keep0 = TRUE)
  
  # exclude 0 and OAG
  use_these <- AgeE > 0 & AgeE < max(AgeE)
  M1.1 <- pclm(x = AgeE[use_these],
               y = ndxE[use_these],
               nlast = 5,
               offset = L1.3[-1])$fitted

  # splice original 1M0 with fitted 1Mx
  M1.1 <- c(ndx[1]/nLx[1], M1.1)

  #
  ndx2 <- M1.1 * L1.2
  # ndx is here rescaled to sum properly, also easier to ensure
  # because it's a count..
  ndx3 <- rescaleAgeGroups(Value1 = ndx2,
                   AgeInt1 = rep(1,length(ndx2)),
                   Value2 = ndxE[-NN],
                   AgeInt2 = AgeIntE[-NN],
                   splitfun = graduate_mono)
  
  ndx4 <-  graduate(Value = ndx3, 
                    Age = names2age(ndx3),
                    AgeInt = age2int( names2age(ndx3),OAG = TRUE, OAvalue = 1),
                    method = "sprague",
                    OAG = FALSE,
                    keep0 = TRUE,
                    constrain = TRUE)
  
  # this version of 1Mx should pass roughly through the middle of the
  # nMx step function implied by the input parameters
  M1.2 <- ndx4 / L1.2
  # Add on OAG
  M1.2 <- c(ndx[1] / nLx[1],M1.2[-1],ndx[NN] / nLx[NN])
  # To see the difference, compare M1.1 with M1.2
  ind <- age1 >= min(extrapFit) & 
    age1 <= (max(extrapFit) + AgeInt[AgeE ==  max(extrapFit)] - 1)
  extrap_fit <- age1[ind]
  # compute life table columns from single year mx
  LT <- lt_single_mx(nMx = M1.2, 
                     Age = age1, 
                     radix = radixKeep,
                     a0rule = a0rule, 
                     Sex = Sex,
                     region = region,
                     IMR = IMR,
                     mod = mod,
                     SRB = SRB,
                     OAG = TRUE,
                     OAnew = OAnew,
                     extrapLaw = extrapLaw,
                     extrapFrom = extrapFrom,
                     extrapFit = extrap_fit) 
  
 return(LT)
  
}

