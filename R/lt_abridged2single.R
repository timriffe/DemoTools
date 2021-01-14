# A life table by single year of age obtained by graduating the abridged lt using ungroup package

#' create a life table by single year of age by graduating an abridged life table
#' @description Computes single year of age life table by graduating the mortality schedule of
#' an abridged life table, using the `ungroup::pclm()` to ungroup binned count
#' data (`ndx` of the abridged life table) offset by `nLx` life table exposures. Returns fitted single year `nMx`.
#' @details Similar to `lt_single_mx()`` details, forthcoming
#' @inheritParams lt_abridged
#' @param ndx numeric. Vector of lifetable deaths at abridged ages.
#' @param nLx numeric. Vector of lifetable exposure at abridged ages.
#' @param ... optional arguments passed to `pclm()`
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
#' @importFrom stats fitted
#' @examples
#'  Mx <- c(.23669,.04672,.00982,.00511,.00697,.01036,.01169,
#'          .01332,.01528,.01757,.02092,.02517,.03225,.04241,.06056,
#'          .08574,.11840,.16226,.23745)
#'  excheckUN <-  c(35.000,42.901,47.190,44.438,
#'                  40.523,36.868,33.691,30.567,27.500,24.485,21.504,18.599,
#'                  15.758,13.080,10.584,8.466,6.729,5.312,4.211)
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
  radix <- sum(ndx)
  stopifnot(is_abridged(Age))
  
  NN    <- length(Age)
  stopifnot(length(nLx) == NN & length(ndx) == NN)

  if (is.null(AgeInt)){
    AgeInt <- age2int(Age, OAvalue = 5)
  }

  # use pclm to ungroup to single year of age from 1 to maxage_closed+5

  # TR:
  # adding in a pre-split of nLx, not because it's really needed,
  # but so that we can do an intermediate rescale step to make
  # sure our single-age mx is on the right scale. Hard to describe,
  # but if we *don't* do this life expectancy is greatly exaggerated.
  L1.1 <- pclm(x = Age[-length(Age)],
               y = nLx[-length(Age)],
               nlast = 5)

  age1 <- 0:length(L1.1$fitted)
  ageint1 <- diff(age1)

  # we can ensure scale of nLx because it's a count
  L1.2 <- rescaleAgeGroups(Value1 = L1.1$fitted,
                           AgeInt1 = ageint1,
                           Value2 = nLx[-NN],
                           AgeInt2 = AgeInt[-NN],
                           splitfun = graduate_mono)

  # exclude 0 and OAG
  use_these <- Age > 0 & Age < max(Age)
  M1.1 <- pclm(x = Age[use_these],
               y = ndx[use_these],
               nlast = 5,
               offset = L1.2[-1])

  # splice original 1M0 with fitted 1Mx
  M1.1 <- c(ndx[1]/nLx[1], M1.1$fitted)

  #
  ndx2 <- M1.1 * L1.2
  # ndx is here rescaled to sum properly, also easier to ensure
  # because it's a count..
  ndx3 <- rescaleAgeGroups(Value1 = ndx2,
                   AgeInt1 = ageint1,
                   Value2 = ndx[-NN],
                   AgeInt2 = AgeInt[-NN],
                   splitfun = graduate_mono)
  
  # this version of 1Mx should pass roughly through the middle of the
  # nMx step function implied by the input parameters
  M1.2 <- ndx3 / L1.2
  # Add on OAG
  M1.2 <- c(M1.2,ndx[NN] / nLx[NN])
  # To see the difference, compare M1.1 with M1.2
  ind <- age1 >= min(extrapFit) & 
    age1 <= (max(extrapFit) + AgeInt[Age ==  max(extrapFit)] - 1)
  extrap_fit <- age1[ind]
  # compute life table columns from single year mx
  LT <- lt_single_mx(nMx = M1.2, 
                     Age = age1, 
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
                     extrapFrom = extrapFrom,
                     extrapFit = extrap_fit) 
  
 return(LT)
  
}

