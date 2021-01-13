# A life table by single year of age obtained by graduating the abridged lt using ungroup package

#' create a life table by single year of age by graduating an abridged life table
#' @description Computes single year of age life table by graduating the mortality schedule of
#' an abridged life table, using the pclm() function of the ungroup package to ungroup binned count
#' data (ndx of the abridged life table) offset by nLx life table exposures. Returns fitted single year nMx.
#' @details Similar to \code{lt_single_mx()} details, forthcoming
#' @inheritParams lt_abridged
#' @param ndx numeric. Vector of lifetable deaths at abridged ages.
#' @param nLx numeric. Vector of lifetable exposure at abridged ages.
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
lt_abridged2single <- function(Age,
                               ndx,
                               nLx,
                               radix = 1e5,
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
                               extrapFit = Age[Age >= 60 & Age < max(Age)]) {
  
  require(ungroup) # for pclm function
  
  stopifnot(is_abridged(Age))
  NN <- length(Age)
  stopifnot(length(nLx) == NN & length(ndx) == NN)
  
  # starting age of last closed age group
  maxage_closed <- ifelse(OAG, max(Age) - 5, max(Age))  
  useAge <- Age > 0 & Age <= maxage_closed

  # use pclm to ungroup to single year of age from 1 to maxage_closed+4
  M <- pclm(x      = Age[useAge],
            y      = ndx[useAge],
            nlast  = 5,
            offset = nLx[useAge])
  
  # splice original 1M0 with fitted 1Mx 
  M <- c(ndx[1]/nLx[1], fitted(M))
  
  # compute life table columns from single year mx
  LT <- lt_single_mx(nMx = M, 
                     Age = 1:length(M) - 1, 
                     radix = radix,
                     a0rule = a0rule, 
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
  
 return(LT)
  
}

