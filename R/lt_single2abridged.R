# An abridged life table that is coherent with an input life table by single year of age

#' calculate an abridged life table that is consistent with a life table by single year of age
#' @description Computes abridged life table columns based on the lx, nLx , and ex values from
#' a single year life table, in accordance with step 2.2 of the Human Life Table Protocol
#' https://www.lifetable.de/methodology.pdf. Output abridged life table has same open age group
#' as input single age life table
#' @details Similar to \code{lt_abridged()} details, forthcoming
#' @param Age integer. Lower bounds of single ages.
#' @param lx numeric. Vector of lifetable survivorship at single ages.
#' @param nLx numeric. Vector of lifetable exposure at single ages.
#' @param ex numeric. Vector of Age-specific remaining life expectancy at single ages.
#' @return Abridged lifetable in data.frame with columns
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
#' 
#' @export
#' 
lt_single2abridged <- function(lx,
                               nLx,
                               ex,
                               Age = 1:length(lx) - 1) {
  
  stopifnot(is_single(Age))
  NN <- length(lx)
  stopifnot(length(nLx) == NN & length(ex) == NN & length(Age) == NN)
  
  # define abridged age groups
  Age5   <- c(0, 1, seq(5, max(Age), 5))
  AgeInt <- age2int(Age = Age5, OAvalue = 5)
  N      <- length(Age5)
  
  # compute abridged lifetable columns
  lx     <- lx[Age %in% Age5]     
  nLx    <- single2abridged(nLx)  
  ex     <- ex[Age %in% Age5]     
  ndx    <- lt_id_l_d(lx)         
  nqx    <- ndx / lx 
  nAx    <- (nLx - (AgeInt * shift.vector(lx,-1,NA))) / ndx
  nAx[N] <- ex[N]
  nMx    <- ndx/nLx
  Tx     <- lt_id_L_T(nLx)
  Sx     <- lt_id_Ll_S(nLx, lx, AgeInt, N = 5)
  
  out <- data.frame(
    Age = Age5,
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

