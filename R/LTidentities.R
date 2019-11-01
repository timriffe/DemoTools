

# Author: tim
# a bunch of identities that assume abridged ages, but will work for single ages
# if AgeInt is a bunch of 1s
###############################################################################


#' Derive nMx from nqx and nax.
#' @description This is the standard identity to derive nMx from nax and nqx.
#'
#' @param nqx numeric. Vector of age specific death probabilities.
#' @param nax numeric. Vector of average time spent in interval by those dying in interval.
#' @param AgeInt integer. Vector of age class widths.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return nMx vector of age specific death rates derived via identity.
#' @export
qxax2mx <- function(nqx, nax, AgeInt = inferAgeIntAbr(vec = nqx)) {
  nqx / (AgeInt - (AgeInt - nax) * nqx)
}

#' Derive nqx from nMx and nax.
#' @description This is the standard identity to derive nqx from nax and nMx.
#' @details qx values calculated as greater than 1 are imputed with 1.
#' @param nMx numeric. Vector of age-specific death rates.
#' @param nax numeric. Vector of average time spent in interval by those dying in interval.
#' @param AgeInt integer. Vector of age class widths.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return nqx vector of age specific death probabilities derived via identity.
#' @export
mx2qx <- function(nMx, nax, AgeInt = inferAgeIntAbr(vec = nMx)) {
  qx <- (AgeInt * nMx) / (1 + (AgeInt - nax) * nMx)
  ind <- qx > 1 | is.na(qx)
  if (sum(ind) > 0) {
    #cat("at least 1 q(x) > 1, imputed as 1")
    qx[ind] <- 1
  }
  qx
}

#' Derive nax from nqx and nMx.
#' @description This is the standard identity to derive nax from nqx and nMx.
#'
#' @param nqx numeric. Vector of age specific death probabilities.
#' @param nMx numeric. Vector of age-specific death rates.
#' @param AgeInt integer. Vector of age class widths.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return nax numeric vector of average time spent in interval by those dying in interval via identity.
#' @export
qxmx2ax <- function(nqx, nMx, AgeInt) {
  1 / nMx - AgeInt / nqx + AgeInt
}

#' Derive nqx from nMx and nax.
#' @description This is the standard identity to derive nqx from nax and nMx.
#' This is a more full-service wrapper of \code{mx2qx()}, with closeout options and optional age 0
#' treatment.
#' @details qx values calculated as greater than 1 are imputed with 1.
#'
#' @param nMx numeric. Vector of age-specific death rates.
#' @param nax numeric. Vector of average time spent in interval by those dying in interval.
#' @param AgeInt integer. Vector of age class widths.
#' @param closeout logical. Set to 1 if TRUE. Identity otherwise. Default \code{TRUE}.
#' @param IMR numeric. Optional q0 to impute, in case available separately.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return nqx vector of age specific death probabilities derived via identity.
#' @export
mxax2qx <- function(nMx, nax, AgeInt, closeout = TRUE, IMR) {
  qx <- mx2qx(nMx, nax, AgeInt)
  if (closeout) {
    qx[length(qx)] <- 1
    if (length(nMx) == 1) {
      warning(
        "only a single nMx was given, and it was treated as age omega, with qx = 1.\nSpecify closeout = FALSE otherwise"
      )
    }
  }
  if (!missing(IMR) & !is.na(IMR)) {
    qx[1] <- IMR
  }
  ind <- qx > 1
  if (sum(ind) > 0) {
    #cat("at least 1 q(x) > 1, imputed as 1")
    qx[ind] <- 1
  }
  qx
}

#' Derive lifetable survivorship (lx) from death probabilities.
#' @description This lifetable identity is the same no matter what kind of lifetable is required.
#'  You can find it in any demography textbook.
#' @details set \code{radix = 1} for the probability of surviving until age x. The vector returned is
#' the same length as \code{nqx}, thereby throwing out the final value of qx, which is usually set to 1.
#'
#' @param nqx numeric. Vector of age specific death probabilities.
#' @param radix numeric. The lifetable starting population. Default 100000.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return lx vector of lifetable survivorship.
#' @export
qx2lx <- function(nqx, radix = 1e5) {
  radix * cumprod(c(1, 1 - nqx[-length(nqx)]))
}

#' Derive lifetable deaths from survivorship.
#' @description This lifetable identity is the same no matter what kind of lifetable is required.
#'  You can find it in any demography textbook.
#' @details The vector returned is the same length as \code{lx} and it sums to the lifetable radix.
#' If the radix is one then this is the discrete deaths distribution.
#'
#' @param lx numeric.  Vector of age-specific lifetable survivorship.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return ndx vector of lifetable deaths.
#' @export
lx2dx <- function(lx) {
  diff(-c(lx, 0))
}

#' Derive lifetable exposure from lx, ndx and nax.
#' @description This is a common approximation of lifetable exposure:
#' All persons surviving to the end of the interval time the interval width, plus all those that died
#' in the interval multiplied by their average time spent in the interval.
#' @details There is no checking of equal vector lengths at this time.
#' @param lx numeric.  Vector of age-specific lifetable survivorship.
#' @param ndx numeric. Vector of lifetable deaths, summing to radix of \code{lx}.
#' @param nax numeric. Vector of average time spent in interval by those dying in interval.
#' @param AgeInt integer. Vector of age class widths.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return nLx numeric vector of lifetable exposure.
#' @export
lxdxax2Lx <- function(lx, ndx, nax, AgeInt) {
  N                   <- length(lx)
  nLx                 <- rep(0, N)
  nLx[1:(N - 1)]      <-
    AgeInt[1:(N - 1)] * lx[2:N] + nax[1:(N - 1)] * ndx[1:(N - 1)]
  nLx[N]		        <- lx[N] * nax[N]
  nLx
}

#' Derive lifetable total person years left to live from exposure.
#' @description A lifetable identity. Tx is interpreted as the total years
#' left to live above age x in the life table stationary population.
#' @details No \code{NA} or other error checking here. This is taken as the numerator in the classic lifetable method
#'  of calculation of life expectancy.
#' @param Lx numeric. Vector of lifetable exposure.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return Tx total years left to live above age x.
#' @export
Lx2Tx <- function(Lx) {
  rev(cumsum(rev(Lx)))
}

#' adjust ax estimate so that qx less than or equal to 1
#' @description Sometimes a given age interval, death rate, and a(x) imply a death probability that is greater than 1. In this case either the interval needs to be extended or a(x) decreased. This especially arises with mid interval a(x) has been assumed in five-year age groups. This backstop reduces a(x) by assuming a constant death rate over the single ages within the interval, assuming mid interval a(x) for each single age, producing nq(x) by identity from the (5) single ages.
#' @details nMx equal to 2 will imply nqx of 1 by this formula. Implied nqx greater than 1 after this procedure are returned as 1. This is not vectorized!
#' @inheritParams mxax2qx
#' @export
#' @examples
#' # implies a qx > 1
#' mx     <- .65
#' AgeInt <- 5
#' ax     <- 2.5
#' # this is a problematic case
#' (mx2qx(mx, ax, AgeInt))
#' # here a workable value
#' (mxax2qx_Backstop(mx, ax, AgeInt))
#' #' # still less than 1
#' mxax2qx_Backstop(1.99, ax, AgeInt) < 1
#' # unity
#' (mxax2qx_Backstop(2, ax, AgeInt))
#' # qx imputed as 1
#' (mxax2qx_Backstop(2.1, ax, AgeInt))

mxax2qx_Backstop <- function(nMx, nax, AgeInt) {
  stopifnot(length(nMx) == 1)
  # sometime assuming mid interval nAx causes the world
  # to turn upside down.
  qx <- (AgeInt * nMx) / (1 + (AgeInt - nax) * nMx)
  # let's assume flat hazard inside interval
  if (qx > 1) {
    qx <- 1 - prod(rep(1 - mx2qx(nMx, .5, 1), ceiling(AgeInt)))
  }
  if (qx > 1) {
    qx <- 1
  }
  qx
}


#' calculate survivor ratios
#' @description An extra lifetable column for use in projections, which require uniform time steps both both age and period. Intervals are either single age (\code{N=1}) or five-year ages (\code{N=5}). Input vectors are assumed to come from either single or standard abridged ages.
#' @details This function does not account for \code{nLx} having been pre-binned into uniform 5-year age widths, which will throw an error. Just leave them in abridged ages instead. Note that in the case of abridged ages, the interpretation for the first and second value don't follow the original abridged age intervals: the first value in the probability of surviving from birth into ages 0-4 in the first five years, and the second value is the probability of surviving from 0-4 to 5-9. This represents a slight misalignment with the rest of the lifetable, user beware.
#' @inheritParams LTabr
#' @param nLx numeric vector of lifetable exposure.
#' @param N integer, the age width for survivor ratios, either 5 or 1. Default 5.
#' @export
Lxlx2Sx           <- function(nLx, lx, AgeInt, N = c(5, 1)) {
  n               <- length(nLx)
  stopifnot(length(lx) == n)
  # either we're in 1 or 5 year age groups
  stopifnot(length(N) == 1 & N %in% c(5, 1))
  ## compute Sx (missing from the LTbr computation
  Sx              <- rep(NA, n)
  # first age group is survival from births to the second age group
  if (N == 5) {
    # double check because assuming abridged nLx is given...
    stopifnot(length(AgeInt) == n)
    ageintcompare <- inferAgeIntAbr(vec = nLx)
    stopifnot(all(ageintcompare[-n] == AgeInt[-n]))
    # birth until 0-4
    Sx[1]         <- (nLx[1] + nLx[2]) / ((AgeInt[1] + AgeInt[2]) * lx[1])
    # second age group is survival age 0-4 to age 5-9
    Sx[2]         <- nLx[3] / (nLx[1] + nLx[2])
    # middle age groups
    mind          <- 3:(n - 2)
    Sx[mind]      <- nLx[mind + 1] / nLx[mind]
  }
  if (N == 1) {
    LLXX          <- c(lx[1], nLx)
    mind          <- 1:(n - 2)
    Sx[mind]      <- LLXX[mind + 1] / LLXX[mind]
  }
  
  # penultimate age group
  Sx[n - 1]       <- nLx[n] / (nLx[n - 1] + nLx[n])
  # closeout
  Sx[n]           <- 0.0
  
  Sx
}
