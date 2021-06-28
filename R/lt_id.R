

# Author: tim
# a bunch of identities that assume abridged ages, but will work for single ages
# if AgeInt is a bunch of 1s
###############################################################################


#' @title Derive nMx from nqx and nax.
#' @description This is the standard identity to derive nMx from nax and nqx.
#'
#' @param nqx numeric. Vector of age specific death probabilities.
#' @param nax numeric. Vector of average time spent in interval by those dying in interval.
#' @param AgeInt integer. Vector of age class widths.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return nMx vector of age specific death rates derived via identity.
#' @export
lt_id_qa_m <- function(nqx, nax, AgeInt) {
  nqx / (AgeInt - (AgeInt - nax) * nqx)
}


#' @title Derive nax from nqx and nMx.
#' @description This is the standard identity to derive nax from nqx and nMx.
#'
#' @param nqx numeric. Vector of age specific death probabilities.
#' @param nMx numeric. Vector of age-specific death rates.
#' @param AgeInt integer. Vector of age class widths.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return nax numeric vector of average time spent in interval by those dying in interval via identity.
#' @export
lt_id_qm_a <- function(nqx, nMx, AgeInt) {
  1 / nMx - AgeInt / nqx + AgeInt
}

# #' Derive nqx from nMx and nax.
# #' @description This is the standard identity to derive nqx from nax and nMx.
# #' This is a more full-service wrapper of \code{lt_id_m_q()}, with closeout options and # optional age 0
# #' treatment.
# #' @details qx values calculated as greater than 1 are imputed with 1.
# #'
# #' @param nMx numeric. Vector of age-specific death rates.
# #' @param nax numeric. Vector of average time spent in interval by those dying in # interval.
# #' @param AgeInt integer. Vector of age class widths.
# #' @param closeout logical. Set to 1 if TRUE. Identity otherwise. Default \code{TRUE}.
# #' @param IMR numeric. Optional q0 to impute, in case available separately.
# #' @references
# #' \insertRef{preston2000demography}{DemoTools}
# #' @return nqx vector of age specific death probabilities derived via identity.
# #' @export
# lt_id_ma_q <- function(nMx, nax, AgeInt, closeout = TRUE, IMR) {
#   if (as.character(match.call()[[1]]) == "mxax2qx") {
#     warning("please use lt_id_ma_q() instead of mxax2qx().", call. = FALSE)
#   }
#
#   qx <- lt_id_m_q(nMx, nax, AgeInt)
#   if (closeout) {
#     qx[length(qx)] <- 1
#     if (length(nMx) == 1) {
#       warning(
#         "only a single nMx was given, and it was treated as age omega, with qx = 1# .\nSpecify closeout = FALSE otherwise"
#       )
#     }
#   }
#   # no Age argument, so strong assumption that if
#   # IMR enters then we're actually starting at age 0!
#   if (!missing(IMR)) {
#     if ( !is.na(IMR)){
#       qx[1] <- IMR
#     }
#   }
#   ind <- qx > 1
#   if (sum(ind) > 0) {
#     #cat("at least 1 q(x) > 1, imputed as 1")
#     qx[ind] <- 1
#   }
#   qx
# }

#' #' @export
#' #' @rdname lt_id_ma_q
#' mxax2qx <- lt_id_ma_q


#' @title Derive lifetable survivorship (lx) from death probabilities.
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
lt_id_q_l <- function(nqx, radix = 1e5) {
  radix * cumprod(c(1, 1 - nqx[-length(nqx)]))
}

#' @title Derive lifetable deaths from survivorship.
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
lt_id_l_d <- function(lx) {
  diff(-c(lx, 0))
}

#' @title Derive lifetable death probabilities from survivorship.
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
lt_id_l_q <- function(lx) {
  dx <- lt_id_l_d(lx)
  dx / lx
}

#' @title Derive survivorship from lifetable deaths
#' @description This lifetable identity is the same no matter what kind of lifetable is required. You can find it in any demography textbook.
#' @details The vector returned is the same length as \code{dx} and it sums to the lifetable radix. If the radix is one then this is the discrete deaths distribution.
#'
#' @param ndx numeric.  Vector of age-specific lifetable deaths.
#' @param radix numeric. 
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return lx vector of lifetable survivorship
#' @export
lt_id_d_l <- function(ndx, radix = sum(ndx)) {
  ndx  <- ndx / sum(ndx)
  N   <- length(ndx)
  CDF <- cumsum(ndx)
  radix * c(1,1 - CDF[-N])
}



#' @title Derive death probabilities from lifetable deaths
#' @description This lifetable identity is the same no matter what kind of lifetable is required.  You can find it in any demography textbook.
#' @details The vector returned is the same length as \code{dx}.
#'
#' @param ndx numeric.  Vector of age-specific lifetable survivorship.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return nqx vector of lifetable death probabilities.
#' @export
lt_id_d_q <- function(ndx) {
  rad <- sum(ndx)
  ndx <- ndx / rad
  N   <- length(ndx)
  CDF <- cumsum(ndx)
  lx  <- c(sum(ndx),1 - CDF[-N])
  ndx / lx
}


#' @title Derive lifetable exposure from lx, ndx and nax.
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
lt_id_lda_L <- function(lx, ndx, nax, AgeInt) {
  N                   <- length(lx)
  nLx                 <- rep(0, N)
  nLx[1:(N - 1)]      <-
    AgeInt[1:(N - 1)] * lx[2:N] + nax[1:(N - 1)] * ndx[1:(N - 1)]
  nLx[N]		        <- lx[N] * nax[N] #open interval
  nLx
}

#' @title Derive lifetable total person years left to live from exposure.
#' @description A lifetable identity. Tx is interpreted as the total years
#' left to live above age x in the life table stationary population.
#' @details No \code{NA} or other error checking here. This is taken as the numerator in the classic lifetable method
#'  of calculation of life expectancy.
#' @param Lx numeric. Vector of lifetable exposure.
#' @references
#' \insertRef{preston2000demography}{DemoTools}
#' @return Tx total years left to live above age x.
#' @export
lt_id_L_T <- function(Lx) {
  rev(cumsum(rev(Lx)))
}

#' @title Calculate conditional death probabilities from nMx and nax
#' @description Sometimes a given age interval, death rate, and a(x) imply a death probability that is greater than 1. In this case either the interval needs to be extended or a(x) decreased. This especially arises with mid interval a(x) has been assumed in five-year age groups. This backstop reduces a(x) by assuming a constant death rate over the single ages within the interval, assuming mid interval a(x) for each single age, producing nq(x) by identity from the (5) single ages.
#' @details nMx equal to 2 will imply nqx of 1 by this formula. Implied nqx greater than 1 after this procedure are returned as 1. This is not vectorized!
#' @inheritParams lt_abridged
#' @param nax numeric. Average time spent in age interval of those that die in the interval.
#' @param closeout logical. Set final age to 1 if TRUE. Identity otherwise. Default \code{TRUE}.
#' @param IMR numeric. Optional q0 to impute, in case available separately.
#' @export
#' @examples
#' # implies a qx > 1
#' mx     <- .65
#' AgeInt <- 5
#' ax     <- 2.5
#' # this is a problematic case
#'  (AgeInt * mx) / (1 + (AgeInt - ax) * mx)
#' # here a workable value
#' (lt_id_ma_q(mx, ax, AgeInt, closeout = FALSE))
#' #' # still less than 1
#' lt_id_ma_q(1.99, ax, AgeInt, closeout = FALSE) < 1
#' # unity
#' (lt_id_ma_q(2, ax, AgeInt, closeout = FALSE))
#' # qx imputed as 1
#' (lt_id_ma_q(2.1, ax, AgeInt, closeout = FALSE))

lt_id_ma_q <- function(nMx, nax, AgeInt, closeout = TRUE, IMR) {
  n <- length(nMx)
  stopifnot(n == length(nax))
  stopifnot(n == length(AgeInt))

  # sometime assuming mid interval nAx causes the world
  # to turn upside down.
  qx <- (AgeInt * nMx) / (1 + (AgeInt - nax) * nMx)

  if (closeout) {
    qx[length(qx)] <- 1
    if (length(nMx) == 1) {
      warning(
        "only a single nMx was given, and it was treated as age omega, with qx = 1.\nSpecify closeout = FALSE otherwise"
      )
    }
  }
  # no Age argument, so strong assumption that if
  # IMR enters then we're actually starting at age 0!
  if (!missing(IMR)) {
    if ( !is.na(IMR)){
      qx[1] <- IMR
    }
  }

  # for cases where qx > 1 let's assume flat hazard inside interval,
  # in essence this overrides the nax given
  ind <- qx > 1
  if (any(ind)){
    ind <- which(ind)
    qx[ind] <-  1 - exp(-nMx[ind] * AgeInt[ind])
    # for (i in ind){
    #   qx[i] <- 1 - ((1 - lt_id_m_q(nMx[i], .5, 1)) ^ ceiling(AgeInt[i]))
    # }
  }
  # final backstop
  qx[qx > 1] <- 1

  qx
}

#' @title Calculate survivor ratios
#' @description An extra lifetable column for use in projections, which require uniform time steps both both age and period. Intervals are either single age or five-year ages. Input vectors are assumed to come from either single or standard abridged ages. Note that the ages of the output Sx are the ages the population would be after the N-year projection.
#' @details This function does not account for \code{nLx} having been pre-binned into uniform 5-year age widths, which will throw an error. Just leave them in abridged ages instead. Note that in the case of abridged ages, the interpretation for the first and second value don't follow the original abridged age intervals: the first value in the probability of surviving from birth into ages 0-4 in the first five years, and the second value is the probability of surviving from 0-4 to 5-9. This represents a slight misalignment with the rest of the lifetable, user beware.
#' @inheritParams lt_abridged
#' @param nLx numeric vector of lifetable exposure.
#' @param Age integer vector of starting ages.
#' @export

lt_id_Ll_S      <- function(nLx, Age) {
  
  # abr or single
  N = ifelse(is_abridged(Age),5,1)
  # infer radix in case was not given
  if(is.null(radix)){
    radix <- ifelse(nLx[1]>1, 10^nchar(trunc(nLx[1])), 1)
  }
  # check validate of nLx
  stopifnot(all(nLx>0 & nLx<radix*N))
  # number ages
  n               <- length(nLx)
  # double check because assuming abridged nLx is given...
  stopifnot(length(Age) == n)
  
  ## compute Sx (missing from the LTbr computation)
  # first age group is survival from births to the second age group
  if (N == 5) {
    Sx              <- rep(NA, n-1)
    # infer always AgeInt
    AgeInt <- inferAgeIntAbr(Age)
    # birth until 0-4
    Sx[1]         <- (nLx[1] + nLx[2]) / ((AgeInt[1] + AgeInt[2]) * radix)
    # second age group is survival age 0-4 to age 5-9
    Sx[2]         <- nLx[3] / (nLx[1] + nLx[2])
    # middle age groups
    mind          <- 3:(n - 2)
    Sx[mind]      <- nLx[mind + 1] / nLx[mind]
    # penultimate age group
    Sx[n - 1]       <- nLx[n] / (nLx[n - 1] + nLx[n])
    # names of ages at arrive
    names(Sx) <- seq(0,Age[length(Age)],5)
  }
  if (N == 1) {
    Sx            <- rep(NA, n)
    LLXX          <- c(radix, nLx)
    mind          <- 1:(n - 1)
    Sx[mind]      <- LLXX[mind + 1] / LLXX[mind]
    # closeout
    Sx[n]           <- nLx[n] / (nLx[n - 1] + nLx[n])
    # names of ages at arrive
    names(Sx) <- 0:(n-1)
  }
  Sx
}

