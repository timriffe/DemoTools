

# Author: tim
###############################################################################
# lots of approaches to a(0). Sometimes called 'separation factors'. These ought to be
# collected, organized, and standardized here. We have two major versions here:
# PAS (mostly uniform) and UN (mostly greville-based).

#' Coale-Demeny a(0) as function of m(0), region, and sex.
#'
#' @description Coale-Demeny a(0) from Manual X Table 164. This is a rule of thumb.
#' In this and some other older texts, a(0) is known as a 'separation factor'.
#'
#' @param M0 numeric. Event exposure infant mortality rate.
#' @param IMR numeric. Optional. {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}}, the death probability in first year of life, in case available separately.
#' @param Sex character. \code{"m"} or \code{"f"}  for male or female,
#' @param region character. \code{"n"}, \code{"e"}, \code{"s"} or \code{"w"} for North, East, South, or West.
#'
#' @details If \code{IMR} is not given, then \code{M0} is converted to q(0) using the following approximation:
#' \enumerate{
#' \item{Find \eqn{\alpha , \beta}.}{ Look up the appropriate slope and intercept for the given sex and region.}
#' \item{calculate \eqn{a} as: }{\ifelse{html}{\out{a = M<sub>0</sub> * &beta;}}{\eqn{a = M_0 * \beta}}}
#' \item{calculate \eqn{b} as: }{\ifelse{html}{\out{b = 1 + M<sub>0</sub> *(1- &alpha;)}}{\eqn{b =  1 + M_0 * (1 - \alpha)}}}
#' \item{approximate {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}} as:}{ \ifelse{html}{\out{q<sub>0</sub> = (b<sup>2</sup>- &radic; [b -4*a*M<sub>0</sub>]) / (2*a)}}{\eqn{q_0 = \frac{ b - sqrt(b^2 - 4 * a * M_0) }{ 2 * a } }}}
#' \item{use {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}} as}{ IMR, and applied directly to the Coale-Demeny piecewise linear formula.}
#' }
#' If \code{IMR} is given, then \code{M0} is disregarded, and transitivity is therefore not guaranteed. In this case, one has the option to use \code{lt_id_qm_a()} to derive \code{a(0)}, however discrepancies between these two parameters could force implausible results in \code{a(0)}, whereas the CD rule always gives something plausible.
#'
#' @references
#' \insertRef{united1983manual}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#'
#' @return The average age at death in the first year of life a(0).
#' @export
#' @examples
#' m0 <- seq(.001, .2, by = .001)
#' \dontrun{
#' plot(m0, sapply(m0, lt_rule_1a0_cd, Sex = "m", region = "e"), ylab = "a0",
#' 		type = 'l', ylim = c(0,.36), lty = 2, col = "blue")
#' lines(m0,sapply(m0, lt_rule_1a0_cd, Sex = "m", region = "w"), col = "blue")
#' lines(m0,sapply(m0, lt_rule_1a0_cd, Sex = "f", region = "e"), lty = 2, col = "red")
#' lines(m0,sapply(m0, lt_rule_1a0_cd, Sex = "f", region = "w"), col = "red")
#' text(.15, lt_rule_1a0_cd(.15,Sex = "m", region = "e"),"males E",font=2)
#' text(.15, lt_rule_1a0_cd(.15,Sex = "m", region = "w"),"males N,W,S",font=2)
#' text(.15, lt_rule_1a0_cd(.15,Sex = "f", region = "e"),"females E",font=2)
#' text(.15, lt_rule_1a0_cd(.15,Sex = "f", region = "w"),"females N,W,S",font=2)
#'
#' # compare with the Preston approximation
#' # constants identical after m0 = .107
#' m0 <- seq(.001,.107,by =.001)
#' a0CDm0 <- sapply(m0, lt_rule_1a0_cd, Sex = "m", region = "w")
#' a0CDpr <- 0.045 + 2.684 * m0
#' plot(m0, a0CDm0, type = 'l', lty = 2, col = "red")
#' lines(m0, a0CDpr)
#' plot(m0, (a0CDm0 - a0CDpr) * 365, main = "difference (days)", ylab = "days")
#'}
# this is called a separation factor in the spreadsheet?
# separate estimate of IMR optional
lt_rule_1a0_cd <- function(M0,
                    IMR = NA,
                    Sex = "m",
                    region = "w") {
 
  # sex can be "m", "f", or "b"
  # region can be "n","e","s","w",or
  Sex    <- match.arg(Sex, choices = c("m","f"))
  region <- match.arg(region, choices = c("w","n","e","s"))


  Age0Const <- matrix(c( 0.33, 0.35,0.33, 0.35, 0.29, 0.31, 0.33, 0.35),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(c("w", "n", "e", "s"), c("m", "f"))
  )

  Intercept <- matrix(
    c( 0.0425, 0.05, 0.0425, 0.05, 0.0025, 0.01, 0.0425, 0.05),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(c("w", "n", "e", "s"), c("m", "f"))
  )
  Slope        <- c(2.875, 	3.000)
  names(Slope) <- c("m", "f")

  Alpha        <- Intercept[region, Sex]
  Beta         <- Slope[Sex]
  # IMR optional here, use approximation
  if (missing(IMR) | is.na(IMR)) {
    # formula from PAS LTPOPDTH
    a          <- M0 * Beta
    b          <- 1 + M0 * (1 - Alpha)
    SQRTmiddle <- b ^ 2 - 4 * a * M0
    if (SQRTmiddle <= 0) {
      IMR      <- .2 # just to trigger constant...
    } else {
      IMR      <- (b - sqrt(SQRTmiddle)) / (2 * a)
    }
  }
  ifelse(
    IMR > .1,
    Age0Const[region, Sex],
    {Alpha + Beta * IMR})
}

# Separate estimate of IMR optional
# TR: I think it's funny that a1-4 doesn't depend at all on m1-4

#' Coale-Demeny 4a1 as function of M(0), region, and sex.
#'
#' @description Coale-Demeny 4a1. This is a rule of thumb. In this and some other older texts, 4a1 is known as a 'separation factor'. These coefficients were pulled from the PAS spreadsheets \code{LTPOPDTH.XLS} and not located in the original Manual X.
#'
#' @details If sex is given as both, \code{"b"}, then we calculate the male and female results separately, then weight them together using SRB. This is bad in theory, but the leverage is trivial, and it's better than using male or female coefs for the total population. If \code{IMR} is not given, then \code{M0} is used in its stead.
#'
#' @param M0 numeric. Event exposure infant mortality rate.
#' @param a0rule character. Either \code{"ak"} (default) or \code{"cd"}.
#' @param IMR numeric. Optional. {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}}, the death probability in first year of life, in case available separately.
#' @param Sex character. \code{"m"}, \code{"f"} or \code{"b"} for male, female, or both.
#' @param region character. \code{"n"}, \code{"e"}, \code{"s"} or \code{"w"} for North, East, South, or West.
#' @param SRB sex ratio at birth.
#'
#' @return The average age at death between ages 1-4, 4a1.
#' @export
#' @references
#' \insertRef{united1983manual}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#' @examples
#' m0 <- seq(.001,.2,by =.001)
#' \dontrun{
#' 
#' # using Andreev-Kingkade for a0 (it makes no difference if you use cd actually...)
#' plot(m0, sapply(m0, lt_rule_4a1_cd, Sex = "m", region = "e"), ylab = "4a1",
#' 		type = 'l', ylim = c(1,2), lty = 2, col = "blue")
#' lines(m0,sapply(m0, lt_rule_4a1_cd, Sex = "m", region = "w"), col = "blue")
#' lines(m0,sapply(m0, lt_rule_4a1_cd, Sex = "m", region = "n"), col = "blue", lty = "8383",lwd=2)
#' lines(m0,sapply(m0, lt_rule_4a1_cd, Sex = "m", region = "s"), col = "blue", lty = "6464",lwd=2)
#' lines(m0,sapply(m0, lt_rule_4a1_cd, Sex = "f", region = "e"), lty = 2, col = "red")
#' lines(m0,sapply(m0, lt_rule_4a1_cd, Sex = "f", region = "w"), col = "red")
#' lines(m0,sapply(m0, lt_rule_4a1_cd, Sex = "f", region = "n"), col = "red", lty = "8383",lwd=2)
#' lines(m0,sapply(m0, lt_rule_4a1_cd, Sex = "f", region = "s"), col = "red", lty = "6464",lwd=2)
#' 
#' text(.05, lt_rule_4a1_cd(.05,Sex = "m", region = "e"),"males E",font=2,pos=4)
#' text(.05, lt_rule_4a1_cd(.05,Sex = "m", region = "w"),"males W",font=2,pos=4)
#' text(.05, lt_rule_4a1_cd(.05,Sex = "m", region = "s"),"males S",font=2,pos=4)
#' text(.05, lt_rule_4a1_cd(.05,Sex = "m", region = "n"),"males N",font=2,pos=4)
#' 
#' text(0, lt_rule_4a1_cd(.01,Sex = "f", region = "e"),"females E",font=2,pos=4)
#' text(0, lt_rule_4a1_cd(.01,Sex = "f", region = "w"),"females W",font=2,pos=4)
#' text(0, lt_rule_4a1_cd(.01,Sex = "f", region = "s"),"females S",font=2,pos=4)
#' text(0, lt_rule_4a1_cd(.01,Sex = "f", region = "n"),"females N",font=2,pos=4)
#'
#' }
lt_rule_4a1_cd <- function(M0,
                           a0rule = "ak",
                           IMR = NA,
                           Sex = "m",
                           region = "w",
                           SRB = 1.05) {

  a0rule   <- match.arg(a0rule, choices =  c("ak","cd"))
  Sex      <- match.arg(Sex, choices = c("m","f","b"))
  region   <- match.arg(region, choices =c("w","n","s","e"))
  
  #
  if (Sex == "b"){
    if (missing(M0)){
      M0 <- NA
    }
    a1f <- lt_rule_4a1_cd(M0 = M0, 
                          a0rule = a0rule,
                          IMR = IMR,
                          Sex = "f",
                          region =region,
                          SRB = SRB)
    a1m <- lt_rule_4a1_cd(M0 = M0, 
                          a0rule = a0rule,
                          IMR = IMR,
                          Sex = "m",
                          region =region,
                          SRB = SRB)
    pm <- SRB / (SRB + 1)
    a1 <- pm * a1m + (1 - pm) * a1f
    return(a1)
  }
  

  Age1_4Const <- matrix(
    c( 1.352, 1.361, 1.558,1.570,1.313, 1.324,1.240, 1.239),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(c("w", "n", "e", "s"), c("m", "f"))
  )

  Intercept <- matrix(
    c( 1.653, 1.524,  1.859, 1.733,  1.614, 1.487, 1.541, 1.402),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(c("w", "n", "e", "s"), c("m", "f"))
  )
  Slope        <- c(3.013, 	1.627)
  names(Slope) <- c("m", "f")
  if (missing(IMR) | is.na(IMR)) {
    a0 <- lt_rule_1a0(
                  rule = a0rule,
                  M0 = M0,
                  IMR = NA,
                  Sex = Sex,
                  region = region,
                  SRB = SRB)
    IMR <-
      lt_id_ma_q(
        nMx = M0,
        nax = a0,
        AgeInt = 1,
        closeout = FALSE,
        IMR = NA
      )
  }
  ifelse(IMR > .1, Age1_4Const[region, Sex], Intercept[region, Sex] - Slope[Sex] * IMR)
}

#' PAS a(x) rule of thumb.
#'
#' @description a(x) is calculated following the Coale-Demeny rules for ages 0 and 1-4, and assumes interval midpoints in higher ages.
#' This is just a rule of thumb. This procedure is as found in the PAS spreadsheet \code{LTPOPDTH.XLS}.
#'
#' @details If sex is given as both, \code{"b"}, then female values are taken for a(0) and 4a1, per the PAS spreadsheet. If IMR is not given, the M(0) is used to estimate a(x) for ages < 5. This function is not vectorized. a(x) closeout assumes constant mortality hazard in the open age group. One safeguard is different from PAS: If assuming the interval midpoint implies a qx greater than 1, then we derive a(x) for the interval by assuming midpoint a(x) for each single age within the interval along with a constant death rate.
#'
#' @param nMx numeric. Event exposure mortality rates.
#' @param AgeInt integer. Vector of age interval widths.
#' @param a0rule character. Either \code{"ak"} (default) or \code{"cd"}.
#' @param IMR numeric. Optional. {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}}, the death probability in first year of life, in case available separately.
#' @param Sex character. \code{"m"}, \code{"f"} or \code{"b"} for male, female, or both.
#' @param region character. \code{"n"}, \code{"e"}, \code{"s"} or \code{"w"} for North, East, South, or West.
#' @param SRB numeric. The sex ratio at birth (boys/girls), default 1.05.
#' @param OAG logical. Whether or not the last element of \code{nMx} is the open age group Default \code{TRUE}.
#'
#' @return nax average contribution to exposure of those dying in the interval.
#' @export
#' @references
#' \insertRef{united1983manual}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#' @examples
#'  Exposures <- c(100958,466275,624134,559559,446736,370653,301862,249409,
#'  		247473,223014,172260,149338,127242,105715,79614,53660,
#'  		31021,16805,8000,4000,2000,1000)
#' 
#'  Deaths <- c(8674,1592,618,411,755,1098,1100,1357,
#'  		1335,3257,2200,4023,2167,4578,2956,4212,
#'  		2887,2351,1500,900,500,300)
#'  # lower age bounds
#'  Age    <- c(0, 1, seq(5, 100, by = 5))
#'  AgeInt <- c(diff(Age), NA)
#'  nMx <- Deaths/Exposures
#'  lt_a_pas(nMx = nMx,AgeInt = AgeInt,Sex = 'm',region = 'n',OAG = TRUE)
#'  lt_a_pas(nMx = nMx,AgeInt = AgeInt,Sex = 'm',a0rule = "cd",region = 'n',OAG = TRUE)
lt_a_pas <-
  function(nMx,
           AgeInt,
           a0rule = "ak",
           IMR = NA,
           Sex = "m",
           region = "w",
           SRB = 1.05,
           OAG = TRUE) {

    a0rule   <- match.arg(a0rule, choices =  c("ak","cd"))
    Sex      <- match.arg(Sex, choices = c("m","f","b"))
    region   <- match.arg(region, choices =c("w","n","s","e"))
    
    N      <- length(nMx)
    ax     <- AgeInt / 2

    # TR: 31-08-2020 updated this. Now AK by default
    ax[1]  <- lt_rule_1a0(
                rule = a0rule,
                M0 = nMx[1],
                IMR = IMR,
                Sex = Sex,
                region = region,
                SRB = SRB)
    
    # TR: would be great to update this tiny method...
    ax[2]  <- lt_rule_4a1_cd(
                M0 = nMx[1],
                IMR = IMR,
                Sex = Sex,
                region = region,
                SRB = SRB)

    # TR 5-12-2019 this two-step returns midpoints IFF they
    # don't imply qx > 1, otherwise, for such pathological
    # values it returns the ax implied by a constant hazard.
    # TR 4-1-2020 removed. caes of qx > 1 can be dealt with
    # downstream if necessary.
    # qx <- lt_id_ma_q(nMx = nMx,
    #                  nax = ax,
    #                  AgeInt = AgeInt,
    #                  closeout = OAG,
    #                  IMR = IMR)
    # # TR: this step vulnerable, see comment below.
    # ax <- lt_id_qm_a(nqx = qx, nMx = nMx, AgeInt = AgeInt)
    ax[N]  <- ifelse(OAG, 1 / nMx[N], ax[N])

    # TR 31-12-2019 at this point, transitivity is guaranteed internally for
    # qx, mx, and ax, however, a0 may have been determied by IMR, and in this
    # case q0 (just created above) might be different from IMR (likely so), AND
    # the a0 might be negative or greater than 1, which simply cannot be. It
    # should be preferred to EITHER override M0 downstream (beyond this function call)
    # OR remove the IMR option altogether.

    # ind <- impliedqx > 1
    # if (sum(ind) > 0) {
    #   for (i in which(ind)) {
    #     qxnew <- lt_id_ma_q_robust(nMx = nMx[i],
    #                        nax = ax[i],
    #                        AgeInt = AgeInt[i])
    #     ax[i] <- lt_id_qm_a(nqx = qxnew,
    #                      nMx = nMx[i],
    #                      AgeInt = AgeInt[i])
    #   }
    # }

    ax
  }

#' UN version of the Greville formula for a(x) from M(x)
#'
#' @description The UN a(x) formula uses Coale-Demeny for ages 0, and 1-4, values of 2.5 for ages 5-9 and 10-14, and the Greville formula thereafter. In the original sources these are referred to as separation factors.
#'
#' @param nMx numeric. Event exposure mortality rates.
#' @param nqx numeric.  Vector of age specific death probabilities in standard abridged age groups.
#' @param lx numeric.  Vector of lifetable survivorship in standard abridged age groups.
#' @param Age integer. Vector of lower bounds of abridged age groups.
#' @param AgeInt integer. Vector of age group intervals.
#' @param a0rule character. Either `"ak"` (default) or `"cd"`.
#' @param IMR numeric. Optional. {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}}, the death probability in first year of life, in case available separately.
#' @param Sex character. `"m"`, `"f"` or `"b"` for male, female, or both.
#' @param region character. `"n",` `"e"`, `"s"` or `"w"` for North, East, South, or West.
#' @param mod logical. Whether or not to use Gerland's modification for ages 5-14. Default \code{TRUE}.
#' @param SRB numeric. The sex ratio at birth (boys/girls), default 1.05.
#' @param closeout logical. Whether or not to estimate open age a(x) via extrapolation. Default \code{TRUE}.
#' @inheritParams lt_a_closeout
#'
#' @details \eqn{a(x)} for age 0 and age group 1-4 are based on Coale-Demeny {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}}-based lookup tables. An approximation to get from M(0) to {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}} for the sake of generating a(0) and 4a1 is used. The final a(x) value is closed out using the `lt_a_closeout()` method (reciprocal and Mortpak methods are deprecated). Age groups must be standard abridged. No check on age groups is done.
#'
#' There are different vectors one can specify for this method: ultimately it's either `nMx` or `nqx`, and the `nax`, but `nax` results should be reasonably close. For full convergence implying transitivity, instead use `lt_a_un()`

#' @seealso
#' \code{\link[DemoTools]{lt_a_closeout}}
#' @references
#' \insertRef{greville1977short}{DemoTools}
#' \insertRef{un1982model}{DemoTools}
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertRef{mortpak1988}{DemoTools}
#' @return nax average contribution to exposure of those dying in the interval.
#' @export

#' @examples
#' #Example witn Mexican data from UN
#' nMx <- c(0.11621,0.02268,0.00409,0.00212,0.00295,0.00418,0.00509,0.00609,
#' 0.00714,0.00808,0.00971,0.0125,0.0175,0.02551,0.03809,0.05595,0.08098,
#' 0.15353,0.2557)
#'
#' nqx <- c(0.10793,0.08554,0.02025,0.01053,0.01463,0.02071,0.02515,0.02999,
#' 0.03507, 0.03958,0.04742,0.0606,0.08381,0.11992,0.17391,0.2454,0.33672,
#' 0.54723,NA)
#'
#' lx  <- c(100000,89207,81577,79924,79083,77925,76312,74393,72162,69631,66875,
#' 63704,59843,54828,48253,39861,30079,19951,9033)
#'
#' Age <- c(0,1,seq(5,85,by = 5))
#' AgeInt <- age2int(Age, OAvalue = 5)
#' # two quite different results depending whether you start with mx or qx
#' lt_id_morq_a_greville(nMx = nMx, 
#'                       Age = Age, 
#'                       AgeInt = AgeInt, 
#'                       Sex = 'f',
#'                       region = 'w')
#' lt_id_morq_a_greville(nqx = nqx, 
#'                       Age = Age, 
#'                       AgeInt = AgeInt, 
#'                       Sex = 'f',
#'                       region = 'w')
#' # same, qx comes from lx (Except OAG, because qx closed above w NA)
#' lt_id_morq_a_greville(lx = lx, 
#'                       Age = Age, 
#'                       AgeInt = AgeInt, 
#'                       Sex = 'f',
#'                       region = 'w')
#' # both qx and lx given, but lx not used for anything = same as qx
#' lt_id_morq_a_greville(nqx = nqx, 
#'                       lx = lx, 
#'                       Age = Age, 
#'                       AgeInt = AgeInt, 
#'                       Sex = 'f',
#'                       region = 'w')
#'#'
#' # nMx taken over lx.
#' lt_id_morq_a_greville(nMx = nMx,  
#'                       lx = lx, 
#'                       Sex = 'f', 
#'                       Age = Age, 
#'                       AgeInt = AgeInt, 
#'                       region = 'w')
#' 
#' # example of transitivity. UN ax method is iterative
#' # and used greville on the inside
#' 
#' nMx <- c(0.11621,0.02268,0.00409,0.00212,0.00295,
#'          0.00418,0.00509,0.00609,0.00714,0.00808,
#'          0.00971,0.0125,0.0175,0.02551,0.03809,
#'          0.05595,0.0809,0.15353,0.2557)
#' AgeInt <- inferAgeIntAbr(vec = nMx)
#' Age    <- int2age(AgeInt)
#' nAx1   <- lt_a_un(nMx, 
#'                   Age = Age,
#'                   AgeInt = AgeInt, 
#'                   a0rule = "ak",
#'                   mod = TRUE, 
#'                   closeout = TRUE)
#' 
#' nqx    <- lt_id_ma_q(nMx = nMx, 
#'                      nax = nAx1, 
#'                      AgeInt = AgeInt, 
#'                      closeout = FALSE)
#' nAx2 <- lt_a_un(nqx = nqx, 
#'                 Age = Age,
#'                 AgeInt = AgeInt, 
#'                 a0rule = "ak", 
#'                 mod= TRUE, 
#'                 closeout = TRUE)
#' 
#' stopifnot(all(abs(nAx1[Age<75] - nAx2[Age<75]) < 1e-7))
#'

lt_id_morq_a_greville <- function(nMx,
                                  nqx,
                                  lx,
                                  Age,
                                  AgeInt = age2int(Age, OAvalue = 5),
                                  a0rule = "ak", 
                                  IMR = NA,
                                  Sex = "m",
                                  region = "w",
                                  mod = TRUE,
                                  SRB = 1.05,
                                  closeout = TRUE,
                                  extrapLaw = "kannisto",
                                  extrapFrom = max(Age),
                                  extrapFit = Age[Age >= 60],
                                  ...) {

  Sex      <- match.arg(Sex, choices = c("m","f","b"))
  a0rule   <- match.arg(a0rule, choices =  c("ak","cd"))
  extrapLaw      <- match.arg(extrapLaw, choices = c("kannisto",
                                         "kannisto_makeham",
                                         "makeham",
                                         "gompertz",
                                         "ggompertz",
                                         "beard",
                                         "beard_makeham",
                                         "quadratic"
  ))
  region   <- match.arg(region, choices =c("w","n","s","e"))
  
  DBL_MIN <- .Machine$double.xmin

  # sort out arguments:
  mxflag  <- missing(nMx)
  qxflag  <- missing(nqx)
  imr_flag <- !is.na(IMR)



  # if no qx, we can get from lx if available
  if (qxflag & !missing(lx)) {
    nqx <- lt_id_l_d(lx) / lx
    qxflag <- FALSE
  }
  if (!qxflag & imr_flag){
    # preserve IMR like so
    nqx[1] <- IMR
  }

  stopifnot(!qxflag | !mxflag)
  # now we have either qx or mx

  if (!mxflag) {
    a0     <- lt_rule_1a0(
                rule = a0rule,
                M0 = nMx[1],
                IMR = IMR,
                Sex = Sex,
                region = region,
                SRB = SRB)
    a1_4   <- lt_rule_4a1_cd(
                M0 = nMx[1],
                IMR = IMR,
                Sex = Sex,
                region = region,
                SRB = SRB)
    # qind slightly different from qxflag?
    qind   <- FALSE
  }

  if (mxflag & !qxflag) {
    # TR: in any case nqx[1] is IMR, if IMR was given.
    a0     <- lt_rule_1a0(
                rule = a0rule,
                q0 = nqx[1],
                IMR = nqx[1],
                Sex = Sex,
                region = region,
                SRB = SRB)
    a1_4   <- lt_rule_4a1_cd(
                M0 = NA,
                IMR = nqx[1],
                Sex = Sex,
                region = region,
                SRB = SRB)
    # here nMx created, but mxflag upheld
    nMx    <- nqx
    qind   <- TRUE
  }


  # some setup
  N           <- length(nMx)

  ## for ages 15-19 onward
  ## AK <- log(QxMx[j+1]/QxMx[j-1])/10 ## we don't have 1/10 ?
  ## ax[j] <- 2.5 - (25.0/12.0) * (QxMx[j] - AK)

  ## improved Greville formula for adolescent ages 5-9 and 10-14
  ## Let the three successive lengths be n1, n2 and n3, the formula for 5a5 is:
  ## ax[i] = 2.5 - (25 / 12) * (mx[i] - log(mx[i + 1] / mx[i-1])/(n1/2+n2+n3/2))
  ## for age 5-9, coefficient should be 1/9.5, because age group 1-4
  ## has only 4 ages (not 5), while the other 5-year age group are 1/10
  ## ax[i] = 2.5 - (25 / 12) * (mx[i] - (1/9.5)* log(mx[i + 1] / mx[i-1]))
  ## Age 20-25, ..., 95-99

  # TR: 6-3-2020 remove pointless loop, add in new back stop

  # constant hazard ax, as filler for pathological cases:
  axConst <- AgeInt + 1 / nMx - AgeInt / (1 - exp(-AgeInt * nMx))
  # inverse of weighted moving avg age interval.
  Abx     <- 1 / (shift.vector(AgeInt,1, fill = NA) / 2 + AgeInt + shift.vector(AgeInt,-1, fill = NA) / 2)
  # Kx is fragile
  Kx      <- log(pmax(shift.vector(nMx,-1, fill = NA),DBL_MIN) /
                   pmax(shift.vector(nMx,1, fill = NA),DBL_MIN)) # / 10?
  # This does nothing special for age 5-9, which has a 4-year age group on the left.
  # age 1-4 is overwritten anyway.
  ax      <- AgeInt / 2 - (AgeInt ^ 2 / 12) * (nMx - Kx * Abx)

  # pick out likely pathological cases and imput w const ax
  ind     <- Age >= 10 & Age < max(Age) & (ax < (.4 * AgeInt) | ax > (.6 * AgeInt)) | is.na(ax)
  ax[ind] <- axConst[ind]

  # TR: remove this, superceded by axConst filler.
  ## (Mortpak used 1 instead of 0.97).
  # if (Age[i] > 35 && ax[i] < 0.97) {
  #   ax[i] <- 0.97
  # }


  # preserve old patch for
  if (qind){
    axfill   <- shift.vector(ax,1, fill = NA) * .8
    ind      <- Age >= 65 & ax < axfill
    ax[ind]  <- axfill[ind]
  }

  # a1_4 only valid for abridged
  ax[1:2] <- c(a0, a1_4)

  # TR: this is questionable now
  if (!mod) {
    ax[3:4] <- AgeInt[3:4] / 2
  }

  # closeout
  if (max(Age) < 130 & closeout) {
    aomega         <- lt_a_closeout(
                         mx = nMx,
                         Age = Age,
                         extrapLaw = extrapLaw,
                         extrapFrom = extrapFrom,
                         extrapFit = extrapFit,
                         ...)

    ax[N]          <- aomega
  } else {
    ax[N]          <- 1 / nMx[N]
  }
  #
  ax
}

#' UN a(x) estimates from either M(x), q(x), or both
#'
#' @description The UN a(x) formula uses Coale-Demeny for ages 0, and 1-4, values of 2.5 for ages 5-9 and 10-14, and the Greville formula for higher ages. In the original sources these are referred to as separation factors.
#'
#' @details a(x) for age 0 and age group 1-4 are based on Coale-Demeny {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}}-based lookup tables. If the main input is \code{nMx}, and if \code{IMR} is not given, we first approximate {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}} for the Coale-Demeny approach before applying the formula. The final a(x) value is closed out using the \code{lt_a_closeout()} method (reciprocal and Mortpak methods are deprecated). For nMx inputs this method is rather direct, but for {\ifelse{html}{\out{q<sub>X</sub>}}{\eqn{q_X}}} or l(x) inputs it is iterative. Age groups must be standard abridged.  No check on age groups are done.
#'
#' @param nMx numeric. Event exposure mortality rates.
#' @param nqx numeric.  Vector of age specific death probabilities in standard abridged age groups.
#' @param lx numeric. Vector of lifetable survivorship in standard abridged age groups.
#' @param a0rule character. Either \code{"ak"} (default) or \code{"cd"}.
#' @param IMR numeric. Optional. {\ifelse{html}{\out{q<sub>0</sub>}}{\eqn{q_0}}}, the death probability in first year of life, in case available separately.
#' @param AgeInt integer. Vector of age interval widths.
#' @param Sex character. \code{"m"}, \code{"f"} or \code{"b"} for male, female, or both.
#' @param region character. \code{"n"}, \code{"e"}, \code{"s"} or \code{"w"} for North, East, South, or West.
#' @param tol numeric. The tolerance for the qx-based iterative method. Default \code{.Machine$double.eps}.
#' @param maxit integer. The maximum number of iterations for the qx-based iterative method. Default 1000.
#' @param mod logical.  Whether or not to use Gerland's modification for ages 5-14. Default \code{TRUE}.
#' @param SRB numeric. The sex ratio at birth (boys/girls), default 1.05.
#' @param extrapLaw character. If extrapolating, which parametric mortality law should be invoked? Options include  \code{"Kannisto", "Kannisto_Makeham", "Makeham","Gompertz", "GGompertz", "Beard",	"Beard_Makeham", "Quadratic"}. Default \code{"Kannisto"}. See details.
#' @inheritParams lt_a_closeout
#'
#' @return nax average contribution to exposure of those dying in the interval.
#' @export
#'
#'
#' @references
#' \insertRef{greville1977short}{DemoTools}
#' \insertRef{un1982model}{DemoTools}
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertRef{mortpak1988}{DemoTools}
#' @examples
#' # example data from UN 1982 Model Life Tables for Developing Countries.
#' # first Latin American model table for males (p. 34).
#' Mx <- c(.23669,.04672,.00982,.00511,.00697,.01036,.01169,
#' 		.01332,.01528,.01757,.02092,.02517,.03225,.04241,.06056,
#' 		.08574,.11840,.16226,.23745)
#' ax <- c(0.330,1.352,2.500,2.500,2.633,2.586,2.528,2.528,
#' 		2.526,2.529,2.531,2.538,2.542,2.543,2.520,2.461,2.386,2.295,4.211)
#'
#' AgeInt     <- inferAgeIntAbr(vec = Mx)
#' Age <- int2age(AgeInt)
#' nAx1       <- lt_a_un(nMx = Mx,
#'                    Age = Age,
#' 		                AgeInt = AgeInt,
#' 		                a0rule = "cd",
#' 		                Sex = "m",
#'					          region = "w",
#' 		                mod = FALSE)
#' nAx2       <- lt_a_un(nMx = Mx,
#'                    Age = Age,
#' 		                AgeInt = AgeInt,
#' 		                a0rule = "cd",
#' 		                Sex = "m",
#'					          region = "w",
#' 		                mod = TRUE)
#' # this is acceptable...
#' round(nAx2,3) - ax # only different in ages 5-9 and 10-14, and last two ages
#' # ignore open age, which is treated differently
#' N <- length(ax)
#' # default unit test...
#' stopifnot(all(round(nAx1[Age<80],3) - ax[Age<80] == 0)) # spot on
#' 
#' # another example:
#' 
#' nMx <- c(0.11621,0.02268,0.00409,0.00212,0.00295,
#'          0.00418,0.00509,0.00609,0.00714,0.00808,
#'          0.00971,0.0125,0.0175,0.02551,0.03809,
#'          0.05595,0.0809,0.15353,0.2557)
#' AgeInt <- inferAgeIntAbr(vec = nMx)
#' Age    <- int2age(AgeInt)
#' nAx1   <- lt_a_un(nMx, 
#'                   Age = Age,
#'                   AgeInt = AgeInt, 
#'                   a0rule = "ak",
#'                   mod = TRUE, 
#'                   closeout = TRUE)
#' 
#' nqx    <- lt_id_ma_q(nMx = nMx, 
#'                      nax = nAx1, 
#'                      AgeInt = AgeInt, 
#'                      closeout = FALSE)
#' nAx2 <- lt_a_un(nqx = nqx, 
#'                 Age = Age,
#'                 AgeInt = AgeInt, 
#'                 a0rule = "ak", 
#'                 mod= TRUE, 
#'                 closeout = TRUE)
#' 
#' stopifnot(all(abs(nAx1[Age<75] - nAx2[Age<75]) < 1e-7))
#'

lt_a_un <- function(nMx,
                 nqx,
                 lx,
                 IMR = NA,
                 Age,
                 AgeInt,
                 a0rule = "ak",
                 Sex = "m",
                 region = "w",
                 SRB = 1.05,
                 tol = .Machine$double.eps,
                 maxit = 1e3,
                 mod = TRUE,
                 extrapLaw = "kannisto",
                 extrapFrom = max(Age),
                 extrapFit = Age[Age >= 60],
                 ...) {

  stopifnot(!missing(nqx) | !missing(nMx))
  smsq    <- 99999
  
  Sex      <- match.arg(Sex, choices = c("m","f","b"))
  a0rule   <- match.arg(a0rule, choices =  c("ak","cd"))
  extrapLaw      <- match.arg(extrapLaw, choices = c("kannisto",
                                         "kannisto_makeham",
                                         "makeham",
                                         "gompertz",
                                         "ggompertz",
                                         "beard",
                                         "beard_makeham",
                                         "quadratic"
  ))
  region   <- match.arg(region, choices =c("w","n","s","e"))
  if (missing(nqx) & !missing(lx)) {
    nqx <- lt_id_l_d(lx) / lx
  }
  # now we have either nqx or nMx

  if (missing(nqx) & !missing(nMx)) {
    # UN (1982) p 31
    # http://www.un.org/esa/population/publications/Model_Life_Tables/Model_Life_Tables.htm
    #		For ages 15 and over, the expression for nQx is derived
    #		from Greville" as ,nax, = 2.5 - (25/12) (nmx) - k), where
    #		k = 1/10 log(nmx+5/nmx-5). For ages 5 and 10, nQx = 2.5
    #		and for ages under 5, nQx values from the Coale and
    #		Demeny West region relationships are used."
    axi <- lt_id_morq_a_greville(
      nMx = nMx,
      Age  = Age,
      AgeInt = AgeInt,
      a0rule = a0rule,
      IMR = IMR,
      Sex = Sex,
      region = region,
      mod  = mod,
      SRB = SRB,
      extrapLaw  = extrapLaw,
      extrapFrom = extrapFrom,
      extrapFit = extrapFit,
      ...
    )
  }
 
  if (!missing(nqx) & missing(nMx)) {
    # UN (1982) p 31
    # http://www.un.org/esa/population/publications/Model_Life_Tables/Model_Life_Tables.htm
    #		With nqx as input, the procedure is identical, except
    #		that an iterative procedure is used to find the nmx and nqx
    #		values consistent with the given nqx and with the Greville
    #		expression.
    axi <- lt_id_morq_a_greville(
      nqx = nqx,
      IMR = nqx[1],
      Age = Age,
      AgeInt = AgeInt,
      a0rule = a0rule,
      Sex = Sex,
      region = region,
      mod = mod,
      SRB = SRB,
      extrapLaw = extrapLaw,
      extrapFrom = extrapFrom,
      extrapFit = extrapFit,
      ...
    )
    # no longer needed
    closeout <- FALSE
    for (i in 1:maxit) {
      mxi   <- lt_id_qa_m(
                    nqx = nqx,
                    nax = axi,
                    AgeInt = AgeInt)
      axi   <- lt_id_morq_a_greville(
        nMx = mxi,
        IMR = nqx[1],
        Age = Age,
        AgeInt = AgeInt,
        a0rule = a0rule,
        Sex = Sex,
        region = region,
        mod = mod,
        # no need to redo extrap in here
        ...
      )
      qxnew <-
        lt_id_ma_q(
          nMx = mxi,
          nax = axi,
          AgeInt = AgeInt,
          IMR = IMR,
          closeout = closeout
        )
      smsq  <- sum((qxnew - nqx) ^ 2)
      if (smsq < tol) {
        break
      }
    }
    # no need for approximate a0 and 4a1 values
    # one last time for nMx
    nMx <- lt_id_qa_m(nqx = nqx,
                   nax = axi,
                   AgeInt = AgeInt)
  }
  # if both given, then we have ax via identity:
  if (!missing(nqx) & !missing(nMx)) {
    axi <- lt_id_qm_a(nqx = nqx,
                   nMx = nMx,
                   AgeInt = AgeInt)
  }

  #	if (closeout == "mortpak" & sum(AgeInt) < 100){
  #		N      <- length(axi)
  #		aomega <- aomegaMORTPAK(mx_or_qx                                                                       = nMx, qind = FALSE)
  #		axi[N] <- aomega
  #	}
  #
  # closeout
  N    <- length(axi)
  if (max(Age) <= 125) {
    aomega         <- lt_a_closeout(
      mx = nMx,
      Age = Age,
      extrapLaw = extrapLaw,
      extrapFrom = extrapFrom,
      extrapFit = extrapFit,
      ...
    )
    axi[N] <- aomega
  } else {
    axi[N] <- 1 / nMx[N]
  }
  # patch
  ind <- is.nan(axi) | axi < 0
  if (any(ind)){
    axi[ind] <- AgeInt[ind] / 2
  }
  
  # if mx, qx, or both are given, then by now we have ax
  axi
}

#' Life expectancy in the open age group.
#'
#' @description Get an estimate of life expectancy in the open age group.
#' @details This method estimates life expectancy in the open age group by fitting one of several potential old-age parametric mortality models, extrapolating rates to age 130, then backing out the implied remaining life expectancy in the open age group. This function replaces \code{aomegaMORTPAK()}.
#' @inheritParams lt_rule_m_extrapolate
#' @param Age integer. A vector of ages of the lower integer bound of the age classes.
#' @param extrapFrom integer. Age from which to impute extrapolated mortality.
#' @param extrapFit integer vector. Ages to include in model fitting. Defaults to all ages \code{>          =60}.
#' @param extrapLaw character. The following options are available: \itemize{
#'   \item{\code{"kannisto"}} -- The Kannisto model;
#'   \item{\code{"kannisto_makeham"}} -- The Kannisto-Makeham model;
#'   \item{\code{"gompertz"}} -- The Gompertz model;
#'   \item{\code{"ggompertz"}} -- The Gamma-Gompertz model;
#'   \item{\code{"makeham"}} -- The Makeham model;
#'   \item{\code{"beard"}} -- The Beard model;
#'   \item{\code{"beard_makeham"}} -- The Beard-Makeham model;
#'   \item{\code{"quadratic"}} -- The Quadratic model.
#'   }
#' @return life expectancy in the open age group
#' @seealso
#' \code{\link[DemoTools]{lt_rule_m_extrapolate}}
#' @export
#' @examples
#' nMx <- c(0.12846,0.02477,0.00603,0.0034,
#'  		0.00417,0.00513,0.00581,0.00645,0.00725,
#'  		0.00813,0.00913,0.01199,0.01647,
#'  		0.0256,0.04047,0.06624,0.10638,0.19611)
#' Age <- c(0,1,seq(5,80,by =5))
#'
#'
#' lt_a_closeout(nMx,Age,"kannisto")
#' lt_a_closeout(nMx,Age,"kannisto_makeham")
#' lt_a_closeout(nMx,Age,"makeham")
#' lt_a_closeout(nMx,Age,"gompertz")
#' lt_a_closeout(nMx,Age,"ggompertz")
#' lt_a_closeout(nMx,Age,extrapLaw ="beard")
#' lt_a_closeout(nMx,Age,"beard_makeham")
#' lt_a_closeout(nMx,Age,"quadratic")

lt_a_closeout <- function(mx,
                          Age,
                          extrapLaw = "kannisto",
                          extrapFrom = max(Age),
                          extrapFit = Age[Age >= 40],
                          ...) {


  extrapLaw      <- match.arg(extrapLaw, choices = c("kannisto",
                                         "kannisto_makeham",
                                         "makeham",
                                         "gompertz",
                                         "ggompertz",
                                         "beard",
                                         "beard_makeham",
                                         "quadratic"
  ))

  OA        <- max(Age)
  x_extr    <- seq(OA, 130, by = .1)

  Mxnew     <- lt_rule_m_extrapolate(
    mx = mx,
    x = Age,
    x_fit = extrapFit,
    x_extr = x_extr,
    law = extrapLaw,
    ...
  )
  mmxx      <- Mxnew$values
  mx        <- mmxx[names2age(mmxx) >= OA]

  # acceptable approximation for small intervals
  lx        <- c(1, exp(-cumsum(mx / 10)), 0)
  # divide by 10 (interval) and 2 (avg), so 20.
  sum(shift.vector(lx, -1) + lx) / 20

}

#' wrapper to invoke PAS or UN ax methods given qx or mx
#' @description Given either mx or qx, call either the \code{lt_a_un()} or \code{lt_a_pas()} functions.
#' @inheritParams lt_a_un
#' @param axmethod character. Either \code{"pas"} or \code{"un"}
#' @param OAG logical. Whether or not the last element of \code{nMx} is the open age group Default \code{TRUE}.
#' @return nax average contribution to exposure of those dying in the interval.
#' @references
#' \insertRef{greville1977short}{DemoTools}
#' \insertRef{un1982model}{DemoTools}
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertRef{mortpak1988}{DemoTools}
#' \insertRef{united1983manual}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#' @export

lt_id_morq_a <- function(nMx,
                      nqx,
                      axmethod = "pas",
                      Age,
                      AgeInt,
                      a0rule = "ak",
                      IMR = NA,
                      Sex ="m",
                      region,
                      OAG = TRUE,
                      mod = TRUE,
                      SRB = SRB,
                      extrapLaw = "kannisto",
                      extrapFrom = max(Age),
                      extrapFit = Age[Age >= 60],
                      ...) {
  Sex      <- match.arg(Sex, choices = c("m","f","b"))
  a0rule   <- match.arg(a0rule, choices = c("ak","cd"))
  extrapLaw      <- match.arg(extrapLaw, choices = c("kannisto",
                                         "kannisto_makeham",
                                         "makeham",
                                         "gompertz",
                                         "ggompertz",
                                         "beard",
                                         "beard_makeham",
                                         "quadratic"
  ))
  region   <- match.arg(region, choices = c("w","n","s","e"))
  
  N <- length(AgeInt)
  if (is.na(AgeInt[N]) | is.infinite(AgeInt[N])) {
    AgeInt[N] <- AgeInt[N - 1]
  }
  if (axmethod == "pas") {
    # what if only qx was given?
    if (missing(nMx)) {
      nAx       <- lt_a_pas(
                     nMx = nqx,
                     AgeInt = AgeInt,
                     a0rule = a0rule,
                     IMR = IMR,
                     Sex = Sex,
                     region = region,
                     OAG = OAG,
                     SRB = SRB)

    } else {
      # if nMx avail, then Open age group
      # closed according to convention.
      nAx       <- lt_a_pas(
                     nMx = nMx,
                     AgeInt = AgeInt,
                     a0rule = a0rule,
                     IMR = IMR,
                     Sex = Sex,
                     region = region,
                     OAG = TRUE,
                     SRB = SRB)
    }
    
  }
  if (axmethod == "un") {
    # UN method just CD west for now, so no region arg
    # no sense calling Abacus here because it gets called later if necessary
    if (missing(nMx)) {
      #fakenMx   <- nqx
      nAx       <- lt_a_un(
                     nqx = nqx,
                     Age = Age,
                     AgeInt = AgeInt,
                     a0rule = a0rule,
                     IMR = IMR,
                     Sex = Sex,
                     region = region,
                     mod = mod,
                     SRB = SRB,
                     extrapLaw = extrapLaw,
                     extrapFrom = extrapFrom,
                     extrapFit = extrapFit,
                     ...)

    } else {
      nAx       <- lt_a_un(
                     nMx = nMx,
                     Age = Age,
                     AgeInt = AgeInt,
                     a0rule = a0rule,
                     IMR = IMR,
                     Sex = Sex,
                     region = region,
                     mod = mod,
                     SRB = SRB,
                     extrapLaw = extrapLaw,
                     extrapFrom = extrapFrom,
                     extrapFit = extrapFit,
                     ...)
    }
  }
  
 
  
  # TR: shall we do ak patch just here at the end?
  # the alternative would be to mesh it in everywhere a0 happens.
  # ergo lt_rule_a0() as a new function 
  
  nAx
}

# deprecated.
##' Life expectancy in the open age group.
##'
##' @description Get the Abacus Mortpak estimate of life expectancy in the open age group.
##' @details Since the Mortpak lifetable just goes to age 100, it only makes sense to call this function if the data have a lower open age group.
##' If the data go to 100 or higher, there is no apparent advantage to closing out with this function. Specify the entire nMx schedule, in standard abridged ages.
##' @details The estimate will be the same for males and females.
##' @param mx_or_qx numeric.  Vector of mortality rates or probabilities in standard abridged age classes.
##' @param qind logical. Default \code{FALSE} (implying Mx used). \code{TRUE} means qx was given.
##' @return Open age groups life expectancy.
##' @references
##' \insertRef{mortpak1988}{DemoTools}
##' @export
#
#aomegaMORTPAK <- function(mx_or_qx,qind                                                                    =FALSE){
#	OA  <- length(mx_or_qx) * 5 - 10
#	if (qind){
#		OUT <- AbacusLIFTB_wrap(qx = mx_or_qx, mx_ind = FALSE, OAnew = OA, Sex = "m")
#	} else {
#		OUT <- AbacusLIFTB_wrap(Mx = mx_or_qx, OAnew = OA, Sex = "m")
#	}
#	OUT[nrow(OUT),"ax"]
#}
