# TODO: add grabill

#' wrapper for \code{ungroup::pclm} method of splitting binned counts
#'
#' @description This is exactly the function \code{pclm()} from the \code{ungroup} package, except with arguments using standard \code{DemoTools} argument names.
#' @details The PCLM method can also be used to graduate rates using an offset if both numerators and denominators are available. In this case \code{Value} is the event count and \code{offset} is person years of exposure. The denominator must match the length of \code{Value} or else the length of the final single age result \code{length(min(Age):OAnew)}.  This method can be used to redistribute counts in the open age group if \code{OAnew} gives sufficient space. Likewise, it can give a rate extrapolation beyond the open age.
#'
#' @inheritParams graduate
#' @param ... further arguments passed to \code{ungroup::pclm()}
#'
#' @references
#' \insertRef{pascariu2018ungroup}{DemoTools}
#' \insertRef{rizzi2015efficient}{DemoTools}
#'
#' @importFrom ungroup pclm
#' @export
#' @seealso \code{\link[ungroup]{pclm}}
#' @examples
#' a5  <- seq(0,100,by=5)
#' p5  <- pop5_mat[, 1]
#' p1  <- graduate_pclm(Value = p5, Age = a5)
#' p1s <- graduate_sprague(Value = p5, Age = a5)
#' \dontrun{
#' plot(a5, p5/5, type = "s",xlim=c(40,60),ylim=c(2000,4000))
#' lines(0:100, p1, lwd = 2, col = "red")
#' lines(0:100, p1s, lwd = 1, col = "blue",lty="8282")
#' }
graduate_pclm <- function(Value, Age, OAnew = max(Age), ...) {
  nlast    <- OAnew - max(Age) + 1
  a1       <- min(Age):OAnew
  DOTS     <- list(...)
  if ("offset" %in% names(DOTS)) {
    # offset could be one or another thing..
    lo     <- length(DOTS$offset)
    o1     <- length(a1) == lo
    o5     <- length(Value) == lo
    stopifnot(o1 | o5)
  }
  
  A        <- pclm(x = Age, y = Value, nlast = nlast, ...)
  B        <- A$fitted
  names(B) <- min(Age):OAnew
  B
}


#' Graduate grouped data
#'
#' @description A wrapper function for several graduation methods, primarily for count data (\code{"sprague"}, \code{"beers(ord)"}, \code{"beers(mod)"}, \code{"grabill"}, \code{"mono"} (Monotonic spline), \code{"uniform"}, \code{"pclm"}), but also with one (\code{"pclm"}) with an option for graduating rates if both event counts and population at risk are available.
#'
#' @details \code{"sprague"}, \code{"beers(ord)"}, \code{"beers(mod)"} methods require original data to be in uniform five-year age groups. If they are not (for example, the infant group is separate) then they are grouped to uniform width prior to splitting. If you want to keep the original infant count in output, then specify \code{keep0 = TRUE}. In this case, it is imputed, and ages 1-4 are rescaled, which may introduce a discontinuity in results from age 4 to 5. \code{keep0 = TRUE} may also be desired along with \code{method = "pclm"}.
#'
#' Some methods are constrained, others not, and others are optionally constrained. If this is required, then this function can be followed up with \code{rescaleAgeGroups()}, which may have the effect of breaking continuity in smooth results. This is inconsequential for downstream demography, but if this aesthetic side effect is undesired, then try one of the constrained methods: \code{"sprague"}, \code{"mono"}, \code{"pclm"} (with \code{control = list(lambda = 1/1e7) specified} or similar).
#'
#' Beers may either be ordinary \code{"beers(ord)"} or modified \code{"beers(mod)"}, and either can pass on the optional argument \code{johnson = TRUE} if desired (this has a different distribution pattern for young ages, \code{FALSE} by default). If \code{method = "beers"} is given, then \code{"beers(ord)"} is used.
#'
#' This wrapper standardizes some inconsistencies in how open ages are dealt with. For example, with the \code{"pclm"} method, the last age group can be redistributed over a specified interval implied by increase \code{OAnew} beyond the range of \code{Age}. To get this same behavior from \code{"mono"}, or \code{"uniform"} specify \code{OAG = FALSE} along with an appropriately high \code{OAnew} (or integer final value of \code{AgeInt}.
#'
#' \code{OAnew} cannot be higher than \code{max(Age)+4} for \code{"sprague"} or \code{"beers"} methods. For \code{"uniform","mono","pclm"} it can be higher than this, and in each case the open age group is completely redistributed within this range, meaning it's not really open anymore.
#'
#' @param Value numeric vector, presumably counts in grouped ages
#' @param Age integer vector, lower bounds of age groups
#' @param AgeInt integer vector, age interval widths
#' @param OAG logical, default = \code{TRUE} is the final age group open?
#' @param OAnew integer, optional new open age, higher than \code{max(Age)}. See details.
#' @param method character, either \code{"sprague"}, \code{"beers(ord)")}, \code{"beers(mod)")}, \code{"mono")}, \code{"uniform")}, or \code{"pclm"}
#' @param keep0 logical. Default \code{FALSE}. If available, should the value in the infant age group be maintained, and ages 1-4 constrained?
#' @param ... extra arguments passed to \code{graduate_beers()} or \code{graduate_pclm()}
#' @seealso \code{\link{graduate_sprague}}, \code{\link{graduate_beers}}, \code{\link{graduate_uniform}}, \code{\link{graduate_mono}}, \code{\link{graduate_pclm}}, \code{\link{graduate_grabill}} 
#' @export
#' @references
#' \insertRef{pascariu2018ungroup}{DemoTools}
#' \insertRef{rizzi2015efficient}{DemoTools}
#' \insertRef{sprague1880explanation}{DemoTools}
#' \insertRef{shryock1973methods}{DemoTools}
#' \insertRef{siegel2004methods}{DemoTools}
#' \insertRef{beers1945modified}{DemoTools}
#'
#' @examples
#' Value <- pop5_mat[, 1]
#' Value <- c(10000,44170,Value[-1])
#' Age   <- sort(c(1,seq(0,100,by=5)))
#'
#' graduate(Value, Age, method = "sprague")
#' graduate(Value, Age, method = "sprague", keep0=FALSE)
#'
#' graduate(Value, Age, method = "beers(ord)")
#' graduate(Value, Age, method = "beers(ord)", keep0=TRUE)
#' graduate(Value, Age, method = "beers(ord)", keep0=TRUE, johnson = TRUE)
#'
#' graduate(Value, Age, method = "beers(mod)")
#' graduate(Value, Age, method = "beers(mod)", keep0=TRUE)
#' graduate(Value, Age, method = "beers(mod)", keep0=TRUE, johnson = TRUE)
#'
#' graduate(Value, Age, method = "mono")
#' graduate(Value, Age, method = "mono", keep0=TRUE)
#'
#' graduate(Value, Age, method = "uniform")
#'
#' graduate(Value, Age, method = "pclm")
#' graduate(Value, Age, method = "pclm", keep0=TRUE)
#' # pclm can also graduate rates if both
#' # numerators and denominators are on hand:
#' Exposures <- c(100958,466275,624134,559559,446736,370653,301862,249409,
#'                247473,223014,172260,149338,127242,105715,79614,53660,
#'                31021,16805,8000,4000,2000,1000)
#'
#' Deaths <- c(8674,1592,618,411,755,1098,1100,1357,
#'             1335,3257,2200,4023,2167,4578,2956,4212,
#'             2887,2351,1500,900,500,300)
#' Age    <- c(0, 1, seq(5, 100, by = 5))
#' AgeInt <- c(diff(Age), NA)
#'
#' # exclude infants for better fit.
#' mx    <- graduate(
#'            Value = Deaths[-1], Age = Age[-1],
#'            AgeInt = AgeInt[-1], OAG = TRUE,
#'            OAnew = 110, offset = Exposures[-1],
#'            method = "pclm")
#' mx_sm <- graduate(
#'            Value = Deaths[-1], Age = Age[-1],
#'            AgeInt = AgeInt[-1], OAG = TRUE,
#'            OAnew = 110, offset = Exposures[-1],
#'            method = "pclm", control = list(lambda = 1e7))
#'
#' \dontrun{
#' plot(Age,
#'      Deaths / Exposures,
#'      type = 's', log = 'y',
#'      main = "Underlying data have differential heaping on 0s and 5s")
#' lines(1:110, mx)
#' lines(1:110, mx_sm, col = "blue")
#' legend("bottomright",
#'        col = c("black","blue"),
#'        lty = c(1, 1),
#'        legend = c("lambda optimized (almost constrained)",
#'                   "higher lambda = smoother")
#'        )
#'   }

graduate <- function(Value,
                     Age,
                     AgeInt = age2int(Age),
                     OAG = TRUE,
                     OAnew = max(Age),
                     method = c("sprague",
                                "beers(ord)",
                                "beers(mod)",
                                "grabill",
                                "pclm",
                                "mono",
                                "uniform"),
                     keep0 = FALSE,
                     ...) {
  method <- tolower(method)
  if (method == "beers") {
    method == "beers(ord)"
  }
  # validate method choice
  method <- match.arg(method)
  
  # TR: those methods that require AgeInt may depend on final value not being NA
  # even if it's an open age, in essence, keep all value inside this same "single"
  # age. NA coding isn't the best choice here, but we anticipate this and assign
  # 1 to AgeInt[length(AgeInt)] IFF OAG & max(Age) == OAnew
  # This is inconsequential for those methods that don't use AgeInt
  N <- length(AgeInt)
  if (OAG & is.na(AgeInt[N])) {
    nlast     <- OAnew - max(Age) + 1
    AgeInt[N] <- nlast
  }
  
  # Sprague in strict 5-year age groups
  if (method == "sprague") {
    out <- graduate_sprague(Value, Age = Age, OAG = OAG)
  }
  
  # Grabill in strict 5-year age groups
  if (method == "grabill") {
    out <- graduate_grabill(Value, Age = Age, OAG = OAG)
  }
  # Beers in strict 5-year age groups
  if (grepl("beers", method)) {
    if (grepl("ord", method)) {
      out <- graduate_beers(Value,
                   Age = Age,
                   OAG = OAG,
                   method = "ord",
                   ...)
    }
    if (grepl("mod", method)) {
      out <- graduate_beers(Value,
                   Age = Age,
                   OAG = OAG,
                   method = "mod",
                   ...)
    }
  }
  
  # inconsistent ways of dealing with top age group.
  # if !OAG, then there is potential for inconsistency
  # between (min(Age) + sum(AgeInt)) and OAnew,
  # namely, OAnew can increase or decrease this age IFF
  # OAG, but if there is an inconsistency then we need
  # to either throw an error or declare a preference.
  # TR: preference should go to AgeInt, because it refers
  # to data observations, whereas OAnew refers to a desired
  # change. So actually no change is needed?
  
  
  # Uniform respects irregular intervals
  if (method == "uniform") {
    OAvalue <- OAnew - max(Age) + 1
    out <- graduate_uniform(
      Value = Value,
      Age = Age,
      AgeInt = AgeInt,
      OAG = OAG,
      OAvalue = OAvalue
    ) # OAvalue only if OAG and
    # extrapolation desired?
  }
  
  # Mono respects irregular intervals
  if (method == "mono") {
    out <- graduate_mono(
      Value = Value,
      Age = Age,
      AgeInt = AgeInt,
      OAG = OAG
    ) # doesn't allow for extension?
    # actually it does IFF !OAG & !is.na(AgeInt[length(AgeInt)])
  }
  
  # PCLM respects irregular intervals
  if (method == "pclm") {
    out <- graduate_pclm(Value = Value,
                         Age = Age,
                         OAnew = OAnew,
                         ...)
  }
  
  # handle infant age for sprague and beers, if required.
  if (keep0) {
    if (Age[1] == 0 & AgeInt[1] == 1) {
      V0       <- Value[1]
      V5       <- sum(Value[Age < 5])
      V4       <- V5 - V0
      a1       <- names2age(out)
      ind      <- a1 < 5 & a1 > 0
      out[ind] <- rescale.vector(out[ind], scale = V4)
      out[1]   <- V0
    }
  }
  
  # TR: TODO detect negatives. Have default option to re
  
  out
}
