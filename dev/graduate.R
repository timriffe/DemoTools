
# TODO make wrapper for graduation methods.

# sprague, beers, splitMono, pclm (import)

library(devtools)
library(here)

load_all()

library(ungroup)
#library(DemoTools)

#' wrapper for the PCLM method of splitting binned counts
#' 
#' @description This is exactly the function \code{pclm()} from the \code{ungroup} package, except with arguments using standard \code{DemoTools} argument names.
#' @details The PCLM method can also be used to graduate rates using an offset if both numerators and denominators are available. The denominators must be in single ages already. Until such time as \code{graduate_rates_pclm()} is written, if denominators are also in grouped ages, then must first be split, and this is therefore a two-step process for users. This method can be used to redistribute counts in the open age group if OAnew gives sufficient space, otherwise resutls may vary.
#' 
#' @inheritParams graduate
#' @param ... further arguments passed to \code{ungroup::pclm()}
#' 
#' @references 
#' \insertRef{pascariu2018ungroup}{DemoTools}
#' \insertRef{rizzi2015efficient}{DemoTools}
#' 
#' @importFrom ungroup pclm
#' @examples 
#' a5 <- seq(0,100,by=5)
#' p5 <- pop5_mat[, 1]
#' p1 <- graduate_counts_pclm(p5, Age = a5)
#' p1s <- sprague(p5,a5)
#' plot(a5, p5/5, type = "s",xlim=c(40,60),ylim=c(2000,4000))
#' lines(0:100, p1, lwd = 2, col = "red")
#' lines(0:100, p1s, lwd = 1, col = "blue",lty="8282")
graduate_counts_pclm <- function(Value, Age, OAnew = max(Age), ...){
  nlast    <- OAnew - max(Age) + 1
  A        <- ungroup::pclm(x = Age, y = Value, nlast = nlast, ...)
  B        <- A$fitted
  names(B) <- min(Age):OAnew
  B
}
args(sprague)
args(beers)
args(splitMono)

Value <- pop5_mat[, 1]
Value <- c(10000,44170,Value[-1])
Age   <- sort(c(1,seq(0,100,by=5)))

graduate(Value, Age, method = "sprague")
graduate(Value, Age, method = "sprague", keep0=FALSE)

graduate(Value, Age, method = "beers(ord)")
graduate(Value, Age, method = "beers(ord)", keep0=TRUE)
graduate(Value, Age, method = "beers(ord)", keep0=TRUE, johnson = TRUE)

graduate(Value, Age, method = "beers(mod)")
graduate(Value, Age, method = "beers(mod)", keep0=TRUE)
graduate(Value, Age, method = "beers(mod)", keep0=TRUE, johnson = TRUE)

graduate(Value, Age, method = "pclm")
graduate(Value, Age, method = "pclm", keep0=TRUE)

graduate(Value, Age, method = "mono")
graduate(Value, Age, method = "mono", keep0=TRUE)

graduate(Value, Age, method = "uniform")

graduate <- function(
  Value, 
  Age, 
  AgeInt = age2int(Age), 
  OAG = TRUE, 
  OAnew = max(Age), 
  method = c("sprague","beers(ord)","beers(mod)","pclm","mono","uniform"), 
  keep0 = FALSE, ...){
  
  # validate method choice
  method <- match.arg(method)
  
  # TR: those methods that require AgeInt may depend on final value not being NA
  # even if it's an open age, in essence, keep all value inside this same "single"
  # age. NA coding isn't the best choice here, but we anticipate this and assign
  # 1 to AgeInt[length(AgeInt)] IFF OAG & max(Age) == OAnew 
  # This is inconsequential for those methods that don't use AgeInt
  N <- length(AgeInt)
  if (OAG & is.na(AgeInt[N])){
    nlast     <- OAnew - max(Age) + 1
    AgeInt[N] <- nlast
  }
  
  # Sprague in strict 5-year age groups
  if (method == "sprague"){
    out <- sprague(Value, Age = Age, OAG = OAG)
  }
  
  # Beers in strict 5-year age groups
  if (grepl("beers", method)){
    if (grepl("ord", method)){
      out <- beers(Value, Age = Age, OAG = OAG, method = "ord", ...)
    }
    if (grepl("mod", method)){
      out <- beers(Value, Age = Age, OAG = OAG, method = "mod", ...)
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
  if (method == "uniform"){
    OAvalue <- OAnew - max(Age) + 1
    out <- splitUniform(
             Value = Value, 
             Age = Age, 
             AgeInt = AgeInt, 
             OAG = OAG,
             OAvalue = OAvalue) # OAvalue only if OAG and
                                # extrapolation desired?
  }
  
  # Mono respects irregular intervals
  if (method == "mono"){
    out <- splitMono(
      Value = Value, 
      Age = Age, 
      AgeInt = AgeInt, 
      OAG = OAG) # doesn't allow for extension?
    # actually it does IFF !OAG & !is.na(AgeInt[length(AgeInt)])
  }
  
  # PCLM respects irregular intervals
  if (method == "pclm"){
    out <- graduate_counts_pclm(Value = Value, Age = Age, OAnew = OAnew, ...)
  }
  
  # handle infant age for sprague and beers, if required.
  if (keep0){
    if (Age[1] == 0 & AgeInt[1] == 1){
      V0       <- Value[1]
      V5       <- sum(Value[Age < 5])
      V4       <- V5 - V0
      a1       <- names2age(out)
      ind      <- a1 < 5 & a1 > 0
      out[ind] <- rescale.vector(out[ind],scale=V4)
      out[1]   <- V0
    }
  }
  
  out
}





