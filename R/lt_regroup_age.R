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
#' @param ... optional args, not currently used.
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


# TODO this needs to be speed profiled. Why is pclm() slow? Is it just my machine?

# A life table by single year of age obtained by graduating the abridged lt using ungroup package

#' create a life table by single year of age by graduating an abridged life table
#' @description Computes single year of age life table by graduating the mortality schedule of an abridged life table, using the `ungroup::pclm()` to ungroup binned count data. Returns complete single-age lifetable.
#' @details Similar to `lt_abridged()` details, forthcoming. 
#' @inheritParams lt_abridged
#' @param ... optional arguments passed to `pclm()`. For example, if you pass an expicit `lambda` parameter via the `control` argument, you can speed up estimation
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
#'  Mx <- c(.23669,.04672,.00982,.00511,.00697,.01036,.01169,
#'          .01332,.01528,.01757,.02092,.02517,.03225,.04241,.06056,
#'          .08574,.11840,.16226,.23745)
#'  Age = c(0,1,seq(5,85,by=5))
#'  AgeInt <- inferAgeIntAbr(vec = Mx)
#'  LTabr <- lt_abridged(nMx = Mx,
#'                       Age = Age, 
#'                       axmethod = "un",
#'                       Sex = "m",
#'                       mod = TRUE)
#'  
#'  LT1 <- lt_abridged2single(nMx = Mx,
#'                     Age = Age, 
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
lt_abridged2single <- function(
  Deaths = NULL, 
  Exposures = NULL, 
  nMx = NULL, 
  nqx = NULL, 
  lx = NULL,
  Age,
  radix = 1e5,
  axmethod = "un",
  a0rule = "ak", 
  Sex = "m",
  region = "w",
  IMR = NA,
  mod = TRUE,
  SRB = 1.05,
  OAG = TRUE,
  OAnew = max(Age),
  extrapLaw = NULL,
  extrapFrom = max(Age),
  extrapFit = NULL,
  ...) {
  
  stopifnot(is_abridged(Age))
  NN <- length(Age)
  #stopifnot(length(nMx) == NN)
  
  if (!is.null(extrapLaw)){
    extrapLaw      <- tolower(extrapLaw)
    extrapLaw      <- match.arg(extrapLaw, choices = c("kannisto",
                                                       "kannisto_makeham",
                                                       "makeham",
                                                       "gompertz",
                                                       "ggompertz",
                                                       "beard",
                                                       "beard_makeham",
                                                       "quadratic"
    ))
  } else {
    extrapLaw <- ifelse(max(Age)>=90, "kannisto","makeham")
  }
  if (is.null(extrapFit)){
    maxAclosed <- ifelse(OAG, Age[which.max(Age)-1],max(Age))
    if (maxAclosed < 85){
      extrapFit  <- Age[Age >= (maxAclosed - 20) & Age <= maxAclosed]
    } else {
      extrapFit  <- Age[Age >= 60 & Age <= maxAclosed]
    }
  } else {
    stopifnot(all(extrapFit %in% Age))
  }
  
  # first extend the abridged life table to OAG = 130 with a big radix so that we don't lose info later when rounding ndx and nLx to integers
  lt_abr <- lt_abridged(Deaths = Deaths, 
                        Exposures = Exposures, 
                        nMx = nMx, 
                        nqx = nqx, 
                        lx = lx, 
                        Age = Age,
                        Sex = Sex, 
                        radix = 1e8, 
                        axmethod = axmethod, 
                        a0rule = a0rule,
                        region = region, 
                        IMR = IMR, 
                        mod = mod, 
                        SRB = SRB, 
                        OAG = OAG, 
                        OAnew = 130,
                        extrapLaw = extrapLaw, 
                        extrapFrom = extrapFrom, 
                        extrapFit = extrapFit)
  
  # use pclm to ungroup to single year of age from 1 to 129
  # need to round ndx and nLx since pclm doesn't perform with values bw 0 and 1
  ndx <- round(lt_abr$ndx)
  nLx <- round(lt_abr$nLx)
  ind <- lt_abr$Age >= 1 & lt_abr$Age <= 125 & ndx>0 & nLx>0
  
  # TR: removed ... because in practice we were passing in a large
  # set of ... indirectly that aren't recognized in pclm
  M <- pclm(x      = lt_abr$Age[ind],
            y      = ndx[ind],
            nlast  = 5,
            offset = nLx[ind],
            ...)
  
  # splice original 1M0 with fitted 1Mx and momega from extended abridged LT
  M <- c(lt_abr$nMx[1], M$fitted)
  
  # TR: handle closeout nMx as well. Should depend on OAnew and Age to 
  # a certain extent.
  
  # redefine Age and extrapFit for single year ages and new maxage
  a1        <- 1:length(M) - 1
  extrapFit  <- a1[a1 >= min(extrapFit, (max(Age)-20)) & a1 <= max(Age)] 
  # always refit from 110 even if extrapFrom > 110
  extrapFrom <- min(max(Age), 110)
  
  # compute life table columns from single year mx
  LT <- lt_single_mx(nMx = M, 
                     Age = a1, 
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

#' calculate an abidged or single age lifetable from abridged or sinlge age data
#' @description This is a wrapper around the other lifetable utilities. We start with either `nMx`, `nqx`, or `lx` in single or abridged ages, and returns a full lifetable in either single or abridged ages. All optional arguments of `lt_abridged()` or `lt_single*()` can be passed in, for instance the `nax` assumptions or the extrapolation arguments.
#' 
#' @param nMx_or_nqx_or_lx numeric vector of either `nMx`, `nqx`, or `lx`
#' @param type character, which variable is `x`?, either `"m"`, `"q"`, or `"l"`. Default `"m"`
#' @param Age integer vector of the lower age bounds of `x`
#' @param Sex character, `"m"`, `"f"`, or `"b"`.
#' @param Single logical, do we want output in single ages?
#' @param ... optional arguments passed to `lt_abridged()` or `lt_single*()` 
#' @export

lt_ambiguous <- function(nMx_or_nqx_or_lx = NULL, 
                         type = "m",
                         Age = NULL, 
                         Sex = NULL, 
                         Single = FALSE,
                         ...){
  
  #extras <- list(...)
  
  xx <- nMx_or_nqx_or_lx
  # TR: adds flexibility when specifying type to reduce user errors
  type                  <- tolower(type)
  possible_types        <- c("m","m","m","q","q","q","l","l")
  names(possible_types) <- c("m","mx","nmx","q","qx","nqx","l","lx")
  stopifnot(type %in% names(possible_types) )
  type                  <- possible_types[type]
  
  if (type == "l"){
    xx = lt_id_l_q(xx)
    type = "q"
  }
  
  # a final catch
  out <- NULL
  # Abridged input lt
  if (is_abridged(Age)){
    
    # If we have nMx
    if (type == "m" & Single){
      
      args_could_have <- formals(lt_abridged2single)
      
      out <- lt_abridged2single(nMx = xx, Age = Age, Sex = Sex, ...)
    }
    if (type == "m" & !Single){
      out <- lt_abridged(nMx = xx, Age = Age, Sex = Sex, ...)  
    }
    # If we have nMx
    if (type == "q" & Single){
      out <- lt_abridged2single(nqx = xx, Age = Age, Sex = Sex,  ...)
    }
    if (type == "q" & !Single){
      out <- lt_abridged(nqx = xx, Age = Age, Sex = Sex,  ...)  
    }
  }
  
  if (is_single(Age)){
    if (type == "m" & Single){
      out <- lt_single_mx(nMx = xx, Age = Age, Sex = Sex,  ...)
    }
    if (type == "m" & !Single){
      out <- lt_single_mx(nMx = xx, Age = Age, Sex = Sex,  ...)
      out <- lt_single2abridged(lx = out$lx,nLx = out$nLx, ex = out$ex) 
    }
    if (type == "q" & Single){
      out <- lt_single_qx(nqx = xx, Age = Age, Sex = Sex,  ...)
    }
    if (type == "q" & !Single){
      out <- lt_single_qx(qx = xx, Age = Age, Sex = Sex,  ...)
      out <- lt_single2abridged(lx = out$lx,nLx = out$nLx, ex = out$ex) 
    }
  }
  
  if (is.null(out)){
    # a final catch
    stop("please check function arguments")
  }  
  return(out)
}
