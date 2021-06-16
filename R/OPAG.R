# TODO / wish-list
# [ ] All choice of Age_fit / AgeInt_fit based on criteria (see PJ Drive folder)
# [ ] ensure age groups are flexible as required.
#    [ ] what happens when one input is in a different age group than another?
# [x] OPAG_simple() should allow for non-single ages
# [ ] add unit tests
#    [ ] check for age group commensurability and warn if necessary
# [ ] add more examples to OPAG?
# [ ] remove rownames message from DownloadLx(), haha I DON'T SEE THIS
# [ ] test OPAG_simple() with non-single ages groups, update documentation if necessary.
# [ ] harmonize args betwenn OPAG_simple and OPAG family.

# [x] make AgeInt not required
# [ ] document Age_fit better
# [x] change open age formula in warp function
# [x] fix continuous = FALSE formula/ then removed 
# [x] ensure Lx is single ages once for use throughout 
# [x] change default method to "mono" 

# Author: tim
###############################################################################
# distribute population in open age group over higher ages.
# The PAS implementation uses stable populations, and it will be added
# here in the future, as well as other options. The main missing piece
# is a good collection of model lifetables.

#' redistripute an open age group count over higher ages proportional to an arbitrary standard
#' @description This method could be useful whenever a reasonable standard is available. At present the standard must be supplied by the user.
#' @details In this implementation both the original population counts and the standard must be in single ages.
#' @param Pop numeric vector of population counts
#' @param Age integer vector of single age lower bounds
#' @param OAnow integer. The lower age bound above which counts will be redistributed
#' @param StPop numeric vector of standard population counts
#' @param StAge integer vector of single age lower bounds for the standard population
#' @param OAnew integer. The desired new open age, must be no higher than \code{max(StAge)}.
#' @export
#' @references
#' \insertRef{PAS}{DemoTools}
#' @examples
#'  Pop        <- c(38129,38382,38824,39275,39500,37304,35152,
#'  34061,33911,32875,31599,30376,29822,29691,28765,
#'  28695,28917,28203,29209,30316,29062,26977,26577,
#'  27727,28599,30513,31774,32347,34093,33736,32085,
#'  30807,28279,26873,25612,23503,22207,21388,20122,
#'  18014,15626,15006,14158,11195,7931,7640,9053,
#'  13276,17226,18918,17697,18424,17723,16706,14410,
#'  13342,14787,15183,15727,16045,14777,14267,13102,
#'  10866,9311,6933,5030,3785,3551,2848,3080,
#'  2874,2368,2681,3165,3010,3009,2721,2705,
#'  2492,2244,1971,1644,1565,1307,5027)
#'  Age        <- 0:85
#'  # standard pop taken from ages 55+
#'  StPop      <- c(6258,6177,6089,5995,5894,5787,5672,5552,
#'  5423,5286,5140,4985,4824,4652,4477,4293,
#'  4107,3912,3712,3502,3282,3055,2823,2591,
#'  2360,2138,1921,1710,1502,1297,1098,910,
#'  741,592,463,353,265,192,137,95,
#'  63,42,26,17,9,13)
#'
#'  StAge      <- 55:100
#'
#' PopExtended <- OPAG_simple(
#'  		Pop = Pop,
#'  		Age = Age,
#'  		StPop = StPop,
#'  		StAge = StAge)
#'
#' \dontrun{
#'  plot(Age, Pop, type = 'l',xlim=c(80,100),ylim=c(0,1e4))
#' lines(0:100, PopExtended, col = "red", lty = 2)
#' }
#' stopifnot((sum(PopExtended[86:101]) - Pop[86]) == 0)

# TR: needs to be generalized to age group widths...
OPAG_simple    <-
  function(Pop,
           Age,
           OAnow = max(Age),
           StPop,
           StAge,
           OAnew = max(StAge)) {
    # # assume single NOT NEEDED See age concordance
    # stopifnot(is_single(Age))
    # stopifnot(is_single(StAge))
    # OAG can be less than or equal to max age
    stopifnot(OAnow %in% Age)
    stopifnot(OAnew %in% StAge)
    # age and pop vectors must match lengths, assume ordered
    stopifnot(length(Pop) == length(Age))
    stopifnot(length(StPop) == length(StAge))
    # age concordance
    minStAge = min(StAge)
    stopifnot(all(Age[Age >= minStAge] %in% StAge))

    # group pop down to OAG
    Pop        <- groupOAG(Pop, Age, OAnow)
    StPop      <- groupOAG(StPop, StAge, OAnew)
    
    # even up lengths
    N          <- length(Pop)
    Age        <- Age[1:N]
    OAtot      <- Pop[N]
    # same for standard
    StN        <- length(StPop)
    StAge      <- StAge[1:StN]
    
    # make standard distribution.
    standard   <- rescale_vector(StPop[StAge >= OAnow], scale = 1)
    # redistribute OAG
    PopUpper   <- OAtot * standard
    # keep lower ages of Pop
    PopLower   <- Pop[1:(N - 1)]
    
    # graft, name, and return
    out        <- c(PopLower, PopUpper)
    Ageout     <- sort(unique(c(Age, StAge)))
    names(out) <- Ageout
    
    out
  }


#' Warps a given stationary population into a stable population
#' @description We take `nLx` as indicative of a stationary population age structure,
#' then subject the population structure to long-term growth by a constant rate, `r`.
#' @details `Lx1` could be any population structure of any scale, as long as you're comfortable
#' assuming it's stationary and can be warped into stable. For the oldest ages, this is probably
#' quite often an acceptable and useful approximation. The transformation is applied at the single-age scale, even if the input `nLx` is in wider (e.g. abridged) age groups. When needed, we reduce to single ages using (default) `graduate_uniform()`, then apply the transformation, then group back. This is innocuous if `nLx` is given in single ages. You may want to change `method` to `"mono"` or `"pclm"`.
#'
#' @param Lx1 numeric vector of stationary population age structure in arbitrary integer age groups
#' @param Age_Lx1 interger vector of lower bounds of age groups of `nLx`
#' @param r stable growth rate
#' @return numeric vector of the transformed `nLx`. Note, this vector sums to `1`.
#' @export
#' @examples
#' Lx  <- downloadnLx(NULL, "Spain","female",1971)
#' Age <- names2age(Lx)
#' r   <- .01
#' ai  <- age2int(Age)
#' \dontrun{
#' plot(Age, Lx/sum(Lx) / ai, type = 's',ylim = c(0,0.015))
#' lines(Age,OPAG_nLx_warp_r(Lx,Age,0.005)/ai,type='s',col = "red")
#' lines(Age,OPAG_nLx_warp_r(Lx,Age,-0.005)/ai,type='s',col = "blue")
#' }
#'

OPAG_nLx_warp_r <- function(Lx1,
                            Age_Lx1,
                            r
){
  a1   <- Age_Lx1
  w1Lx  <- exp(-r * (a1 + .5)) * Lx1
  nAges <- length(Age_Lx1)
  
  if (r == 0 ) {
    w1Lx[nAges]  <- Lx1[nAges]  
  } else {
    Tlast <- Lx1[nAges]
    Tprev <- Tlast + Lx1[nAges-1]
    # PJ's closeout
    abar <- -((log(Lx1[nAges-1])-r*a1[nAges-1])-log(Tprev*exp(r)-Tlast)) / r
    w1Lx[nAges]  <- exp(-r * abar) * Tlast
  }
  
  w1Lx  <- w1Lx / sum(w1Lx)
  
  w1Lx
}

#' calculates residual for optimizing growth rate r for OPAG family
#' @description For a given set of age groups to fit against, and a given stable growth rate, $r$,
#' what is the error implied given the current $r$ and stationary standard?
#' @details This is a utility function for `OPAG()`, which needs to optimize $r$ for a
#' given population vector and stationary standard.
#' @param r given stable growth rate
#' @param Pop_fit numeric vector of at least two population counts to use for fitting
#' @param Age_fit integer vector of lower bounds for age groups of `Pop_fit`
#' @param AgeInt_fit integer vector of widths of age groups of `Pop_fit`
#' @param Lx1 numeric vector of stable population standard by single ages
#' @param Age_Lx1  integer vector of lower bounds for age groups of `Lx1`
#' @return numeric. A residual that you're presumably trying to minimize.
#' @export

#' @examples
#' # Make up some population data to fit to:
#' Pop_fit    <- c(85000,37000)
#' Age_fit    <- c(70,80)
#' AgeInt_fit <- c(10,10)
#' nLx        <- downloadnLx(NULL, "Spain","female",1971)
#' # graduate(nLx, Age_nLx, method = method, constrain = TRUE)
#' Ageab    <- names2age(nLx)
#' Lx1 <- graduate(c(nLx), Ageab, method = "mono", constrain = TRUE)
#' Age_Lx1 <- 0:100
#' r          <- .01
#'
#' OPAG_r_min(r,
#'            Pop_fit = Pop_fit,
#'            Age_fit = Age_fit,
#'            AgeInt_fit = AgeInt_fit,
#'            Lx1 = Lx1,
#'            Age_Lx1 = Age_Lx1)
#'
#' (r_opt <- optimize(OPAG_r_min,
#'          Pop_fit = Pop_fit,
#'          Age_fit = Age_fit,
#'          AgeInt_fit = AgeInt_fit,
#'          Lx1 = Lx1,
#'          Age_Lx1 = Age_Lx1,
#'          interval = c(-0.05,.05))$min)
#'

OPAG_r_min <- function(r,
                       Age_fit,
                       Pop_fit,
                       AgeInt_fit, # necessary
                       Lx1,
                       Age_Lx1
){
  AgeInt_nLx   <- age2int(Age_Lx1, OAvalue = 1)
  
  # 1) This is the standard we want to match to Pop,
  # which has presumably been cut down / grouped to the
  # ages we want to calibrate to.
  w1Lx    <- OPAG_nLx_warp_r(
    Lx1 = Lx1,
    Age_Lx1 = Age_Lx1,
    r = r
  )
  
  # 2) now need to get it to the same age groups as Pop
  # so that we can get a residual
  
  w1Lx_fit <- rep(NA, length(Age_fit))
  
  for (i in 1:length(Age_fit)){
    ind <- Age_Lx1 >= Age_fit[i] & Age_Lx1 < (Age_fit[i] + AgeInt_fit[i])
    w1Lx_fit[i] <- sum(w1Lx[ind])
  }
  
  # 3) rescale standard and Pop_fit to sum to 1  
  stand   <- rescale_vector(w1Lx_fit, scale = 1)
  Pop_fit <- rescale_vector(Pop_fit, scale = 1)
  
  # 4) return the residual
  sum(abs(stand - Pop_fit))
}


#' creates stable standard based on optimizing the growth rate
#' @description The stationary standard, `nLx` is transformed into a stable standard by optimizing a growth rate, `r` such that the stable standard matches observed population counts in selected age groups. Usually the ages used for fitting are wide age groups in older ages preceding the open age group. The standard output by this function is used by `OPAG` to create the standard used to redistribute counts over older age groups up to a specified open age group, such as 100.
#' @details The argument `method` don't have much leverage on the result. In short, the stable population transformation is done by ungrouping `nLx` to single ages (if it isn't already), and `method` controls which graduation method is used for this, where `"uniform"`, `"mono"`, `"pclm"` are the reasonable choices at this writing. 
#' 
#' 
#' @inheritParams OPAG_r_min
#' @return
#' list constaining
#' 1.  `Standard` numeric vector, the transformed `nLx` to be used for
#'     redistribution in `OPAG()`
#' 2. r_opt the output of `optimize()`, where `min` is the growth parameter, `r`
#' @export
#' @importFrom stats optimize
#'
#' @examples
#' Pop_fit    <- c(85000,37000)
#' Age_fit    <- c(70,80)
#' nLx        <- downloadnLx(NULL, "Spain","female",1971)
#' Age_nLx    <- names2age(nLx)
#' Lx1 <- graduate(nLx,Age=Age_nLx,method = "mono")
#' Age_Lx1 <- 0:100
#' # India Males, 1971
#' Pop        <- smooth_age_5(pop1m_ind,
#'                            Age = 0:100,
#'                            method = "Arriaga")
#' Pop80      <- groupOAG(Pop, names2age(Pop), 80)
#' Age        <- names2age(Pop80)
#' 
#' nLx        <- downloadnLx(NULL, "India","male",1971)
#' Age_nLx    <- names2age(nLx)
#' 
#' 
#' Pop_fit    <- groupAges(Pop80, Age, N = 10)[c("60","70")]
#' Age_fit    <- c(60,70)
#' AgeInt_fit <- c(10,10)
#' 
#' Standard <- OPAG_fit_stable_standard(
#'   Pop_fit = Pop_fit,
#'   Age_fit = Age_fit,
#'   AgeInt_fit = AgeInt_fit,
#'   Lx1=Lx1,
#'   Age_Lx1 = Age_Lx1
#' )
#' 
#' # A visual comparison:
#' nL60 <- rescale_vector(nLx[Age_nLx >= 60])
#' St60p <- rescale_vector( Standard$Standard[c(0:100) >= 60] )
#' ages_plot <- seq(60,100,by=5)
#' \dontrun{
#'   plot(ages_plot,nL60, type = 'l')
#'   lines(60:100, St60p, col = "blue")
#' }

OPAG_fit_stable_standard <- function(Pop_fit,
                                     Age_fit,
                                     AgeInt_fit,
                                     Lx1,
                                     Age_Lx1
){
  
  
  # optimize the parameter r
  r_opt <- optimize(OPAG_r_min,
                    Pop_fit = Pop_fit,
                    Age_fit = Age_fit,
                    AgeInt_fit = AgeInt_fit,
                    Lx1 = Lx1,
                    Age_Lx1 = Age_Lx1, 
                    interval = c(-0.02, .05)) # changed interval
  
  
  standard <- OPAG_nLx_warp_r(Lx1 = Lx1,
                              Age_Lx1 = Age_Lx1,
                              r = r_opt$min
  )
  # return both stable standard and the optimization output,
  # which will let us know if r is simply unreasonable or similar.
  out <- list(Standard = standard,
              r_opt = r_opt)
  out
}

#' Redistribute population over a specified age based on a stable standard fit to the data
#' @description This can be used as an external check of population counts
#' in older ages, assuming the stable population standard is representative enough, or it can be used to redistribute population in ages above a
#' specified ages `Redistribute_from`. This is handy, for instance, for
#' ensuring all censuses extend to a specified maximum age (e.g. 100+)
#' prior to intercensal interpolations. The assumption is that, at least in
#'  ages including `Age_fit` and higher ages, the population should follow
#'  a stable pattern proportional to a given survival curve subject to
#'  constant growth, `r`.
#' @details It may be helpful to try more than one fitting possibility,
#' and more than one `Redistribute_from` cut point, as results may vary.
#'
#' `Redistribute_from` can be lower than your current open age group,
#' and `OAnew` can be higher, as long as it is within the range of `Age_nLx`.
#' If `Age_nLx` doesn't go high enough for your needs, you can extrapolate
#' it ahead of time. For this, you'd want the `nMx` the underly it, and you
#' can use `lt_abridged()`, specifying a higher open age, and then
#' extracting `nLx` again from it.
#'
#' @inheritParams OPAG_r_min
#' @param Pop numeric vector of population counts
#' @param Age_Pop integer vector of the lower bounds of the population age groups
#' @param nLx numeric vector of stationary population age structure in arbitrary integer age groups
#' @param Age_nLx interger vector of lower bounds of age groups of `nLx`
#' @param Redistribute_from integer lower age bound that forms the cutoff, above which we redistribute counts using the stable standard.
#' @param OAnew integer. Desired open age group in the output (must being element of `Age_nLx`)
#' @param method character, graduation method used for intermediate graduation. Default `"mono"`. Other reasonable choices include `"pclm"` or `"uniform"`.
#' @export
#' @importFrom utils tail
#' @examples
#' # India Males, 1971
#' Pop            <- smooth_age_5(pop1m_ind,
#'                          Age = 0:100,
#'                          method = "Arriaga")
#' Age_Pop        <- names2age(Pop)
#' AgeInt_Pop     <- age2int(Age_Pop, OAvalue = 1)
#'
#' nLx            <- downloadnLx(NULL, "India","male",1971)
#' Age_nLx        <- names2age(nLx)
#' AgeInt_nLx     <- age2int(Age_nLx, OAvalue = 1)
#'
#' Pop_fit <- OPAG(Pop,
#'     Age_Pop = Age_Pop,
#'     nLx = nLx,
#'     Age_nLx = Age_nLx,
#'     Age_fit =  c(60,70),
#'     AgeInt_fit = c(10,10),
#'     Redistribute_from = 80)
#'
#' \dontrun{
#' # look at 75+
#' ind <- Age_Pop >= 75
#' plot(Age_Pop[ind], Pop[ind])
#' lines(Age_Pop[ind], Pop_fit$Pop_out[ind], col = "blue")
#'
#' # relative differences in ages 80+
#' ind <- Age_Pop >= 80
#' plot(Age_Pop[ind],  (Pop_fit$Pop_out[ind] - Pop[ind]) / Pop[ind])
#'}

OPAG <- function(Pop,
                 Age_Pop,
                 nLx, 
                 Age_nLx, 
                 Age_fit = NULL,
                 AgeInt_fit = NULL,
                 Redistribute_from = max(Age_Pop),
                 OAnew = max(Age_nLx),
                 method = "mono" 
){
  
  # ensure OAnew is possible
  stopifnot(OAnew <= max(Age_nLx))
  
  # TB: if OAnew < min(Age_nLx) that's an error
  
  method <- match.arg(method, choices = c("uniform","pclm","mono"))
  
  #TB: checking if pop and nLx have different intervals and warning users - still working on it
  if(!identical(as.integer(unique(diff(Age_Pop))), as.integer(unique(diff(Age_nLx))))){ # put a different
    cat("\nAge_Pop and Age_nLx age intervals are different!\n")
  }
  
  # PJ adds this. Note final age group not assigned a width
  AgeInt_Pop <- diff(Age_Pop)
  # AgeInt_nLx <- diff(Age_Pop)
  
  # setup, prelims:
  # 0) if Age_fit isn't given assume last two 10-year age groups.
  
  if (is.null(Age_fit)){
    OA         <- max(Age_Pop)
    Age_fit    <- OA - c(20,10)
    AgeInt_fit <- c(10,10)
    stopifnot(Age_fit %in% Age_Pop)
  }
  if (is.null(AgeInt_fit)){
    # assume age intervals are age differences, and repeat last one
    AgeInt_fit <- diff(Age_fit)
    AgeInt_fit <- c(AgeInt_fit, tail(AgeInt_fit, n=1))
    # if Age_fit includes pop OA then set last fit age int to Inf
    if (tail(Age_fit,1) == tail(Age_Pop,1)) {
      AgeInt_fit[length(AgeInt_fit)] <- Inf
    }
  }
  if (any(!Age_fit %in% Age_Pop)){
    ind <- Age_fit %in% Age_Pop
    Age_fit <- Age_fit[ind]
    AgeInt_fit <- AgeInt_fit[ind]
    stopifnot(length(Age_fit) > 1)
  }
  
  # 1) get Pop_fit
  
  # TR: note: this calls for a special age utility function I think
  # earmarking this code chunk to swap it out in the future.
  Pop_fit <- rep(NA, length(Age_fit))
  for (i in 1:length(Age_fit)){
    ind <- Age_Pop >= Age_fit[i] & Age_Pop < (Age_fit[i] + AgeInt_fit[i])
    Pop_fit[i] <- sum(Pop[ind])
  }
  
  # 2) make sure Lx is single ages
  Lx1  <- graduate(nLx, Age_nLx, method = method, constrain = TRUE)
  Age_Lx1 <- as.integer(names(Lx1))
  
  Stab_stand <- OPAG_fit_stable_standard(Pop_fit,
                                         Age_fit,
                                         AgeInt_fit,
                                         Lx1,
                                         Age_Lx1
  )
  StPop <- Stab_stand$Standard
  
  # 3) get total to redistribute:
  OAG_total <- sum(Pop[Age_Pop >= Redistribute_from])
  
  # 4) select standard in those age groups.
  StPop_sel <- StPop[Age_Lx1 >= Redistribute_from]
  StPop_sel <- rescale_vector(StPop_sel, scale = 1)
  
  # 5) redistribute
  Pop_redistributed <- StPop_sel * OAG_total
  
  # 5a) regroup into original pop age grouping
  # TR: possibly not the bes tuse of AgeInt_Pop..
  if (tail(AgeInt_Pop, n=2)[-1] == 5) {
    Pop_redistributed <- groupAges(Pop_redistributed, N = 5)
  }
  
  # 6) graft together
  Pop_grafted <- c(Pop[Age_Pop < Redistribute_from],
                   Pop_redistributed)
  Age_grafted <- c(Age_Pop[Age_Pop < Redistribute_from],
                   Age_nLx[Age_nLx >= Redistribute_from])
  
  names(Pop_grafted) <- Age_grafted
  # 7) potentially group down OAG
  Pop_out <- groupOAG(Value = Pop_grafted,
                      Age = Age_grafted,
                      OAnew = OAnew)
  Age_out <- names2age(Pop_out)
  
  # 8) compose list for output
  out <- list(
    Pop_out = Pop_out,
    Age_out = Age_out,
    Pop_in = Pop,
    Standard = StPop,
    r_opt = Stab_stand$r_opt)
  
  out
}
