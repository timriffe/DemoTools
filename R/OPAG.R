# TODO / wish-list
# [ ] All choice of Age_fit / AgeInt_fit based on criteria (see PJ Drive folder)
# [ ] ensure age groups are flexible as required. 
#    [ ] what happens when one input is in a different age group than another?
# [ ] OPAG_simple() should allow for non-single ages
# [ ] add unit tests
# [ ] add more examples to OPAG?
# [ ] remove rownames message from DownloadLx(), haha

# Author: tim
###############################################################################
# distribute population in open age group over higher ages.
# The PAS implementation uses stable populations, and it will be added 
# here in the future, as well as other optiond. The main missing piece 
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
    # assume single
    stopifnot(is_single(Age))
    stopifnot(is_single(StAge))
    # OAG can be less than or equal to max age
    stopifnot(OAnow %in% Age)
    # age and pop vectors must match lengths, assume ordered
    stopifnot(length(Pop) == length(Age))
    # age concordance
    #stopifnot(all(Age %in% StAge))
    
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
    
    # make stadnard distribution.
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
#' @details `nLx` could be any population structure of any scale, as long as you're comfortable 
#' assuming it's stationary and can be warped into stable. For the oldest ages, this is probably 
#' quite often an acceptable and useful approximation. The transformation is applied at the single-age scale, even if the input `nLx` is in wider (e.g. abridged) age groups. When needed, we reduce to single ages using (default) `graduate_uniform()`, then apply the transformation, then group back. This is innocuous if `nLx` is given in single ages. You may want to change `method` to `"mono"` or `"pclm"`.
#' 
#' @param nLx numeric vector of stationary population age structure in arbitrary integer age groups
#' @param Age interger vector of lower bounds of age groups of `nLx`
#' @param r stable growth rate
#' @param AgeInt optional integer vector of widths of age groups, inferred if not given.
#' @param continous logical. If `TRUE` we use the growth adjustment. `e^(-age*r)`. If `FALSE` we assume `r` is geometric growth, and we use `(1+r)^age` for the growth adjustment.
#' @param method character, graduation method used for intermediate graduation. Default `"uniform"`. Other reasonable choices include `"mono"` or `"pclm"`.
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

OPAG_nLx_warp_r <- function(nLx, 
                            Age, 
                            r, 
                            AgeInt = NULL, 
                            continuous = TRUE, 
                            method = "uniform"){
  # Let's do this in single ages :-)
  # for now, just uniform, but could pass in args to graduate of course
  # if that is preferred.
  Lx1  <- graduate(nLx, Age, method = method, constrain = TRUE)
  a1   <- names2age(Lx1) 
  if (continuous){
     wLx  <- exp(-r * (a1 + .5)) * Lx1
  } else {
    # then geometric
    w     <- (1 + r) ^ (a1 + .5)
    wLx   <- w * nLx
  }
  wLx  <- wLx / sum(wLx)
  if (is.null(AgeInt)){
    AgeInt   <- age2int(Age, OAvalue = 1)
  }
  a12A <- rep(Age, AgeInt)
  nwLx <- groupAges(wLx, Age = a1, AgeN = a12A)
  nwLx
}

#' calculates residual for optimizing growth rate r for OPAG family
#' @description For a given set of age groups to fit against, and a given stable growth rate, $r$,
#' what is the error implied given the current $r$ and stationary standard?
#' @details This is a utiltiy function for `OPAG()`, which needs to optimize $r$ for a 
#' given population vector and stationary standard.
#' @param r given stable growth rate
#' @param Pop_fit numeric vector of at least two population counts to use for fitting
#' @param Age_fit integer vector of lower bounds for age groups of `Pop_fit`
#' @param AgeInt_fit integer vector of widths of age groups of `Pop_fit`
#' @param nLx numeric vector of stable population standard
#' @param Age_nLx  integer vector of lower bounds for age groups of `nLx`
#' @param AgeInt_nLx optional integer vector of widths of age groups of `nLx`, inferred if not given. 
#' @param continous logical. If `TRUE` we use the growth adjustment. `e^(-age*r)`. If `FALSE` we assume `r` is geometric growth, and we use `(1+r)^age` for the growth adjustment.
#' @param method character. Graduation method, default `"uniform"`. `"mono"` or `"pclm"` would also be good choices.
#' @return numeric. A residual that you're presumably trying to minimize.
#' @export

#' @examples
#' Make up some population data to fit to:
#' Pop_fit    <- c(85000,37000)
#' Age_fit    <- c(70,80)
#' AgeInt_fit <- c(10,10)
#' nLx        <- downloadnLx(NULL, "Spain","female",1971)
#' Age_nLx    <- names2age(nLx)
#' r          <- .01
#' 
#' OPAG_r_min(r, 
#'            Pop_fit, 
#'            Age_fit, 
#'            AgeInt_fit,
#'            nLx, 
#'            Age_nLx)
#' 
#' (r_opt <- optimize(OPAG_r_min, 
#'          Pop_fit = Pop_fit,
#'          Age_fit = Age_fit,
#'          AgeInt_fit = AgeInt_fit,
#'          nLx = nLx,
#'          Age_nLx = Age_nLx,
#'          interval = c(-0.05,.05))$min)
#' 
#' ai  <- age2int(Age_nLx)
#' 
#' # Note the whole age range is being scaled to 1 here, but in practice
#' # you'd only be doing this in the highest ages. If only two fitting
#' # ages are given, then we can get an r that matches them perfectly,
#' # at least within reasonable bounds.
#' \dontrun{
#' plot(Age, nLx / sum(nLx) / ai, type = 's')
#' lines(Age,OPAG_nLx_warp_r(Lx,Age,r=r_opt)/ai,type='s',col = "red")
#' }

OPAG_r_min <- function(r, 
                       Pop_fit, 
                       Age_fit, 
                       AgeInt_fit, # necessary
                       nLx, 
                       Age_nLx, 
                       AgeInt_nLx = NULL, 
                       continuous = TRUE,
                       method = "uniform"){
  if (is.null(AgeInt_nLx)){
    AgeInt_nLx   <- age2int(Age_nLx, OAvalue = 1)
  }
  # This is the standard we want to match to Pop,
  # which has presumably been cut down / grouped to the
  # ages we want to calibrate to.
  wnLx    <- OPAG_nLx_warp_r(
               nLx = nLx, 
               Age = Age_nLx, 
               r = r, 
               AgeInt = AgeInt_nLx, 
               continuous = continuous)
  
  # now need to get it to the same age groups as Pop 
  # so that we can get a residual
  
  # 1) Move stable pop to single ages
  w1Lx    <- graduate(
               wnLx, 
               Age = Age_nLx, 
               AgeInt = AgeInt_nLx,
               method = method)
  a1t     <- names2age(w1Lx)
  a1t     <- as.integer(a1t)
  
  # 2) which single ages implied by Pop?
  N       <- length(AgeInt_fit)
  a1match <- Age_fit[1]:(max(Age_fit) + AgeInt_fit[N] - 1)
  a1match <- as.integer(a1match)
  
  # 3) select down to just those ages:
  ind     <- a1t %in% a1match
  w1Lx    <- w1Lx[ind]
  
  # 4) group w1Lx to same as Pop_fit
  ageN    <- rep(Age_fit, times = AgeInt_fit)
  stand   <- groupAges(w1Lx, Age = a1match, AgeN = ageN)
  
  # 5) rescale standard and Pop_fit to sum to 1
  stand   <- rescale_vector(stand, scale = 1)
  Pop_fit <- rescale_vector(Pop_fit, scale = 1)
  
  # 6) return the residual
  sum(abs(stand - Pop_fit))
}


#' creates stable standard based on optimizing the growth rate
#' @description The stationary standard, `nLx` is transformed into a stable standard by optimizing a growth rate, `r` such that the stable standard matches observed population counts in selected age groups. Usually the ages used for fitting are wide age groups in older ages preceding the open age group. The standard output by this function is used by `OPAG` to creat the standard used to redistribute counts over older age groups up to a specified open age group, such as 100.
#' @details The arguments `method` and `continous` don't have much leverage on the result. In short, the stable population transformation is done by ungrouping `nLx` to single ages (if it isn't already), and `method` controls which graduation method is used for this, where `"uniform"`, `"mono"`, `"pclm"` are the reasonable choices at this writing. In single ages, the difference between using a geometric `r` versus continuous `r` are quite small for this task.
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
#' AgeInt_fit <- c(10,10)
#' nLx        <- downloadnLx(NULL, "Spain","female",1971)
#' Age_nLx    <- names2age(nLx)
#' 
#' # India Males, 1991
#' Pop        <- smooth_age_5(pop1m_ind,
#'                            Age = 0:100,
#'                            method = "Arriaga")
#' Pop80      <- groupOAG(Pop, names2age(Pop), 80)
#' Age        <- names2age(Pop80)
#' AgeInt     <- age2int(Age, OAvalue = 1) 
#' 
#' nLx        <- downloadnLx(NULL, "India","male",1991)
#' Age_nLx    <- names2age(nLx)
#' AgeInt_nLx <- age2int(Age_nLx,OAvalue = 1)
#' 
#' Pop_fit    <- groupAges(Pop80, Age, N = 10)[c("60","70")]
#' Age_fit    <- c(60,70)
#' AgeInt_fit <- c(10,10)
#'  
#' Standard <- OPAG_fit_stable_standard(
#'              Pop_fit,
#'              Age_fit,
#'              AgeInt_fit,
#'              nLx = nLx,
#'              Age_nLx = Age_nLx,
#'              AgeInt_nLx = AgeInt_nLx,
#'              method = "uniform",
#'              continuous = TRUE)
#' 
#' # A visual comparison:
#' nL60 <- rescale_vector(nLx[Age_nLx >= 60]) 
#' St60p <- rescale_vector( Standard$Standard[Age_nLx >= 60] )
#' ages_plot <- seq(60,100,by=5)
#' \dontrun{
#' plot(ages_plot,nL60, type = 'l')
#' lines(ages_plot, St60p, col = "blue")
#' }

OPAG_fit_stable_standard <- function(Pop_fit,
                                     Age_fit,
                                     AgeInt_fit,
                                     nLx,
                                     Age_nLx,
                                     AgeInt_nLx,
                                     method = "uniform",
                                     continuous = TRUE){


  # optimize the parameter r
  r_opt <- optimize(OPAG_r_min, 
                    Pop_fit = Pop_fit,
                    Age_fit = Age_fit,
                    AgeInt_fit = AgeInt_fit,
                    nLx = nLx,
                    Age_nLx = Age_nLx,
                    interval = c(-0.05, .05))
  

  standard <- OPAG_nLx_warp_r(nLx = nLx,
                              Age = Age_nLx,
                              r = r_opt$min,
                              AgeInt = AgeInt_nLx,
                              continuous = continuous,
                              method = method)
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
#' The argument `"method"` refers to which graduation method (see `?graduate`)
#' is only relevant if input data are in grouped ages. This is innocuous if 
#' ages are single to begin with. The choice of whether to assume 
#' `continuous = TRUE` constant growth versus geometric (`FALSE`) growth 
#' has little leverage. 
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
#' @param AgeInt_Pop integer vector of the population age group interval widths, using `Inf` for the open age group.

#'India Males, 1991
#' Pop            <- smooth_age_5(pop1m_ind,
#'                          Age = 0:100,
#'                          method = "Arriaga")
#' Age_Pop        <- names2age(Pop)
#' AgeInt_Pop     <- age2int(Age_Pop, OAvalue = 1) 
#' 
#' nLx            <- downloadnLx(NULL, "India","male",1991)
#' Age_nLx        <- names2age(nLx)
#' AgeInt_nLx     <- age2int(Age_nLx, OAvalue = 1)
#' 
#' Pop_fit <- OPAG(Pop, 
#'     Age_Pop = Age_Pop, 
#'     AgeInt_Pop = AgeInt_Pop,
#'     nLx = nLx,
#'     Age_nLx = Age_nLx,
#'     AgeInt_nLx,
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
                 AgeInt_Pop,
                 nLx, 
                 Age_nLx, 
                 AgeInt_nLx = NULL,
                 Age_fit = NULL,
                 AgeInt_fit = NULL,
                 Redistribute_from = max(Age_Pop),
                 OAnew = max(Age_nLx),
                 method = "uniform", 
                 continuous = TRUE){
  
  # setup, prelims:
  # 0) if Age_fit isn't given assume last two 10-year age groups.
  if (is.null(Age_fit)){
    OA         <- max(Age_Pop)
    Age_fit    <- OA - c(20,10)
    AgeInt_fit <- c(10,10)
    stopifnot(Age_fit %in% Age_Pop)
  }
  
  # 1) get Pop_fit
  
  # TR: note: this calls for a special age utility function I think
  # earmarking this code chunk to swap it out in the future.
  Pop_fit <- rep(NA, length(Age_fit))
  for (i in 1:length(Age_fit)){
    ind <- Age_Pop >= Age_fit[i] & Age_Pop < (Age_fit[i] + AgeInt_fit[i])
    Pop_fit[i] <- sum(Pop[ind])
  }
  
  # 2) get the standard
  Stab_stand <- OPAG_fit_stable_standard(Pop_fit,
                                         Age_fit,
                                         AgeInt_fit,
                                         nLx,
                                         Age_nLx,
                                         AgeInt_nLx,
                                         method = method,
                                         continuous = continuous)
  StPop <- Stab_stand$Standard
  
  # 3) get total to redistribute:
  OAG_total <- sum(Pop[Age_Pop >= Redistribute_from])
  
  # 4) select standard in those age groups.
  StPop_sel <- StPop[Age_nLx >= Redistribute_from]
  StPop_sel <- rescale_vector(StPop_sel, scale = 1)
  
  # 5) redistribute
  Pop_redistributed <- StPop_sel * OAG_total
  
  # 6) graft together
  Pop_grafted <- c(Pop[Age_Pop < Redistribute_from], 
                   Pop_redistributed)
  Age_grafted <- c(Age_Pop[Age_Pop < Redistribute_from],
                   Age_nLx[Age_nLx >= Redistribute_from])
  
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


