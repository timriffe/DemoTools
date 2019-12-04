

#' Graduate age groups using a monotonic spline.
#' @description Take the cumulative sum of \code{Value} and then run a monotonic spline through it. The first differences split back single-age estimates of \code{Value}. Optionally keep the open age group untouched.
#'
#' @details The \code{"monoH.FC"} method of \code{stats::splinefun()} is used to fit the spline because 1) it passes exactly through the points, 2) it is monotonic and therefore guarantees positive counts, and 3) it seems to be a bit less wiggly (lower average first differences of split counts). Single-age data is returned as-is. If you want to use this function as a smoother you first need to group to non-single ages.
#' @inheritParams graduate
#' @return Numeric. vector of single smoothed age counts.
#' @importFrom stats splinefun
#' @references
#' \insertRef{fritsch1980monotone}{DemoTools}
#' @export
#' @examples
#' Value                <- structure(c(88623, 90842, 93439, 96325, 99281, 102051, 104351,
#'				 106555, 109170, 112188, 113582, 112614, 108904, 102622, 95867,
#'				 80874, 60196, 37523, 17927, 5642, 1110), .Names = c("0", "5",
#'				 "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60",
#'				 "65", "70", "75", "80", "85", "90", "95", "100"))
#'
#' # overwrite open age group with a single age estimate for that age
#' # (doesn't extrapolate)
#' graduate_mono(Value)
#' # or respect open age group
#' graduate_mono(Value, OAG = TRUE)
#'
#' # Also accepts single ages:
#' Value                <- structure(pop1m_ind, .Names = 0:100)
#'
#' 		\dontrun{
#' 	ages                <- seq(0,100,5)
#'  plot(graduate_mono(Value),xlab = 'Age', ylab = 'Counts', type = 'l',main = 'Ungrouped counts')
#' 		}

graduate_mono   <- function(
  Value, 
  AgeInt, 
  Age, 
  OAG = TRUE) {
  
  if (as.character(match.call()[[1]]) == "splitMono") {
    warning("please use graduate_mono() instead of splitMono().", call. = FALSE)
  }
  
  if (missing(Age) & missing(AgeInt)) {
    Age                 <- names2age(Value)
  }
  if (missing(AgeInt)) {
    # give 1 to final interval to preserve
    AgeInt              <- age2int(Age, OAG = OAG, OAvalue = 1)
  }
  if (missing(Age)) {
    Age                 <- int2age(AgeInt)
  }
  
  # if age is single return as-is
  if (is_single(Age)) {
    return(Value)
  }
  
  # if last age Open, we preserve it
  if (OAG) {
    N                   <- length(Value)
    OAvalue             <- Value[N]
    Value               <- Value[-N]
    Age                 <- Age[-N]
    AgeInt              <- AgeInt[-N]
  }
  # if the final age is Open, then we should remove it and then
  # stick it back on
  
  AgePred               <- c(min(Age), cumsum(AgeInt))
  y                     <- c(0, cumsum(Value))
  AgeS                  <- min(Age):sum(AgeInt)
  y1                    <- splinefun(y ~ AgePred, method = "monoH.FC")(AgeS)
  out                   <- diff(y1)
  names(out)            <- AgeS[-length(AgeS)]
  
  # The open age group is maintained as-is.
  if (OAG) {
    out                 <- c(out, OAvalue)
    names(out)          <- AgeS
  }
  
  out
}

#' @export
#' @rdname graduate_mono
splitMono <- graduate_mono

#' blend the Sprague upper boundary age estimates into monotonic spline estimates
#'
#' @description A simple monotonic spline on the cumulative sum of population counts may return more convincing single age count estimates than the Sprague or other splitting methods. This function blends the given single age population estimates starting at \code{pivotAge}.
#'
#' @inheritParams graduate
#' @param pops optional numeric vector of single age population counts derived from \code{Value}.
#' @param pivotAge integer (default 90). Age at which to switch to spline-based estimates.
#' @param splitfun optional. The function used to create pops. Default \code{graduate_sprague}. Could also be \code{graduate_grabill}, or any other function that similarly transforms.
#' @param OAG logical (default \code{FALSE}). Would we like to re-impute the last
#' element of \code{Value} as the open age group?
#' @param ... arguments to be optionally passed to \code{splitfun()}.
#' @return numeric matrix of age by year estimates of single-age counts.
#'
#' @details The \code{pivotAge} must be at least 10 years below the maximum age detected from
#' \code{rownames(popmat)}, but not lower than 75. In the exact \code{pivotAge}, we may either take the Sprague estimates or the spline estimates, depending on which is larger, then the single-age estimates for this 5-year age group are rescaled to sum to the original total in \code{Value}. Higher ages are taken from the spline-based age splits. The spline results are derive from the \code{"monoH.FC"} method of \code{splinefun()} on the cumulative sum of the original age grouped data. One could use this function to perform the same closeout to Grabill estimates, if these are given via the \code{pops} argument. See examples. Note that the Grabill split method mixed with this closeout will not necessarily preserve the annual totals, and this function performs to rescaling. The open age group is preserved (and must be included in \code{Value}).
#'
#' @export
#'
#' @examples
#' a5 <- as.integer(rownames(pop5_mat))
#' popvec               <- pop5_mat[,1]
#' closed.out           <- graduate_mono_closeout(Value = popvec, Age = a5, OAG = TRUE)
#' sum(closed.out) - sum(popvec)
#' graduate_mono_closeout(Value = popvec, pivotAge = 85, Age = a5, OAG = TRUE)
#' # giving a different single-age split to close out this way:
#' popg                 <- graduate_grabill(Value = popvec, Age = a5, OAG = TRUE)
#' grabill.closed.out   <- graduate_mono_closeout(Value = popvec, Age = a5, pops = popg)
#' # totals not necessarily preserved if mixed w Grabill
#' # I wouldn't recommend a rescale of the total, since the
#' # only part we mess with here is the old age section. Ergo,
#' # one may wish to instead rescale results colSums() of
#' # popg at age pivotAge and higher.
#' sum(grabill.closed.out) - sum(popvec)
#' # also works on an age-labelled vector of data

#' closed.vec           <- graduate_mono_closeout(popvec, Age = a5, OAG = TRUE)
#' # let's compare this one with sprague()
#' simple.vec           <- graduate_sprague(popvec, Age = a5, OAG = TRUE)
#' # and with a simple monotonic spline
#' mono.vec             <- graduate_mono(popvec, Age = a5, OAG = TRUE)
#' \dontrun{
#' plot(85:100,simple.vec[86:101], type = 'l',
#'  main = "In this case graduate_sprague() is the smoothest")
#' lines(85:100,closed.vec[86:101], col = "red", lwd = 2)
#' lines(85:100,mono.vec[86:101], col = "blue", lty = 2)
#' legend("topright",lty=c(1,2,2), col = c("black","red","blue"),lwd = c(1,2,1),
#' 		legend = c("graduate_sprague()","monoCloseout()", "graduate_mono()"))
#' }

graduate_mono_closeout <-
  function(Value,
           Age,
           pops,
           pivotAge = 90,
           splitfun = graduate_sprague,
           OAG = TRUE,
           ...) {
    if (as.character(match.call()[[1]]) == "monoCloseout") {
      warning("please use graduate_mono_closeout() instead of monoCloseout().", call. = FALSE)
    }
    names(Value)        <- Age
    
    if (missing(pops)) {
      pops              <- splitfun(Value, Age = Age, OAG = OAG, ...)
    }
    # get the spline population split
    AgeIn               <- Age
    # this does not smooth single ages, it only splits to single
    popmono             <- graduate_mono(Value, OAG = OAG, Age = AgeIn)

    # some age pars
    Age                 <- as.integer(names(popmono))
    
    # some checks on pivotAge...
    if (!(max(Age) - 10) >= pivotAge) {
      pivotAge          <- max(Age) - 10
      if (pivotAge < 75) {
        warning(
          "pivotAge wasn't in rownames(popmat), moved it to 3rd ",
          "from bottom row of popmat, but appears to be < 75 ",
          "so returning sprague() output as-is, no extra closeout performed."
        )
        return(pops)
      }
      warning("pivotAge moved to ", pivotAge, ", continued.")
    }
    # -----------------------------
    # now begin the closeout blend.
    p.i                 <- which(Age == pivotAge)
    ## substitute Sprague interpolation if > pchip for better blending of the two series
    pop.c               <- popmono[p.i:(p.i + 4)]
    ind                 <- pops[p.i] > pop.c[1]
    if (ind){
      pop.c[1]       <- pops[p.i]
    }
   
    
    ## adjust back on initial pop 5x5 for age 90-94
    ## proportional distribution
    pop.c[is.na(pop.c)] <- 0
    prop                <- rescale.vector(pop.c, scale = 1)
    pivot5              <- Value[which(AgeIn == pivotAge)]
    pop.c               <- prop * pivot5
    ## append the remaining of the age groups (except last open age)
    ## 95-99 onward
    m                   <- length(pops)
    pop.c               <- c(pop.c, popmono[(p.i + 5):m])
    ## append Sprague interpolation before age 90
    pop.c               <- c(pops[1:(p.i - 1)], pop.c)
    
    ## deal with negative values if applicable (but in principle should not be happening)
    pop.c[pop.c < 0]    <- 0
    
    # label and return
    #dimnames(pop.c)    <- list(Age1, colnames(popmat))
    names(pop.c)     <- Age

    pop.c
  }

#' @export
#' @rdname graduate_mono_closeout
monoCloseout <- graduate_mono_closeout
