
# Author: tim
###############################################################################


#' An oscillatory average of age splits.
#' @description Single ages can be grouped into 5-year age groups in 5 ways by staggering terminal digits. This method is a bit smoother than the standard Sprague or Beers methods, but not as smooth as \code{grabill()}.
#'
#' @details This function works on a single vector of single-age counts, not on a matrix. Results are not constrained to any particular age group, but are constrained to the total count. Negatives, \code{NA}, or \code{NaN} values are ignored in averaging. This can happen in older ages . It is recommended to run \code{monoCloseout()} or
#' similar after the oscillatory split in such situations.
#'
#' @param Value numeric. Vector of single age counts.
#' @param Age integer. Vector of single ages.
#' @param OAG logical. Whether or not the last value is the open age group. Default \code{TRUE}.
#' @param splitfun function used to split at each digit grouping. Default \code{sprague()}.
#' @param closeout logical or character. Default \code{"mono"}.
#' @param pivotAge integer. Age to start blending in closeout values.
#' @param ... optional arguments passed to \code{splitfun()}.
#'
#' @return Numeric vector of smoothed counts.
#' @references
#' \insertRef{booth2015demographic}{DemoTools}
#' @export
#' @examples
#' # code currently breaking, needs to be revisited and updates completed, sorry
#' \dontrun{
#'
#' Value     <- structure(pop1m_ind, .Names = 0:100)
#' #barplot(Value, main = "yup, these have heaping!")
#' # this is the basic case we compare with:
#' pop0      <- graduate_sprague(groupAges(Value), Age = 0:100,  OAG = TRUE)
#' # note: this function needs single ages to work because
#' # ages are grouped into 5-year age groups in 5 different ways.
#' # breaks
#' #pop1     <- splitOscillate(Value, OAG = TRUE, splitfun = graduate_sprague)
#' pop2      <- splitOscillate(Value, OAG = TRUE, splitfun = beers)
#' # what's smoother, splitOscillate() or graduate_grabill()?
#' # note, same closeout problem, can be handled by graduate_mono_closeout()
#' pop3      <- graduate_grabill(Value, OAG = TRUE)
#' # and technically you could give grabill as splitfun too
#' pop4      <- splitOscillate(Value, OAG = TRUE, splitfun = graduate_grabill)
#'
#' Age       <- 0:100
#' plot(Age, Value)
#' lines(Age, pop0, col = "blue")
#' # slightly smoother (also shifted though)
#' lines(Age, pop1)
#' # only different at very high ages, small nrs
#' lines(Age, pop2, col = "red", lty = 2, lwd = 2)
#' lines(Age, pop3, col = "magenta")
#' lines(Age, pop4, col = "orange", lty = 2)
#' legend("topright",
#' lty = c(1,1,2,1,2),
#' lwd = c(1,1,2,1,1),
#' col = c("blue","black","red","magenta","orange"),
#' 		legend = c("graduate_sprague()",
#'                 "splitOscillate(splitfun = graduate_sprague)",
#' 				   "splitOscillate(splitfun = beers)",
#' 				   "graduate_grabill()",
#'                 "splitOscillate(splitfun = graduate_grabill)"))
#'
#' # index of dissimilarity
#' ID(Value, pop0) # original vs sprague
#' ID(pop0,pop1) # sprague vs sprague osc
#' ID(pop1,pop2) # sprague osc vs beers osc
#' ID(pop2,pop3) # beers osc vs grabill
#' ID(pop3,pop4) # grabill vs grabill osc
#' # measure of smoothness:
#' mean(abs(diff(Value)))
#' mean(abs(diff(pop0)))
#' mean(abs(diff(pop1)))
#' mean(abs(diff(pop2)))
#' mean(abs(diff(pop3)))
#' mean(abs(diff(pop4)))
#' }
splitOscillate <- function(
  Value,
  Age = 1:length(Value) - 1,
  OAG = TRUE,
  splitfun = graduate_sprague,
  closeout = "mono",
  pivotAge = 90,
  ...) {
  N                              <- length(Value)
  if (OAG) {
    open                         <- Value[N]
    OA                           <- Age[N]
    Value                        <- Value[-N]
    Age                          <- Age[-N]
    N                            <- N - 1
  }
  TOT                            <- sum(Value)
  # select which ages to keep:
  p1x1                           <- matrix(nrow = N, ncol = 5)
  rownames(p1x1)                 <- Age
  for (i in 0:4) {
    # regroup ages
    Age.i.5                      <- calcAgeN(Age, shiftdown = i)
    # only use age groups w 5 single ages represented
    keep.i                       <-
      rep(rle(Age.i.5)$leng, rle(Age.i.5)$leng) == 5
    # cut vector down to those cases
    Age.i.5                      <- Age.i.5[keep.i]
    # cut counts down to those cases
    Val.i                        <- Value[keep.i]
    # group ages into said 5-year age groups
    Val.i.5                      <- groupAges(Val.i, AgeN = Age.i.5)
    
    # get first run estimate
    pop.est                      <- splitfun(Val.i.5, 
                                             Age = as.integer(names(Val.i.5)), 
                                             OAG = FALSE, ...)
    #        a                   <- rownames(pop.est)
    #		if (closeout){
    #			a.fake                   <- (1:nrow(pop.est) - 1) * 5
    #			pop.est                  <- monoCloseout(Val.i.5, Age = a.fake, pops = pop.est, OAG = FALSE)
    #		}
    
    pop.est[pop.est < 0]         <- 0
    p1x1[keep.i, i + 1]          <- pop.est
  }
  # take average per age
  p.out                          <- rowMeans(p1x1, na.rm = TRUE)
  # rescale toa proper total
  p.out                          <- rescale_vector(p.out, TOT)
  # re-append the open age group if needed
  if (OAG) {
    Age                          <- c(Age, OA)
    p.out                        <- c(p.out, open)
    names(p.out)                 <- Age
  }
  if (is.logical(closeout)) {
    if (!closeout) {
      return(p.out)
    }
    closeout                     <- "mono"
  }
  if (closeout == "mono") {
    p.out                        <-
      graduate_mono_closeout(
        popmat = Value,
        pops = p.out,
        OAG = OAG,
        pivotAge = 90
      )
  }
  
  p.out
}
