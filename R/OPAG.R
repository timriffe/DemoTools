# Author: tim
###############################################################################
# distribute population in open age group over higher ages.
# The PAS implementation uses stable populations, and it will be added here in the future,
# as well as other optiond. The main missing piece is a good collection of model lifetables.

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
    standard   <- rescale.vector(StPop[StAge >= OAnow], scale = 1)
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
