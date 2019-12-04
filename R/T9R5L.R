# Author: Juan Galeano
# handle OAG with care.
###############################################################################

#' Feeney'S formula on 9 years to correct for heaping on multiples of 5.

#' @description  Fenney's technique for correcting age distributions for heaping on multiples of five.

#' @param Value numeric. A vector of demographic counts in single age groups.
#' @param Age numeric or character. A vector with ages in single years.
#' @param maxit integer. Maximum number of iterations.
#' @param OAG logical. Is the final age group open? Default \code{FALSE}.

#' @details \code{Value} can be given in single or 5-year age groups.

#' @return a vector of adjusted counts in 5-year age groups
#'
#' @export
#' @references
#' \insertRef{feeney1979}{DemoTools}

#' @examples
#' # data from feeney1979, Table 1, page 12: Population of Indonesia, 22 provinces,
#' # by single year of age: Census of 24 September 1971.
#'  Pop        <- c(2337,3873,3882,3952,4056,3685,3687,3683,3611,3175,
#'          3457,2379,3023,2375,2316,2586,2014,2123,2584,1475,
#'          3006,1299,1236,1052,992,3550,1334,1314,1337,942,
#'          3951,1128,1108,727,610,3919,1221,868,979,637,
#'          3409,887,687,533,313,2488,677,426,524,333,
#'          2259,551,363,290,226,1153,379,217,223,152,
#'          1500,319,175,143,89,670,149,96,97,69,
#'          696,170,60,38,23,745)
#'  Ages       <- c(0:75)
#'  result     <- T9R5L(Pop, Ages, OAG = TRUE)
#'  A5         <- names2age(result)
#'  V5         <- groupAges(Pop,Ages)
#'  \dontrun{
#'  plot(Ages, Pop, type= 'l')
#'  segments(A5,
#'		  result/5,
#'		  A5+5,
#'		 result/5,
#' 		 col = "red")
#' segments(A5,
#'		 V5/5,
#'		 A5+5,
#'		 V5/5,
#'		 col = "blue")
#'  legend("topright",col=c("black","blue","red"),
#'    lty=c(1,1,1),
#'    legend=c("recorded 1","recorded 5","corrected 5"))
#' }

T9R5L          <- function(Value,
                  Age,
                  maxit = 200,
                  OAG = FALSE) {
  # ages need to be single to use this method.
  stopifnot(is_single(Age))
  TOT          <- sum(Value, na.rm = TRUE)
  
  # handle OAG with care
  if (OAG) {
    NN         <- length(Value)
    OAvalue    <- Value[NN]
    OA         <- Age[NN]
    Value      <- Value[-NN]
    Age        <- Age[-NN]
    
  }
  
  V5           <- groupAges(Value, Age, N = 5)
  A5           <- names2age(V5)
  i5           <- Age %% 5 == 0
  
  # TR: is this length-vulnerable?
  V50          <- Value[i5]
  V54          <- V5 - V50
  
  
  # need N anyway
  N            <- length(V5)
  # internal function used iteratively
  f_adjust     <- function(v5, v4) {
    N          <- length(v4)
    Bup        <- c(1, shift.vector(v4,-1, fill = 0) + v4)[1:N]
    FC         <- (8 / 9) * ((v5 + Bup) / Bup)
    FC[1]      <- 1
    POB1       <- c(v5 - (FC - 1) * Bup)
    POB2       <- (c(c(shift.vector(FC,-1, fill = 0) + FC)[-N],
                  (FC[N] + 1)) - 1) * v4
    # return a list of 3 components
    aj         <- list(FC, POB1, POB2)
    aj
  }
  
  
  # adjustment loop
  for (i in 1:maxit) {
    adjust     <- f_adjust(v5 = V50, v4 = V54)
    V50        <- adjust[[2]]
    V54        <- adjust[[3]]
    if (all(abs(adjust[[1]]) < 1e-8)) {
      break
    }
  }
  
  G            <- (V50 * .6)[c(2:(N - 1))]
  H            <- (V50 * .4)[c(3:N)]
  I            <- G + H #f(x+2.5)
  
  # corrected, but unknowns still need to be redistributed
  #I5          <- rescale.vector(I * 5,sum(V5[2:(N-1)]))
  I5           <- I * 5
  out          <- c(V5[1], I5, V5[N])
  
  if (OAG) {
    A5         <- c(A5, OA)
    out        <- c(out, OAvalue)
  }
  names(out)   <- A5
  # rescale to sum, inclusing open age group and boudned tails
  out          <- rescale.vector(out, TOT)
  out
}
