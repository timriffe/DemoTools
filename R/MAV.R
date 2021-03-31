# Author: JG
# Edited 9-Dec-2017 by TR
# Edited Aug 2018 by TR

###############################################################################

#' Calculate the moving average (mav) over 3 or 5 years.
#' @description  This arithmetic smoothing technique aims to eliminate irregularities of the population pyramid by averaging values in a moving window of user-defined width.
#' @details
#'
#' The moving window is applied symmetrically. Data endpoints are imputed with \code{NA}s in output: the is nothing under 0 or over the highest closed age group to average with. The open age group is imputed with \code{NA} prior to calculations, since it cannot be averaged into the next lowest group. For example, for \code{n=3}, age 0 will be \code{NA}, as will the open age group and next lowest age. Age intervals are assumed uniform. This function could be used with either single or 5-year age groups.
#'
#'
#' If \code{tails} is set to \code{TRUE}, it returns the same result as
#' \code{mav} except for both tails where a cascading smoothing is applied.

#' @param Value numeric. A vector of demographic counts in single age groups.
#' @param n integer. A single number, (often 3 or 5), indicating the number of years taken to smooth the population distribution by single ages.
#' @param Age integer. A vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}.
#' @param tails logical. If set to TRUE, a smoothing is applied on both tails
#' such that all values are non-NA. If FALSE (default), tails are set to NA
#' due to the lag of moving averages.
#'
#' @details Ages may be single or grouped, but all age intervals are assumed equal.

#' @return Vector with the smoothed demographic counts.
#'
#' @export
#' @author Juan Galeano

#' @references
#' \insertRef{GDA1981IREDA}{DemoTools}

#' @examples
#'Pop  <-c(303583,390782,523903,458546,517996,400630,485606,325423,471481,189710,
#'		385442,143205,270890,145105,138078,157444,153035,91566,247160,73115,
#'		384222,83551,198555,111347,129851,477510,149272,100814,178465,50684,
#'		577167,51878,97788,55544,58011,393200,85048,51131,80336,31246,
#'		454698,34864,51810,31146,26618,228718,38504,23616,40836,15589,
#'		339158,21349,26997,17590,17513,119763,22704,12336,17635,8485,
#'		323263,9535,13906,9063,8294,90459,9817,6376,8884,3773,160609)
#'Age  <- 0:70
#'# final age group assumed open
#' mav(Pop, n = 3, Age = Age)
#'
#'\dontrun{
#'	nwindows <- sapply(seq(3, 11, by = 2),mav, Value = Pop, Age = Age)
#'	cols     <- gray(seq(.8, 0, length = 5))
#'	lwds     <- seq(3, 1, length = 5)
#'	plot(Age,Pop, col = "red", xlab = "Age", ylab = "The counts", pch=16,
#'			main = "Moving average windows and extreme heaping")
#'	matplot(Age,nwindows,type='l',col=cols, lwd=lwds, add=TRUE,lty=1)
#'	legend("topright",
#'			lty=1,
#'			col = cols,
#'			lwd = lwds,
#'			legend = paste0("n=",seq(3,11,by=2)))
#'}
#'
#' # For cascading smoothing on the tails:
#' mav(Pop, Age, tails = TRUE)
#'
mav <- function(Value, Age, n = 3, OAG = TRUE, tails = FALSE) {
  In <- Value
  if (missing(Age)) {
    Age <- as.integer(names(Value))
  }
  if (OAG) {
    OrigOAGpop <- Value[length(Value)]
    Value[length(Value)] <- NA
  }

  # TR: not sure why n needs to be hard coded
  Out <- ma(Value, n)
  Out <- if (tails) mav_tails(Value, Age, Out, n, OAG) else Out

  if (tails && OAG) Out[length(Value)] <- OrigOAGpop

  structure(Out, names = Age)
}

# Not exported since it can be called with tails = FALSE on mav.
mav_tails <- function(Value, Age, MavOut, n = 3, OAG=TRUE) {
  NewMavOut <- MavOut

  #Last should point to last age group to use
  Last <- length(Age)
  if (OAG) {
    Last <- Last - 1
    NewMavOut[Last+1] <- Value[Last+1]
  }

  NewMavOut[1] <-Value[1]
  NewMavOut[Last] <- Value[Last]

  MavLev <- c(1,2,4,6,8)

  if (n >= 2)   {
    for(i in 2:(as.integer(n/2))) {
        NewMavOut[i] <- ma(Value[1:(MavLev[i]+1)], MavLev[i])[i]
        # subscripts right and select just the correct age
        NewMavOut[Last - i + 1] <- ma(Value[(Last - MavLev[i] ):Last], MavLev[i])[ i ]

    }
 }

  NewMavOut
}


#Pop  <-c(303583,390782,523903,458546,517996,400630,485606,325423,471481,189710,
#		385442,143205,270890,145105,138078,157444,153035,91566,247160,73115,
#		384222,83551,198555,111347,129851,477510,149272,100814,178465,50684,
#		577167,51878,97788,55544,58011,393200,85048,51131,80336,31246,
#		454698,34864,51810,31146,26618,228718,38504,23616,40836,15589,
#		339158,21349,26997,17590,17513,119763,22704,12336,17635,8485,
#		323263,9535,13906,9063,8294,90459,9817,6376,8884,3773,160609)
#Age  <-c(0:70)
#mav(Pop,Age,OAG=FALSE)
#groupAges(mav(Pop,Age,OAG=FALSE))-
#smooth_age_5(Pop,Age,method = "un")
