
#' Calculate the moving average (mav) over 3 or 5 years.
#' @description  This arithmetic smoothing technique aims to eliminate irregularities of the population pyramid by averaging values in a moving window of user-defined width.
#' @details
#' The moving window is applied symmetrically. By default (`tails = FALSE`) data endpoints are imputed with `NA`s in output: the is nothing under 0 or over the highest closed age group to average with. The open age group is not used in averaging, and it is returned as-is. Age intervals are assumed uniform. This function could be used with either single or 5-year age groups.
#'
#' If `tails` is set to `TRUE`, then tails have been imputed using moving averages with successively smaller values of `n`, the cascade method.

#' @param Value numeric. A vector of demographic counts in single age groups.
#' @param n integer. A single number, (often 3 or 5), indicating the number of years taken to smooth the population distribution by single ages.
#' @param Age integer. A vector of ages corresponding to the lower integer bound of the counts.
#' @param OAG logical. Whether or not the top age group is open. Default `TRUE`.
#' @param tails logical. If set to `TRUE`, smaller-n moving averages are applied on both tails
#' such that all values are non-NA. If `FALSE` (default), tails are set to NA
#' due to the lag of moving averages.

#' @return Vector with the smoothed demographic counts.
#'
#' @export

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
#'  odds <- seq(3, 11, by = 2)
#'	nwindows <- sapply(odds,
#'	                   mav, 
#'	                   Value = Pop, 
#'	                   Age = Age,
#'	                   OAG = TRUE,
#'	                   tails = FALSE)
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
#'\dontrun{
#'# Compare
#' nwindows_tails <- sapply(odds,
#'                          mav, 
#'                          Value = Pop, 
#'                          Age = Age, 
#'                          OAG = TRUE, 
#'                          tails = TRUE)
#' 
#' colnames(nwindows)        <- odds
#' colnamaes(nwindows_tails) <- odds
#' 
#' # NA triangles are completed with
#' # successively smaller ns.
#' head(nwindows)
#' head(nwindows_tails)
#' 
#' tail(nwindows)
#' tail(nwindows_tails)
#' }

mav <- function(Value, Age, n = 3, OAG = TRUE, tails = FALSE) {
  In <- Value
  if (missing(Age)) {
    Age <- as.integer(names(Value))
  }
  
  # save OAG 
  if (OAG) {
    OrigOAGpop           <- Value[length(Value)]
    Value[length(Value)] <- NA
  }

  Out <- ma(Value, n)
  
  # apply cascading tails if needed
  if (tails){
    Out <- mav_tails(Value = Value, 
                     Age = Age, 
                     MavOut = Out, 
                     n = n, 
                     OAG = OAG) 
  }
  # plug OAG back in
  if (OAG){
    Out[length(Value)] <- OrigOAGpop 
  }

  structure(Out, names = Age)
}

# Not exported since it can be called with tails = FALSE on mav.
mav_tails <- function(Value, Age, MavOut, n = 3, OAG = TRUE) {
  NewMavOut <- MavOut

  #Last should point to last age group to use
  Last <- length(Age)
  if (OAG) {
    Last <- Last - 1
    NewMavOut[Last+1] <- Value[Last+1]
  }

  NewMavOut[1]    <- Value[1]
  NewMavOut[Last] <- Value[Last]

  MavLev <- c(1,2,4,6,8)

  if (n >= 2)   {
    for(i in 2:(as.integer(n/2))) {
      
      # TR: why not just calculate the whole thing and once and pick out
      # the two values as needed?
        NewMavOut[i] <- ma(Value[1:(MavLev[i]+1)], n = MavLev[i])[i]
        # subscripts right and select just the correct age
        NewMavOut[Last - i + 1] <- ma(Value[(Last - MavLev[i] ):Last], 
                                      n = MavLev[i])[ i ]

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
