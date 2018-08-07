# Author: Juan Galeano
# Edited 9-Dec-2017 by Tim Riffe
# redundant with utils.R ma(), which takes any n

# TODO: make Age detected from names if not given. nyear default value.
# add OAG arg (discard value if TRUE). Merge with ma(). This belongs in utils.R
###############################################################################

#' Calculate the moving average (mav) over 3 or 5 years.
#' @description  This arithmetic smoothing technique aims to eliminate irregularities of the population pyramid by averaging values in a moving window of user-defined width. 

#' @param Value numeric. A vector of demographic counts in single age groups.
#' @param n integer. A single number, (often 3 or 5), indicating the number of years taken to smooth the population distribution by single ages.
#' @param Age integer. A vector of ages corresponding to the lower integer bound of the counts. 
#' @param OAG logical. Whether or not the top age group is open. Default \code{TRUE}. 
#' @details Ages may be single or grouped, but all age intervals are assumed equal.

#' @return Vector with the smoothed demographic counts.
#' @export
#' @author Juan Galeano

#' @references 
#' \insertRef{GDA1981IREDA}{DemoTools}

#' @examples 
#' # data from GDA1981IREDA, Table 11, page 361-363: Morocco (Census of 1960)
#' Pop  <-c(303583,390782,523903,458546,517996,400630,485606,325423,471481,189710,
#'          385442,143205,270890,145105,138078,157444,153035,91566,247160,73115,
#'          384222,83551,198555,111347,129851,477510,149272,100814,178465,50684,
#'          577167,51878,97788,55544,58011,393200,85048,51131,80336,31246,
#'          454698,34864,51810,31146,26618,228718,38504,23616,40836,15589,
#'          339158,21349,26997,17590,17513,119763,22704,12336,17635,8485,
#'          323263,9535,13906,9063,8294,90459,9817,6376,8884,3773,160609)
#' Age  <-c(0:70)  
#' # final age group assumed open
#' mav(Pop, n = 3,Age=Age)      
mav <- function(Value, Age, n = 3, OAG = TRUE){
  In <- Value
  if (missing(Age)){
	  Age <- as.integer(names(Value))
  }	
  if (OAG){
	  Value[length(Value)] <- NA
  }
  # TR: not sure why n needs to be hard coded
  Out <- ma(Value, n)
  
  structure(Out, names = Age)
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
#agesmth(Pop,Age,method = "un")
