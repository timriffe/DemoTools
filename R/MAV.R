# Author: Juan Galeano
# Edited 9-Dec-2017 by Tim Riffe
# redundant with utils.R ma(), which takes any n
###############################################################################

#' calculate the moving average (mav) over 3 or 5 years
#' @description  This arithmetic smoothing aims to eliminate irregularities of the population pyramid 
#' by single ages due to errors (attractions or repulsions for certain ages). 

#' @param Value numeric. A vector of demographic counts in single age groups.
#' @param nyears numeric. A single number, (3 or 5), indicating the number of years taken to smooth the population distribution by single ages.
#' @param Age integer lower age bounds
#' @details Single year age groups are assumed.

#' @return a named vector with the smoothed values
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
#' Age  <-c(0:69,'70+')  
#' # final age group assumed open
#' mav(Pop, n = 3,Age=Age)      
mav <- function(Value, nyears, Age){
  # TR: not sure why nyears needs to be hard coded
  stopifnot(nyears %in% c(3, 5))
  result <- as.numeric(stats::filter(Value,rep(1 / nyears, nyears)))
  structure(result, names = Age)
}      


