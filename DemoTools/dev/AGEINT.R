# Author: sean
###############################################################################
#
#EarlyPop      <- c(100958, 466275, 624134, 559559, 446736, 370653, 301862, 249409, 247473, 223014, 172260, 149338, 127242, 105715, 79614, 53660, 31021, 34596)
#LaterPop    <- c(201916, 932550, 1248268, 1119118, 893472, 741306, 603724, 498818, 494946, 446028, 344520, 298676, 254484, 211430, 159228, 107320, 62042, 69192)
#EarlyDate       <- "1980-04-07"
#LaterDate     <- "1990-10-06"
#DesiredDate    <- "1990-07-01"

#' Interpolate between two population age distributions.
#' @description The interpolation is done by age (not cohort) using a linear or exponential function. This comes from the PAS spreadsheet called ADJINT.

#' @param EarlyValue   numeric. A vector of demographic population counts for an earlier date.
#' @param LaterValue   numeric. A vector of demographic population counts for a later date.
#' @param EarlyDate date (See details for ways to express it). The date corresponding to the earlier population age distribution.
#' @param LaterDate date (See details for ways to express it). The date corresponding to the later population age distribution.
#' @param DesiredDate date (See details for ways to express it). The desired date of the output population age distribution.
#' @param method string. The method to use for the interpolation, either "linear" or "exponential". Default "Linear".

#' @details The age group structure of the output is the same as that of the input. Ideally, the DesiredDate should be between the EarlyDate and LaterDate. Dates can be given in three ways 1) a \code{Date} class object, 2) an unambiguous character string in the format \code{"YYYY-MM-DD"}, or 3) as a decimal date consisting in the year plus the fraction of the year passed as of the given date.

#' @return a vector of the interpolated population for the requested date.
#' 
#' @export
#' 
#' @examples 
#' EarlyPop      <- c(100958, 466275, 624134, 559559, 446736, 370653, 301862, 249409, 247473, 223014, 172260, 149338, 127242, 105715, 79614, 53660, 31021, 34596)
#' LaterPop    <- c(201916, 932550, 1248268, 1119118, 893472, 741306, 603724, 498818, 494946, 446028, 344520, 298676, 254484, 211430, 159228, 107320, 62042, 69192)
#' EarlyDate       <- "1980-04-07"
#' LaterDate     <- "1990-10-06"
#' DesiredDate    <- "1990-07-01"
#' 
#' interpolatePop(EarlyPop, LaterPop, EarlyDate, LaterDate, DesiredDate)
#' interpolatePop(EarlyPop, LaterPop, EarlyDate, LaterDate, DesiredDate, method = "exponential")

interpolatePop <- function(EarlyValue, LaterValue, EarlyDate, LaterDate, DesiredDate, method = "linear"){
  
  earlyDateDec <- dec.date(EarlyDate)
  laterDateDec <- dec.date(LaterDate)
  desireDateDec <- dec.date(DesiredDate)
  
  interpolateFactor <- (desireDateDec - earlyDateDec)/(laterDateDec - earlyDateDec)
  
  if (method == "exponential"){
    adjustedPop <- round(EarlyValue*exp(interpolateFactor*log(LaterValue/EarlyValue)))
  }
  else {
    adjustedPop <- round(EarlyValue + (LaterValue-EarlyValue)*interpolateFactor)
  }
  
  return(adjustedPop)
}