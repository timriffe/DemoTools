# Author: sean
###############################################################################
#
#Males      <- c(423, 1654, 1523, 1475, 1482, 1356, 1321, 1298, 1121, 1023, 974, 968, 859, 824, 756, 723, 652, 622)
#Females    <- c(401, 1608, 1498, 1450, 1465, 1345, 1338, 1295, 1156, 1098, 1036, 987, 867, 882, 755, 746, 823, 875)
#DesiredTotalMale       <- 23657
#DesiredTotalFemale     <- 24764

#' proportionally adjust a population age distribution to a new total.
#' @description The adjustment factor for each age is the ratio of the desired population total divided by the initial population total. This adjustment factor is used to adjust each age group to the new total.
#' This comes from the PAS spreadsheet called ADJAGE.

#' @param Value   numeric. A vector of demographic population counts.
#' @param DesiredTotal integer. An integer giving the desired total population.

#' @details The age group structure of the output is the same as that of the input. This function does not adjust for unknown age groups.

#' @return a vector of the adjusted population.
#' 
#' @export
#' 
#' @examples 
#' Males      <- c(423, 1654, 1523, 1475, 1482, 1356, 1321, 1298, 1121, 1023, 974, 968, 859, 824, 756, 723, 652, 622)
#' Females    <- c(401, 1608, 1498, 1450, 1465, 1345, 1338, 1295, 1156, 1098, 1036, 987, 867, 882, 755, 746, 823, 875)
#' DesiredTotalMale       <- 23657
#' DesiredTotalFemale     <- 24764
#' adjustAge(Males, DesiredTotalMale)
#' adjustAge(Females, DesiredTotalFemale)
#' adjustAge(Males+Females, DesiredTotalMale+DesiredTotalFemale)

adjustAge <- function(Value, DesiredTotal){
  
  valueSum <- sum(Value)
  
  adjustFactor <- DesiredTotal/valueSum
  
  adjustValues <- round(Value*adjustFactor)
  
  return(adjustValues)
}