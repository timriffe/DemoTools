# Author: sean
# TR: redundant with utils.R rescale.vector(), which is depended on elsewhere
###############################################################################
#
#Males      <- c(423, 1654, 1523, 1475, 1482, 1356, 1321, 1298, 1121, 1023, 974, 968, 859, 824, 756, 723, 652, 622)
#Females    <- c(401, 1608, 1498, 1450, 1465, 1345, 1338, 1295, 1156, 1098, 1036, 987, 867, 882, 755, 746, 823, 875)
#DesiredTotalMale       <- 23657
#DesiredTotalFemale     <- 24764

#' Proportionally adjust population or death counts by age to a given total.
#' @description The adjustment factor for each age is the ratio of the
#' desired counts total divided by the initial counts total. This
#' adjustment factor is used to adjust each age group to the new given total.
#' This comes from the PAS spreadsheet called ADJAGE.

#' @param Value   numeric. A vector of population or death counts by age.
#' @param DesiredTotal integer. An integer giving the desired total counts.

#' @details The age group structure of the output is the same as that of the
#' input. This function does not adjust for unknown age groups.

#' @return A vector of the adjusted counts by age.
#' @references
#' \insertRef{PAS}{DemoTools}
#' @export
#'
#' @examples
#' Males      <- c(423, 1654, 1523, 1475, 1482, 1356, 1321, 1298, 1121,
#'                 1023, 974, 968, 859, 824, 756, 723, 652, 622)
#' Females    <- c(401, 1608, 1498, 1450, 1465, 1345, 1338, 1295, 1156,
#'                 1098, 1036, 987, 867, 882, 755, 746, 823, 875)
#' DesiredTotalMale       <- 23657
#' DesiredTotalFemale     <- 24764
#' adjustAge(Males, DesiredTotalMale)
#' adjustAge(Females, DesiredTotalFemale)
#' adjustAge(Males+Females, DesiredTotalMale+DesiredTotalFemale)
#' \dontrun{
#'plot(adjustAge(Males, DesiredTotalMale), t ='l', col='red',
#'     ylim = c(400,2100), ylab = 'The counts', xlab = 'Age groups')
#'lines(Males, t = 'l', col='black')
#'legend(15,2000, legend = c('Original','Adjusted'),
#'       col=c('black','red'),
#'       lty = 1)
#' }

adjustAge <- function(Value, DesiredTotal) {
  adjustValues <- rescale.vector(Value, DesiredTotal)
  
  return(adjustValues)
}
