# Author: Juan Galeano
###############################################################################

#' Zelnik 11-term moving average. Adjusting for digit preference

#' @description  Zelnik method is used to adjust distributions by year of age for digit preference.This comes from 
#' Gray_1987_The Missing Ages: Adjusting for Digit Preference

#' @param Value numeric. A vector of demographic counts in single age groups.
#' @param q to be defined
#' @param Age integer lower bound of age classes

#' @details Single year age groups are assumed.

#' @return a named vector with the adjusted values
#' @export

#' @references 
#' \insertRef{gray1987missingages}{DemoTools}

#' @examples 
#' # data from gray1987missingages, Table 2, page 21: Aplication of Q1 and 
#' # Q2 linear operators, Bangladesh Census, 1 March 1974-Males.
#' Pop<-c(941307,1041335,1237034,1411359,1383853,1541942,1321576,1285877,1563448,886705,
#'       1623998,562924,1485173,543216,771219,903496,686431,370007,942999,250820,
#'       1023667,200131,688640,222011,281738,1239965,288363,263326,483143,78635,
#'       1349886,68438,415127,101596,100758,1392434,178633,126351,286520,50836,
#'       1331036,48995,251153,58393,54995,1033812,68792,72766,175943,28254,
#'       1038747,32894,136179,37667,38230,596049,52602,36493,74106,16759,
#'       790643,20596,70109,18044,19891,357491,15253,17489,31057,8481,
#'       429816,7951,35583,8612,6589,454645)
#' Age  <-c(0:75) 
#' zelnik(Pop,1,Age)
zelnik <- function(Value,q,Age){
  
  # Table 11, pp.20 Linear operators for use in deriving unbiased estimates of single-age distributions
  q1 <-c(-0.0025,-0.01,-0.02,-0.03,-0.04,0.05,0.14,0.13,0.12,0.11,
         0.105,0.11,0.12,0.13,0.14,0.05,-0.04,-0.03,-0.02,-0.01,-0.0025)
  
  q2<-c(-0.00025,-0.0015,-0.0045,-0.0095,-0.0165,-0.018,-0.0065,0.0105,0.0255,0.0385,
        0.05025,0.063,0.079,0.099,0.123,0.136,0.123,0.099,0.079,0.063,0.05025,0.0385,
        0.0255,0.0105,-0.0065,-0.018,-0.0165,-0.0095,-0.0045,-0.0015,-0.00025)
  
  # Multiply the vector of single ages by the linear operators Q1 or Q2
  multi<-Value %*% t(if(q==1){q1}else 
  {if(q==2) {q2}})
  
  # Get the diagonal for each age and split the vector
  diagonal <- row(multi) - col(multi)
  dflist<-split(multi, diagonal)
  
  # Get the unbiased estimates of single-age distributions in a named vector
  Q<-if(q==1) {as.vector(sapply(dflist, sum)[22:length(Age)-1])}else 
    {if(q==2) {as.vector(sapply(dflist, sum)[32:length(Age)-1])}}
  
  structure(Q, names=if(q==1){(min(Age)+10):(max(Age)-11)}else 
                    {if(q==2){(min(Age)+15):(max(Age)-16)}})
}
