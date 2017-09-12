# Author: Juan Galeano
###############################################################################

#' smoothing of an age structure by single years using the formula of Spencer

#' @description  Spencer's method is used to adjust distributions by year of age based on cumulative series.This comes from 
#' pp. 373-377 in GDA_1981_Structures-par-sexe-et-age-en-Afrique_[IREDA]

#' @param Value numeric. A vector of demographic counts in single age groups.
#' @param Age numeric or character. A vector with ages (in single years).
#' @details Single year age groups are assumed.

#' @return a named vector with the smoothed values
#' @export

#' @references 
#' \insertRef{GDA1981IREDA}{DemoTools}

#' @examples 
#' # data from GDA1981IREDA, Table 15 and 16, page 375-376: Ghana (1970 Males + Females), 
#' # population between 10 and 40 years by single age.
#' Pop  <-c(286487,282182,319718,345610,329133,324683,330555,269536,271593,253798,
#'          242297,165755,233159,185812,175647,187138,160864,128099,178568,123386,
#'          194999,113372,139088,106190,127482,178375,123286,95637,151514,82614,
#'          239243,74788,112354,61746,72366,144828,92563,53924,94635,52351,178832)
#' Age  <-c(0:40) 
#' spencer(Pop,Age=Age)
spencer <- function(Value,Age){
  
  # ni3= The sum of the ages of three consecutive ages, centered on age i.
  ni3<-c(NA,stats::filter(Value, rep(1, 3L))[-c(1, length(Value))],NA)
  
  # ni5= The sum of the ages of five consecutive ages, centered on age i.
  ni5<-c(NA,stats::filter(Value, rep(1, 5L))[-c(1, length(Value))],NA)
  
  # ni7= The sum of the ages of seven consecutive ages, centered on age i.
  ni7<-c(NA,stats::filter(Value, rep(1, 7L))[-c(1, length(Value))],NA)
  
  # Ni= 
  Ni<- Value+ni3+ni5-ni7
  
  # n5i= The sum of ni population for five consecutive ages, centered on age i.
  n5i<-c(NA,stats::filter(Ni, rep(1, 5L))[-c(1, length(Ni))],NA)
  
  # n55i= The sum of the n51 population for five consecutive ages, centered on age i
  n55i<-c(NA,stats::filter(n5i, rep(1, 5L))[-c(1, length(n5i))],NA)
  
  # n75i= The sum of the n55i population for seven consecutive ages, centered on age i
  n75i<-c(NA,stats::filter(n55i, rep(1, 7L))[-c(1, length(n55i))],NA)
  
  #Ani=
  Ani=1/350*n75i
  structure(Ani, names=Age)
}
