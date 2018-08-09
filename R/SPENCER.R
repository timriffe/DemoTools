# Author: Juan Galeano
###############################################################################
# TODO add ref to bibtex
#' Smoothing of an age structure by single years using Spencer's formula.

#' @description  Spencer's method is used to adjust distributions by year of age based on cumulative series. 
#' 
#' @param Value numeric. A vector of demographic counts in single age groups.
#' @param Age numeric or character. A vector with ages in single years.
#' @details Single year age groups are assumed.

#' @return A named vector with the smoothed age structure.
#' @importFrom stats filter
#' @export
#' @references 
#' \insertRef{GDA1981IREDA}{DemoTools}
#' Roger et al (1981, pp. 373-377).
#' @examples 
#' # data from GDA1981IREDA, Table 15 and 16, page 375-376: Ghana (1970 Males + Females), 
#' # population between 10 and 40 years by single age.
#' Pop  <-c(286487,282182,319718,345610,329133,324683,330555,269536,271593,253798,
#'          242297,165755,233159,185812,175647,187138,160864,128099,178568,123386,
#'          194999,113372,139088,106190,127482,178375,123286,95637,151514,82614,
#'          239243,74788,112354,61746,72366,144828,92563,53924,94635,52351,178832)
#' Age  <-c(0:40) 
#' sp_smoothed <- spencer(Value=Pop,Age=Age)
#'  \dontrun{
#' plot(Age,Pop, col = "red", xlab = "Age", ylab = "Population",type ='l')
#' lines(Age,spencer(Pop,Age), col = "blue")
#' legend("topright", 
#' lty=1,
#' col = c("red","blue"), 
#' legend =  c("Original","Spencer"))
#' }

#' # on the fly unit test:
#' Tab16answer <- c(rep(NA,10),233847,216876,201938,189873,179734,171529,164280,158214,
#'  152477,147735,143099,139301,135930,133889,132458,132495,
#'  132660,132894,131463,128560,123150,116725,109180,102386,
#'  96394,92586,89962,88681,86812,84209,79819)
#' stopifnot(all(na.omit(abs(sp_smoothed - Tab16answer) < 0.5)))

spencer <- function(Value, Age){
  
  # ni3= The sum of the ages of three consecutive ages, centered on age i.
  ni3  <- c(NA,stats::filter(Value, rep(1, 3L))[-c(1, length(Value))], NA)
  
  # ni5= The sum of the ages of five consecutive ages, centered on age i.
  ni5  <- c(NA,stats::filter(Value, rep(1, 5L))[-c(1, length(Value))], NA)
  
  # ni7= The sum of the ages of seven consecutive ages, centered on age i.
  ni7  <- c(NA,stats::filter(Value, rep(1, 7L))[-c(1, length(Value))], NA)
  
  # Ni= 
  Ni   <- Value + ni3 + ni5 - ni7
  
  # n5i= The sum of ni population for five consecutive ages, centered on age i.
  n5i  <- c(NA,stats::filter(Ni, rep(1, 5L))[-c(1, length(Ni))], NA)
  
  # n55i= The sum of the n51 population for five consecutive ages, centered on age i
  n55i <- c(NA,stats::filter(n5i, rep(1, 5L))[-c(1, length(n5i))],NA)
  
  # n75i= The sum of the n55i population for seven consecutive ages, centered on age i
  n75i <- c(NA,stats::filter(n55i, rep(1, 7L))[-c(1, length(n55i))],NA)
  
  #Ani=
  Ani <- 1 / 350 * n75i
  structure(Ani, names = Age)
}
