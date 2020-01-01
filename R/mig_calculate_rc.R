
# Author: MJA
###############################################################################


#' Calculate Rogers-Castro migration age schedule

#' @description Given a set of ages and parameters, calculate the migration age schedule based on the Rogers and Castro formula. 
#' Choose between a 7,9,11 or 13 parameter model. 

#' @param ages numeric. A vector of ages for migration rates to be calculated. 
#' @param pars numeric. A named list of parameters. Must have 7, 9, 11 or 13 values. 
#' @export

#' @details In the full 13 parameter model, the migration rate at age x, \eqn{m(x)} is defined as
#' \deqn{m(x) = a1*exp(-1*alpha1*x) + a2*exp(-1*alpha2*(x - mu2) - exp(-1*lambda2*(x - mu2))) + a3*exp(-1*alpha3*(x - mu3) - exp(-1*lambda3*(x - mu3))) + a4*exp(lambda4*x) + c}
#' 
#' The first, second, third and fourth pieces of the equation represent pre-labour force, working age, retirement and post retirement age patterns, respectively. Models with less parameters gradually remove terms at the older ages. 
#' @references
#' \insertRef{rogers1981model}{DemoTools}
#' @examples 
#' pars <- c(a1= 0.09, alpha1= 0.1, a2= 0.2, 
#' alpha2= 0.1, mu2= 21, lambda2= 0.39, a3= 0.001, 
#' alpha3= 1, mu3= 67, lambda3= 0.6, c= 0.01)
#' ages <- 0:75
#' mx <- mig_calculate_rc(ages = ages, pars = pars)
#' \dontun{
#' plot(ages, mx, type = 'l')
#'}
mig_calculate_rc <- function(ages,
                             pars){
  
  # parameter name groups
  comp1 <- c("a1", "alpha1")
  comp2 <- c("a2", "alpha2", "lambda2", "mu2")
  comp3 <- c("a3", "alpha3", "lambda3", "mu3")
  comp4 <- c("a4", "lambda4")
  
  
  # simple check
  stopifnot(length(pars) %in% c(7, 9, 11, 13))
  # check for specific parameter groups
  if (any(comp1 %in% names(pars))){
    stopifnot(all(comp1 %in% names(pars)))
  }
  if (any(comp2 %in% names(pars))){
    stopifnot(all(comp2 %in% names(pars)))
  }
  if (any(comp3 %in% names(pars))){
    stopifnot(all(comp3 %in% names(pars)))
  }
  if (any(comp4 %in% names(pars))){
    stopifnot(all(comp4 %in% names(pars)))
  }
  
  pars_blank <- c(a1 = 0, alpha1 = 0, 
                  a2 = 0, alpha2 = 0, mu2 = 0, lambda2 = 0, 
                  a3 = 0, alpha3 = 0, mu3 = 0, lambda3 = 0, 
                  a4 = 0, lambda4 = 0, 
                  c = 0)
  
  pars_blank[names(pars)] <- pars
  pars       <- pars_blank

  x  <- ages
  mx <- 
    # pre labor
    pars[["a1"]]*exp(-1 * pars[["alpha1"]]*x) + 
    
    # working
    pars[["a2"]]*exp(-1 * pars[["alpha2"]] * (x - pars[["mu2"]]) - 
                     exp(-1 * pars[["lambda2"]] * (x - pars[["mu2"]]))) + 
    
    # retirement
    pars[["a3"]] * exp(-1 * pars[["alpha3"]] * (x - pars[["mu3"]]) - 
                       exp(-1 * pars[["lambda3"]] * (x - pars[["mu3"]]))) + 
    
    # post-retirement
    pars[["a4"]] * exp(pars[["lambda4"]] *x ) + 
    
    # intensity parameter
    pars[["c"]]
  
  return(mx)
}