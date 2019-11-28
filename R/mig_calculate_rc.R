
# Author: MJA
###############################################################################


#' Calculate Rogers-Castro migration age schedule

#' @description Given a set of ages and parameters, calculate the a model migration age schedule based on the Rogers and Catsro formula. 
#' Choose between a 7,9,11 or 13 parameter model (TODO, only does 11 for now). 

#' @param ages numeric. A vector of ages for migration rates to be calculated. 
#' @param pars numeric. A named list of parameters. 
#' @param num_pars integer. The number of parameters in the model (7,9,11, or 13)
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
#' mx <- mig_calculate_rc(ages = ages, pars = pars, num_pars = 11)
#' plot(ages, mx, type = 'l')

mig_calculate_rc <- function(ages,
                             pars,
                             num_pars = 11){
  pars_blank <- c(a1 = 0, alpha1 = 0, a2 = 0, 
               alpha2 = 0, mu2 = 0, lambda2 = 0, a3 = 0, 
               alpha3 = 0, mu3 = 0, lambda3 = 0, c = 0)
  pars_blank[names(pars)] <- pars
  pars       <- pars_blank
  # TR: We can now do the whole eq even if all pars aren't given.
  # so maybe no need for this check? Is there some other test of 
  # validity that makes sense? ranges?
  if(length(pars) != num_pars){
    stop("Incorrect number of parameters specified.")
  }
  # TR: this won't work in a package framework, won't pass tests.
  # Will extract pars the old fashioned way?
  # for (i in 1:length(parameters)){
  #   assign(names(parameters)[i], parameters[[i]])
  # }
  
  x  <- ages
  mx <- 
    # pre labor
    pars["a1"]*exp(-1 * pars["alpha1"]*x) + 
    
    # working
    pars["a2"]*exp(-1 * pars["alpha2"] * (x - pars["mu2"]) - 
                     exp(-1 * pars["lambda2"] * (x - pars["mu2"]))) + 
    
    # retirement
    pars["a3"] * exp(-1 * pars["alpha3"] * (x - pars["mu3"]) - 
                       exp(-1 * pars["lambda3"] * (x - pars["mu3"]))) + 
    
    # post retirement?
    pars["c"]
  
  return(mx)
}