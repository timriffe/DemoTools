
# Author: MJA
###############################################################################


#' Calculate Rogers-Castro migration age schedule

#' @description Given a set of ages and parameters, calculate the a model migration age schedule based on the Rogers and Catsro formula. 
#' Choose between a 7,9,11 or 13 parameter model (TODO, only does 11 for now). 

#' @param ages numeric. A vector of ages for migration rates to be calculated. 
#' @param parameters numeric. A named list of parameters. 
#' @param num_parameters integer. The number of parameters in the model (7,9,11, or 13)
#' @export

#' @details In the full 13 parameter model, the migration rate at age x, \eqn{m(x)} is defined as
#' \deqn{m(x) = a1*exp(-1*\alpha1*x) + a2*exp(-1*\alpha2*(x - mu2) - exp(-1*\lambda2*(x - mu2))) + a3*exp(-1*\alpha3*(x - mu3) - exp(-1*\lambda3*(x - mu3))) + a4*exp(\lambda4*x) + c}
#' The first, second, third and fourth pieces of the equation represent pre-labour force, working age, retirement and post retirement age patterns, respectively. Models with less parameters gradually remove terms at the older ages. 
#' @references
#' \insertRef{rogers1981model}{DemoTools}
#' @examples 
#' parameters <- list(a1= 0.09, alpha1= 0.1, a2= 0.2, alpha2= 0.1, mu2= 21, lambda2= 0.39, a3= 0.001, alpha3= 1, mu3= 67, lambda3= 0.6, c= 0.01)
#' ages <- 0:75
#' mx <- mig_calculate_rc(ages = ages, parameters = parameters, num_parameters = 11)
#' plot(ages, mx, type = 'l')

mig_calculate_rc <- function(ages,
                             parameters,
                             num_parameters = 11){
  
  if(length(parameters)!=num_parameters){
    stop("Incorrect number of parameters specified.")
  }
  for (i in 1:length(parameters)){
    assign(names(parameters)[i], parameters[[i]])
  }
  
  x <- ages
  mx <- a1*exp(-1*alpha1*x) + a2*exp(-1*alpha2*(x - mu2) - exp(-1*lambda2*(x - mu2))) + a3*exp(-1*alpha3*(x - mu3) - exp(-1*lambda3*(x - mu3))) + c
  return(mx)
}