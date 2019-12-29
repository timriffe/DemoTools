# Author: MJA
###############################################################################


#' Estimate Rogers-Castro migration age schedule

#' @description Given a set of ages and observed age-specific migration rates, estimate the parameters of a Rogers-Castro model migration schedule. 
#' Choose between a 7,9,11 or 13 parameter model. 

#' @param ages numeric. A vector of ages. 
#' @param mx numeric. A vector of observed age-specific migration rates. 
#' @param num_pars integer. Number of parameters to be estimated in the model. 
#' @param ... additional inputs to stan, see ?rstan::stan for details. 
#' @importFrom rstan stan extract
#' @importFrom stats quantile
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom tidyr gather
#' @importFrom rlang sym
#' @export


mig_estimate_rc <- function(ages, 
                            mx, 
                            num_pars,
                            ...){
  
  stopifnot(num_pars %in% c(7, 9, 11, 13))
  
  # data for model input
  y <- mx
  x <- ages
  
  mig_data <- list(
    N = length(x),
    y = y,
    x = x
  )
  
  # stan model
  rc_model <- .return_rc_stan_model(num_pars = num_pars)
  
  # fit the model
  rc_fit <- rstan::stan(
    model_code = rc_model,
    data = mig_data,    # named list of data
    ...
  )
  
  # extract the posterior samples
  list_of_draws <- rstan::extract(rc_fit)
  
  # create a matrix to store fitted values
  y_hat      <- matrix(nrow = length(list_of_draws[[1]]), ncol = length(x))
  these_pars <- list()
  parnames   <- names(list_of_draws)[grep("alpha|a[0-9]|mu[0-9]|lambda|^c$",names(list_of_draws))]
  for(j in 1:length(list_of_draws[[1]])){
    for(i in 1:length(parnames)){
      these_pars[[names(list_of_draws)[i]]] <- list_of_draws[[names(list_of_draws)[i]]][j]
    }
    y_hat[j,] <- mig_calculate_rc(ages = ages, pars = these_pars)
  }
  
  dfit <- tibble::tibble(age = x, 
                 data = y, median = apply(y_hat, 2, median),
                 lower = apply(y_hat, 2, quantile,0.025),
                 upper = apply(y_hat, 2, quantile, 0.975),
                 diff_sq = (!!sym("median") - !!sym("data"))^2)

  pars_df <- rc_fit %>% tidybayes::spread_draws(`a[0-9]`,
                                                `alpha[0-9]`,
                                                `mu[0-9]`,
                                                `lambda[0-9]`,
                                                `^c$`,
                                                regex = TRUE) %>%
    # TR: can we plz switch to pivot_longer and add @imporFrom dplyr pivot_longer ?
    # I didn't switch this cuz I didn't have live data to test it on.
    gather(!!sym("variable"), !!sym("value"), -!!sym(".chain"), -!!sym(".iteration"), -!!sym(".draw")) %>%
    group_by(!!sym("variable")) %>%
    summarise(median = median(!!sym("value")),
              lower = quantile(!!sym("value"), 0.025),
              upper = quantile(!!sym("value"), 0.975))
  
  return(list(pars_df = pars_df, fit_df = dfit))
}
