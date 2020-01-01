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
#' @import Rcpp
#' @importFrom stats quantile
#' @importFrom dplyr group_by summarise rename
#' @importFrom tidyr gather
#' @importFrom rlang sym
#' @export
#' @examples 
#' # define ages and migration rates
#' ages <- 0:75
#' mig_rate <- c(0.1014,0.0984,0.0839,0.0759,0.0679,0.0616,
#' 0.055,0.0518,0.0438,0.0399,0.0363,0.0342,0.0307,0.0289,
#' 0.0267,0.0237,0.0283,0.0294,0.0392,0.0547,0.068,0.0933,
#' 0.1187,0.1282,0.1303,0.1293,0.1247,0.1163,0.1059,0.0976,
#' 0.09,0.0817,0.0749,0.0683,0.0626,0.058,0.0536,0.0493,
#' 0.0459,0.0418,0.0393,0.0358,0.0343,0.0322,0.0298,0.0281,
#' 0.0266,0.0248,0.0235,0.0222,0.0208,0.0191,0.0178,0.0171,
#' 0.0157,0.0149,0.0139,0.0132,0.012,0.0112,0.0116,0.0106,
#' 0.0102,0.0109,0.0107,0.0143,0.0135,0.0134,0.0116,0.0099,
#' 0.0093,0.0083,0.0078,0.0067,0.0069,0.0054)
#' # fit the model
#' \dontrun{
#' res <- mig_estimate_rc(ages, mig_rate, num_pars = 7)
#' # plot the results and data
#' plot(ages, mig_rate, ylab = "migration rate", xlab = "age")
#' lines(ages, res[["fit_df"]]$median, col = "red")
#' legend("topright", legend=c("data", "fit"), col=c("black", "red"), lty=1, pch = 1)
#' }
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
  pars_df <- rc_fit %>% tidybayes::gather_draws(!!sym("a[0-9]"),
                                                !!sym("alpha[0-9]"),
                                                !!sym("mu[0-9]"),
                                                !!sym("lambda[0-9]"),
                                                !!sym("^c$"),
                                                regex = TRUE) %>%
    group_by(!!sym(".variable")) %>%
    summarise(median = median(!!sym(".value")),
              lower = quantile(!!sym(".value"), 0.025),
              upper = quantile(!!sym(".value"), 0.975)) %>% 
    dplyr::rename("variable" = !!sym(".variable"))
  
  return(list(pars_df = pars_df, fit_df = dfit))
}
