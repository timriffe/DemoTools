# Author: MJA
###############################################################################


#' Estimate Rogers-Castro migration age schedule

#' @description Given a set of ages and observed age-specific migration rates, estimate the parameters of a Rogers-Castro model migration schedule. 
#' Choose between a 7,9,11 or 13 parameter model. 

#' @param ages numeric. A vector of ages. 
#' @param mx numeric. A vector of observed age-specific migration rates. 
#' @export

mig_estimate_rc <- function(ages, 
                            mx, 
                            num_pars = 11,
                            chains = 4,
                            warmup= 1000,
                            iter = 3000,
                            cores = 2,
                            adapt_delta = 0.8, 
                            max_treedepth = 10){
  
  # data for model input
  y <- mx
  x <- ages
  
  mig_data <- list(
    N = length(x),
    y = y,
    x = x
  )
  
  # stan model
  rc_model <- "
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
real<lower=0> alpha1;
real<lower=0> alpha2;
real<lower=0> alpha3;
real<lower=0, upper=1> a1;
real<lower=0, upper=1> a2;
real<lower=0, upper=1> a3;
real<lower=0, upper=1> c;
real<lower=0> mu2;
real<lower=mu2+1, upper=max(x)> mu3;
real<lower=0> lambda2;
real<lower=0> lambda3;
real<lower=0> sigma;
}

transformed parameters {
vector[N] mu_rc;

mu_rc = a1*exp(-alpha1*x) + a2*exp(-alpha2*(x - mu2) - exp(-lambda2*(x - mu2))) + a3*exp(-alpha3*(x - mu3) - exp(-lambda3*(x - mu3))) + c;
}

model {
// likelihood
y ~ normal(mu_rc, sigma);

//priors
alpha1 ~ normal(0,1);
alpha2 ~ normal(0,1);
alpha3 ~ normal(0,1);
a1 ~ normal(0,1);
a2 ~ normal(0,1);
a3 ~ normal(0,1);
c ~ normal(0,1);
mu2 ~ normal(25,1);
mu3 ~ normal(65,1);
lambda2 ~ normal(0,1);
lambda3 ~ normal(0,1);
sigma ~ normal(0,1);
}
"
  
# fit the model
  rc_fit <- rstan::stan(
    model_code = rc_model,
    data = mig_data,    # named list of data
    chains = chains,             # number of Markov chains
    warmup = warmup,          # number of warmup iterations per chain
    iter = iter,            # total number of iterations per chain
    cores = cores,              # number of cores 
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    refresh = 2000          # show progress every 'refresh' iteration
  )
  
  # extract the posterior samples
  list_of_draws <- rstan::extract(rc_fit)
  
  # create a matrix to store fitted values
  y_hat <- matrix(nrow = length(list_of_draws[[1]]), ncol = length(x))
  for(j in 1:length(list_of_draws[[1]])){
      these_pars <- list(a1 = list_of_draws[["a1"]][j],
                         alpha1 = list_of_draws[["alpha1"]][j],
                         a2 = list_of_draws[["a2"]][j],
                         alpha2 = list_of_draws[["alpha2"]][j],
                         a3 = list_of_draws[["a3"]][j],
                         alpha3 = list_of_draws[["alpha3"]][j],
                         mu2 = list_of_draws[["mu2"]][j],
                         mu3 = list_of_draws[["mu3"]][j],
                         lambda2 = list_of_draws[["lambda2"]][j],
                         lambda3 = list_of_draws[["lambda3"]][j],
                         c = list_of_draws[["c"]][j]
      )
      y_hat[j,] <- mig_calculate_rc(ages = ages, pars = these_pars)
  }
  
  
  dfit <- tibble::tibble(age = x, 
                 data = y, median = apply(y_hat, 2, median),
                 lower = apply(y_hat, 2, quantile,0.025),
                 upper = apply(y_hat, 2, quantile, 0.975),
                 diff_sq = (median - data)^2)
  
  
  
  pars_df <- rc_fit %>% tidybayes::spread_draws(a1, a2, a3,
                                             alpha1, alpha2, alpha3,
                                             mu2, mu3, lambda2, lambda3,
                                             c, sigma) %>%
    gather(variable, value, -.chain, -.iteration, -.draw) %>%
    group_by(variable) %>%
    summarise(median = median(value),
              lower = quantile(value, 0.025),
              upper = quantile(value, 0.975))
  
  return(list(pars_df = pars_df, fit_df = dfit))
}