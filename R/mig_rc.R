 # Functions to calculate and estimate Rogers-Castro migration age schedules

 # Author: MJA
 ###############################################################################
 
 
 #' Calculate Rogers-Castro migration age schedule
 
 #' @description Given a set of ages and parameters, calculate the migration age schedule based on the Rogers and Castro formula. 
 #' Choose between a 7,9,11 or 13 parameter model. 
 
 #' @param ages numeric. A vector of ages for migration rates to be calculated. 
 #' @param pars numeric. A named list of parameters. See below for details.
 #' @export
 
 #' @details In the full 13 parameter model, the migration rate at age x, \eqn{m(x)} is defined as
 #' \deqn{m(x) = a1*exp(-1*alpha1*x) + a2*exp(-1*alpha2*(x - mu2) - exp(-1*lambda2*(x - mu2))) + a3*exp(-1*alpha3*(x - 3) - exp(-1*lambda3*(x - mu3))) + a4*exp(lambda4*x) + c}
 #' 
 #' The first, second, third and fourth pieces of the equation represent pre-working age, working age, retirement and post-retirement age patterns, respectively. 
 #' Models with less parameters gradually remove terms at the older ages. Parameters in each family are:
 #' \itemize{
 #' \item pre-working age: {a1, alpha1}
 #' \item working age: {a2, alpha2, mu2, lambda2}
 #' \item retirement: {a3, alpha3, mu3, lambda3}
 #' \item post retirement: {a4, lambda4}
 #' }
 #' For a specific family to be included, values for all parameters in that family must be specified. 
 #' 
 #' @references
 #' \insertRef{rogers1981model}{DemoTools}
 #' @examples 
 #' pars <- c(a1= 0.09, alpha1= 0.1, a2= 0.2, 
 #' alpha2= 0.1, mu2= 21, lambda2= 0.39, a3= 0.001, 
 #' alpha3= 1, mu3= 67, lambda3= 0.6, c= 0.01)
 #' ages <- 0:75
 #' mx <- mig_calculate_rc(ages = ages, pars = pars)
 #' \dontrun{
 #' plot(ages, mx, type = 'l')
 #'}
 mig_calculate_rc <- function(ages,
                              pars){
   
   # parameter name groups
   comp1 <- c("a1", "alpha1")
   comp2 <- c("a2", "alpha2", "lambda2", "mu2")
   comp3 <- c("a3", "alpha3", "lambda3", "mu3")
   comp4 <- c("a4", "lambda4")
   
   
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
     # pre working age
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
 

# Author: MJA
###############################################################################

#' Estimate Rogers-Castro migration age schedule
 
#' @description Given a set of ages and observed age-specific migration rates, estimate the parameters of a Roger-Castro model migration schedule. 
#' Choose between a 7,9,11 or 13 parameter model. 

#' @param ages numeric. A vector of ages. 
#' @param mx numeric. A vector of observed age-specific migration rates. 
#' @param pre_working_age logical (TRUE/FALSE). Whether or not to include pre working age component. 
#' @param working_age logical (TRUE/FALSE). Whether or not to include working age component. 
#' @param retirement logical (TRUE/FALSE). Whether or not to include retirement age component. 
#' @param post_retirement logical (TRUE/FALSE). Whether or not to include post retirement age component. 
#' @param ... additional inputs to stan, see ?rstan::stan for details. 
#' @importFrom rstan stan extract
#' @import Rcpp
#' @importFrom stats quantile
#' @importFrom dplyr group_by summarise rename mutate 
#' @importFrom rlang sym
#' @importFrom tibble tibble
#' @importFrom tibble as.tibble
#' @importFrom tidybayes gather_draws
#' @importFrom rstan extract
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
#' 
#' res <- mig_estimate_rc(ages, mig_rate, 
#' pre_working_age = TRUE, 
#' working_age = TRUE, 
#' retirement = FALSE, 
#' post_retirement = FALSE)
#' \dontrun{
#' # plot the results and data
#' plot(ages, mig_rate, ylab = "migration rate", xlab = "age")
#' lines(ages, res[["fit_df"]]$median, col = "red")
#' legend("topright", legend=c("data", "fit"), col=c("black", "red"), lty=1, pch = 1)
#' }
mig_estimate_rc <- function(ages, 
                           mx, 
                           pre_working_age,
                           working_age,
                           retirement,
                           post_retirement,
                           ...){
  
 stopifnot(any(pre_working_age, working_age, retirement, post_retirement))
 
 # data for model input
 y <- mx
 x <- ages
 
 mig_data <- list(
   N = length(x),
   y = y,
   x = x, 
   pre_working_age = as.numeric(pre_working_age),
   working_age = as.numeric(working_age),
   retirement = as.numeric(retirement),
   post_retirement = as.numeric(post_retirement)
 )
 
 # model
 
 rc_flexible <- 'data {
     int<lower=0,upper=1> pre_working_age; // 0 = no, 1 = yes
 int<lower=0,upper=1> working_age; // 0 = no, 1 = yes
 int<lower=0,upper=1> retirement; // 0 = no, 1 = yes
 int<lower=0,upper=1> post_retirement; // 0 = no, 1 = yes
 int<lower=0> N;
 vector[N] x;
 vector[N] y;
}
parameters {
real<lower=0> alpha1[1*pre_working_age];
real<lower=0> alpha2[1*working_age];
real<lower=0> alpha3[1*retirement];
real<lower=0, upper=1> a1[1*pre_working_age];
real<lower=0, upper=1> a2[1*working_age];
real<lower=0, upper=1> a3[1*retirement];
real<lower=0, upper=1> a4[1*post_retirement];
real<lower=0> mu2[1*working_age];
real<lower=0, upper=max(x)> mu3[1*retirement];
real<lower=0> lambda2[1*working_age];
real<lower=0> lambda3[1*retirement];
real<upper=0.05> lambda4[1*post_retirement];
real<lower=0, upper=1> c;
real<lower=0> sigma;
}
 transformed parameters {
 vector[N] mu_rc;
 vector[N] mu_rc_1;
 vector[N] mu_rc_2;
 vector[N] mu_rc_3;
 vector[N] mu_rc_4;
 vector[N] zero;
 
 for(i in 1:N){
 zero[i] = 0;
 }
 
 
 mu_rc_1 = pre_working_age==1?a1[1]*exp(-alpha1[1]*x):zero;
 mu_rc_2 = working_age==1?a2[1]*exp(-alpha2[1]*(x - mu2[1]) - exp(-lambda2[1]*(x - mu2[1]))):zero;
 mu_rc_3 = retirement==1?a3[1]*exp(-alpha3[1]*(x - mu3[1]) - exp(-lambda3[1]*(x - mu3[1]))):zero;
 mu_rc_4 = post_retirement==1?a4[1]*exp(lambda4[1]*(x)):zero;
 mu_rc = mu_rc_1 + mu_rc_2 + mu_rc_3 + mu_rc_4 + c;
 }
 model {
 // likelihood
 y ~ normal(mu_rc, sigma);
 
 //priors
 
 if(pre_working_age==1){
 alpha1 ~ normal(0,1);
 a1 ~ normal(0,0.1);
 }
 if(working_age==1){
 alpha2 ~ normal(0,1);
 a2 ~ normal(0,0.1);
 mu2 ~ normal(25,1);
 lambda2 ~ normal(0,1);
 }
 if(retirement==1){
 alpha3 ~ normal(0,1);
 a3 ~ normal(0,0.1);
 mu3 ~ normal(65,1);
 lambda3 ~ normal(0,1);
 }
 if(post_retirement==1){
 a4 ~ normal(0,0.05);
 lambda4 ~ normal(0,0.01);
 }
 c ~ normal(min(y),0.1);
 sigma ~ normal(0,1);
 }
 '
 
 # fit the model
 #rc_fit <- rstan::sampling(stanmodels$rc_flexible, data = mig_data, ...)
 rc_fit <- rstan::stan(model_code = rc_flexible, data = mig_data, ...)
 
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
 
 dfit <- tibble(age = x, 
                data = y, median = apply(y_hat, 2, median),
                lower = apply(y_hat, 2, quantile,0.025),
                upper = apply(y_hat, 2, quantile, 0.975),
                diff_sq = (!!sym("median") - !!sym("data"))^2)
 
 #TR: experimenting rm pipes re segfault error on osx...
 pars_df <- gather_draws(rc_fit, !!sym("a[0-9]\\[1\\]"),
                         !!sym("alpha[0-9]\\[1\\]"),
                         !!sym("mu[0-9]\\[1\\]"),
                         !!sym("lambda[0-9]\\[1\\]"),
                         !!sym("^c$"),
                         regex = TRUE) %>% 
            group_by(!!sym(".variable")) %>%
            summarise(median = median(!!sym(".value")),
                      lower = quantile(!!sym(".value"), 0.025),
                      upper = quantile(!!sym(".value"), 0.975)) %>% 
            dplyr::rename("variable" = !!sym(".variable")) %>% 
            mutate("variable" = gsub("\\[1\\]", "", "variable"))
 
 return(list(pars_df = pars_df, fit_df = dfit))
 
 # for sake of R CMD checks
 # .value <- .variable <- NULL
 # dt <- as.data.table(pars_df) 
 # dt <- 
 #   dt[, list(median = median( .value ),
 #             lower = quantile(.value, 0.025),
 #             upper = quantile(.value, 0.975)),
 #      by = list( .variable )] %>% 
 #   setnames(".variable","variable") %>% 
 #   as.tibble()
 
 return(list(pars_df = pars_df, fit_df = dfit))
}
