
# author: IW --------------------------------------------------------------

#' HMD pattern for adult ages.
#' @description Adjust rates in oldest ages using HMD pattern, based on log-quad method.
#' @details One possible scenario when mortality data on last ages is not reliable, is to use a mortality pattern with some known index in previous ages.
#' This function gives a HMD pattern based on 5q0 and 45q15, using log-quad model. Additionally, a value on mortality between 60 and 75 can be included to make a better adjustment in level. 
#' @param nMx numeric. Vector of mortality rates in abridged age classes.
#' @param Age integer. Single ages (abridged not allowed).
#' @param Sex character. Either male \code{"m"}, female \code{"f"}, or both \code{"b"}.
#' @param q0_5 numeric. Probability of death from born to age 5. By default implicit values in `nMx` should be entered.
#' @param q15_45 numeric. Probability of death from age 15 to age 60. By default implicit values in `nMx` should be entered.
#' @param fitted_logquad Optional, defaults to \code{NULL}. An object of class
#'  \code{wilmoth}. If full HMD is not enough, one
#'  can fit a Log-Quadratic (\url{https://github.com/mpascariu/MortalityEstimate}) model
#'  based on any other collection  of life tables;
#' @param q60_15 numeric. Probability of death from age 60 to age 75. When external information on those ages level is available, 
#' can be included to increase parameter `ax` from log-quad model in last ages (Li, 2003).
#' @param Age_transition integer. Form which age should transition to HMD pattern starts.
#' @param window_transition integer. Number of ages to the left and to the right of `Age_transition` to do a log-linear transition in rates.
#' @param plot_comparison Show or not a plot with the result.
#' @param ... Other arguments to be passed on to the \code{lt_single} function.
#' @export
#' @return life table as in \code{lt_single} function.
#' @examples
#' # Mortality rates from UN Chilean with e0=70. Wat would be the rates based on HMD pattern? 
#' # In this case making a transition of 10 years at age 80, and returning an OAG=100.
#' \dontrun{
#' lt <- DemoToolsData::modelLTx1
#' lt <- lt[lt$family == "Chilean" & lt$sex == "female" & lt$e0 == 70,]
#' chilean70_adjHMD <- HMD_old_logquad(nMx = lt$mx1, 
#'                                          Age = lt$age, 
#'                                        Sex = "f", 
#'                                        q0_5 = 1 - lt$lx1[lt$age==5]/lt$lx1[lt$age==0], 
#'                                        q15_45 = 1 - lt$lx1[lt$age==60]/lt$lx1[lt$age==15],
#'                                        Age_transition = 80, 
#'                                        window_transition = 10, 
#'                                        plot_comparison = TRUE,
#'                                        OAnew = 100)
#' # We know (as an example) that q60_15 is .5 higher than  what HMD pattern would be.
#' chilean70_adjHMD_augm <- HMD_old_logquad(nMx = lt$mx1, 
#'                                        Age = lt$age, 
#'                                        Sex = "f", 
#'                                        q0_5 = 1 - lt$lx1[lt$age==5]/lt$lx1[lt$age==0], 
#'                                        q15_45 = 1 - lt$lx1[lt$age==60]/lt$lx1[lt$age==15],
#'                                         q60_15 = (1 - lt$lx1[lt$age==75]/lt$lx1[lt$age==60]) * 1.5, 
#'                                         Age_transition = 80, window_transition = 10,
#'                                        OAnew = 100, plot_comparison = TRUE)
#' }

HMD_old_logquad <- function(nMx, Age = NULL, 
                                  Sex = "b",  
                                  q0_5 = NULL, q15_45 = NULL, 
                                  q60_15 = NULL,
                                  Age_transition = 80,
                                  window_transition = 3,
                                  plot_comparison = FALSE,
                                  fitted_logquad = NULL,
                                  ...){
     
     # check if age is not complete
     if(!is_single(Age)) stop("Need rates by single age.") 
     if(is.null(Age)) Age <- 0:length(nMx)
     
     # check if an optional fitted_logquad is specified
     if(is.null(fitted_logquad)){
       if(Sex == "b"){
         fitted_logquad <- DemoTools::fitted_logquad_b
       }
       if(Sex == "f"){
         fitted_logquad <- DemoTools::fitted_logquad_f
       }
       if(Sex == "m"){
         fitted_logquad <- DemoTools::fitted_logquad_m
       }
     }

     # load the log-quad model already fitted in the package (different from MortalityEstimate) and estimate with input data
     logquad_model <- lt_model_lq(Sex = Sex, q0_5 = q0_5, q15_45 = q15_45, fitted_logquad = fitted_logquad)
     mx_logquad_5q0_45q15 <- logquad_model$lt
     
     # augmented method if was asked
     if(!is.null(q60_15)){
          mx_logquad_5q0_45q15 <- tryCatch({
               logquad_augmented(coeffs = fitted_logquad$coefficients, 
                                 k = logquad_model$values$k, 
                                 Age = logquad_model$lt$Age,
                                 q0_5 = q0_5, q60_15 = q60_15, Sex = Sex)},
               error = function(e) {
                    warning("Augmented log-quad was not possible. Revise q60_15. Returned basic log-quad.")
                    mx_logquad_5q0_45q15})
     }
     
     # make it single
     mx_logquad_5q0_45q15 <- lt_abridged2single(nMx = mx_logquad_5q0_45q15$nMx, 
                                                Age = mx_logquad_5q0_45q15$Age, 
                                                lx = mx_logquad_5q0_45q15$lx, 
                                                Sex = Sex,
                                                ...)
     
     # smooth transition with a given length window
     Ages <- 0:max(mx_logquad_5q0_45q15$Age)
     Age_smooth <- (Age_transition-window_transition):(Age_transition+window_transition)
     Age_not_smooth <- Ages[!Ages %in% Age_smooth]
     nMx_to_smooth_transition <- data.frame(Age = Age_not_smooth,
                                            nMx = c(nMx[Age<min(Age_smooth)], 
                                                    mx_logquad_5q0_45q15$nMx[mx_logquad_5q0_45q15$Age>max(Age_smooth)]))
     nMx_interpolated <- exp(stats::approx(x = nMx_to_smooth_transition$Age, 
                                           y = log(nMx_to_smooth_transition$nMx), 
                                           xout = Age_smooth)$y)
     smooth_transtition_nMx <- dplyr::bind_rows(
          nMx_to_smooth_transition,
          data.frame(Age = Age_smooth, nMx = nMx_interpolated)) %>% 
          dplyr::arrange(Age)
     
     # plot diagnostic
     if(plot_comparison){
       df <- 
       rbind(
         data.frame(Age = Age, nMx = nMx, Type = "Input"),
         data.frame(Age = smooth_transtition_nMx$Age, nMx = smooth_transtition_nMx$nMx, Type = "Adjusted"),
         data.frame(Age = Age_smooth, nMx = nMx_interpolated, Type = "Transition")) %>% 
       ggplot2::ggplot(ggplot2::aes(x = Age, y = nMx, color = Type)) +
         ggplot2::geom_line() +
         ggplot2::geom_vline(xintercept = Age_transition, linetype = "dashed", color = "grey") +
         ggplot2::scale_y_log10() +
         ggplot2::theme_bw()
     }

     # rebuild lt and finish, with input OAnew if was not defined as additional input argument
     lt_extra_arguments <- list(...)
     if(!("OAnew" %in% names(lt_extra_arguments))){OAnew <- max(Age)} else {OAnew <- lt_extra_arguments$OAnew} 
     mx_logquad_5q0_45q15 <- lt_single_mx(nMx = smooth_transtition_nMx$nMx, 
                                          Age = smooth_transtition_nMx$Age, 
                                          Sex = Sex,
                                          OAnew = OAnew)
     return(mx_logquad_5q0_45q15)
}

#' Augmented logquad
#' @description Adjust rates in oldest ages that comes from a HMD model, using an external estimate of 15q60 (Li, 2014). As an example see\code{\link[DemoTools]{HMD_old_logquad}}.
#' @details Parameter \code{a(x)} is augmented based on an external estimate of 15q60. 
#' @param coeffs data.frame. Columns \code{a(x)}, \code{b(x)}, \code{c(x)} and \code{v(x)} from fitted logquad model. See \code{fitted_logquad_b}.
#' @param k numeric. Adult mortality related value from log-quad estimatation based on one or two input parameters. See \code{lt_model_lq}.
#' @param Sex character. Either male \code{"m"}, female \code{"f"}, or both \code{"b"}.
#' @param Age integer. Abridged lower bound ages. Same length than rows in `coeffs`.
#' @param q0_5 numeric. Probability of death from born to age 5.
#' @param q60_15 numeric. Probability of death from age 60 to age 75.
#' @param ... Other arguments to be passed on to the \code{lt_abridged} function.
#' @export
#' @return life table as in \code{lt_abridged} function.
#' @references See [Li (2014)](https://www.un.org/development/desa/pd/content/estimating-life-tables-developing-countries).

logquad_augmented <- function(coeffs, k, q0_5, q60_15, Sex = "b", Age, ...){
     
     # arrange parameters and get rates from the model
     params <- c(1, log(q0_5), log(q0_5)^2, k)
     mx_logquad <- exp(rowSums(as.matrix(coeffs) %*% diag(params)))
     
     # adjust ax with the ratio between input value and implicit value from model
     q60_15_hat <- q60_15  
     age_q60_15 <- which(Age == 60)
     q60_15 <- 1 - exp(-5*(sum(mx_logquad[c(age_q60_15, age_q60_15+1, age_q60_15+2)])))
     coeffs_hat <- coeffs
     coeffs_hat$ax[age_q60_15:nrow(coeffs_hat)] <- coeffs_hat$ax[age_q60_15:nrow(coeffs_hat)] + log(log(1-q60_15_hat)/log(1-q60_15)) 
     
     # new rates with changed ax
     mx_logquad_augm <- exp(rowSums(as.matrix(coeffs_hat) %*% diag(params)))
     lt_logquad_augm <- lt_abridged(nMx = mx_logquad_augm, Age = Age, Sex = Sex, ...)
     
     # output
     return(lt_logquad_augm)
}