#===============================================================================
# Wrap model life tables methods
# Ilya Kashnitsky, ilya.kashnitsky@gmail.com
#===============================================================================

# Log-Quad family (Wilmoth) -----------------------------------------------
# Author: Marius D. Pascariu
# License: GNU General Public License v3.0
# Last update: Tue Dec  4 22:17:06 2018
# https://github.com/mpascariu/MortalityEstimate/blob/master/R/fun_Wilmoth.R


#' Estimate Wilmoth Model Life Table
#' 
#' Construct model life tables based on the Log-Quadratic (wilmoth) estimates
#' with various choices of 2 input parameters: 
#' \code{q0_5, q0_1, q15_45, q15_35} and \code{e0}. There are 8 possible 
#' combinations (see examples below). 
#' 
#' @details Due to limitations of the R language the notation for probability 
#' of dying \code{nqx} is written \code{qx_n}, where \code{x} and \code{n} are 
#' integers. For example \code{45q15} is represented as \code{q45_15}.
#' @note This function is ported from \code{MortalityEstimate::wilmothLT} experimental package by Marius Pascariu. The package is no longe maintained. The latest version can be found here: \url{https://github.com/mpascariu/MortalityEstimate}
#' @param Sex Choose the sex of the population. This choice defines the use
#' of a corresponding Log-Quadratic (\code{wilmoth})
#'  model fitted for the whole Human Mortality Database (as of Dec 2019, 
#'  there are 968 life tables for each sex).
#' The following options are available: \itemize{
#'   \item{\code{"b"}} -- Both sex;
#'   \item{\code{"f"}} -- Females;
#'   \item{\code{"m"}} -- Males.
#'   }
#' @param fitted_logquad Optional, defaults to \code{NULL}. An object of class
#'  \code{wilmoth}. If full HMD is not enough, one 
#'  can fit a Log-Quadratic (\url{https://github.com/mpascariu/MortalityEstimate}) model 
#'  based on any other collection  of life tables;
#' @param q0_5 5q0. The probability that a new-born will die during the 
#' subsequent 5 years;
#' @param q0_1 1q0. The probability that a life aged 0 will die during the 
#' following year;
#' @param q15_45 45q15. The probability that a life aged 15 will die during 
#' the subsequent 45 years;
#' @param q15_35 35q15. The probability that a life aged 15 will die during 
#' the subsequent 35 years;
#' @param e0 Life expectancy at birth;

#' @param radix Life table radix. Default: 10^5;
#' @param tol Tolerance level for convergence. The tolerance level, is relevant 
#' for case 7 and 8 (e0 and 45q15 or 35q15 are known);
#' @param maxit Maximum number of iterations allowed. Default: 100;
#' @param ... Additional arguments affecting the predictions produced.
#' @return The output is of class \code{lt_model_lq} with the components:
#'  \item{lt}{ Life table matching given inputs}
#'  \item{values}{ Associated values of \code{q0_5, q0_1, q15_45, q15_35} 
#' and \code{e0}.}
#' @importFrom stats uniroot
#' @examples 
#' 
#' # Build life tables with various choices of 2 input parameters
#'  
#' # case 1: Using 5q0 and e0
#' L1 <- lt_model_lq(Sex = "b", q0_5 = 0.05, e0 = 65)
#' L1
#' ls(L1)
#' 
#' L1f <- lt_model_lq(Sex = "f", q0_5 = 0.05, e0 = 65)
#' L1m <- lt_model_lq(Sex = "m", q0_5 = 0.05, e0 = 65)
#' 
#' # case 2: Using 5q0 and 45q15
#' L2 <- lt_model_lq(Sex = "b", q0_5 = 0.05, q15_45 = 0.2)
#' 
#' # case 3: Using 5q0 and 35q15
#' L3 <- lt_model_lq(Sex = "b", q0_5 = 0.05, q15_35 = 0.125)
#' 
#' # case 4: Using 1q0 and e0
#' L4 <- lt_model_lq(Sex = "b", q0_1 = 0.01, e0 = 65)
#' 
#' # case 5: Using 1q0 and 45q15
#' L5 <- lt_model_lq(Sex = "b", q0_1 = 0.05, q15_45 = 0.2)
#' 
#' # case 6: Using 1q0 and 35q15
#' L6 <- lt_model_lq(Sex = "b", q0_1 = 0.05, q15_35 = 0.125)

#' 
#' # case 7: Using 45q15 and e0
#' L7 <- lt_model_lq(Sex = "b", q15_45 = 0.125, e0 = 65)
#' 
#' # case 8: Using 35q15 and e0
#' L8 <- lt_model_lq(Sex = "b", q15_35 = 0.15, e0 = 65)
#' 
#' @export
lt_model_lq <- function(
    Sex = c("b", "f", "m"), # has to be specified always
    fitted_logquad = NULL, 
    q0_5 = NULL, 
    q0_1 = NULL, 
    q15_45 = NULL, 
    q15_35 = NULL, 
    e0 = NULL, 
    
    radix = 1e5, 
    tol = 1e-9, 
    maxit = 200, 
    ...
) {
    # securing Sex input
    Sex <- tolower(Sex)
    if(Sex%in%c("total", "t")){
        Sex <- "b"
    }
    if(!nchar(Sex)==1){
        Sex <- substr(Sex,1,2)
    }
    
    
    # check if an optional fitted_logquad is specified
    
    if(is.null(fitted_logquad)){
        
        if(Sex == "b"){
            fitted_logquad <- fitted_logquad_b
        }
        if(Sex == "f"){
            fitted_logquad <- fitted_logquad_f
        }
        if(Sex == "m"){
            fitted_logquad <- fitted_logquad_m
        }
    }
    
    
    my_case <- find.my.case(q0_5, q0_1, q15_45, q15_35, e0)
    cf      <- coef(fitted_logquad)
    x       <- fitted_logquad$input$x
    
    # Cases 1-3:  5q0 is known, plus e0, 45q15 or 45q15

    if (my_case %in% c("C1", "C2", "C3")) {
        if (my_case == "C1") fun.k <- function(k) {
            lthat.logquad(cf, x, q0_5, k, radix, Sex = Sex)$lt$ex[1] - e0
        }
        if (my_case == "C2") fun.k <- function(k) { 
            lt <- lthat.logquad(cf, x, q0_5, k, radix, Sex = Sex)$lt
            return(1 - lt[lt$Age == 60, "lx"] / lt[lt$Age == 15, "lx"] - q15_45)
        }
        if (my_case == "C3") fun.k <- function(k) { 
            lt <- lthat.logquad(cf, x, q0_5, k, radix, Sex = Sex)$lt
            return(1 - lt[lt$Age == 50, "lx"] / lt[lt$Age == 15, "lx"] - q15_35)
        }
        
        root <- uniroot(f = fun.k, interval = c(-10, 10))$root
        tmp  <- lthat.logquad(cf, x, q0_5, k = root, radix, Sex = Sex) 
    }
    
    # Cases 4-6: 1q0 is known, plus e0, 45q15 or 35q15;
    # after finding 5q0 (assume k=0, but it doesn't matter), these become Cases 1-3
    if (my_case %in% c("C4","C5","C6") ) {
        fun.q0_5 <- function(q0_5) lthat.logquad(cf, x, q0_5, k = 0, radix, Sex = Sex)$lt$nqx[1] - q0_1
        root <- uniroot(f = fun.q0_5, interval = c(1e-5, 0.8))$root
    }

    if (my_case == "C4") tmp <- lt_model_lq(fitted_logquad = fitted_logquad, q0_5 = root, e0 = e0, ...)
    if (my_case == "C5") tmp <- lt_model_lq(fitted_logquad = fitted_logquad, q0_5 = root, q15_45 = q15_45, ...)
    if (my_case == "C6") tmp <- lt_model_lq(fitted_logquad = fitted_logquad, q0_5 = root, q15_35 = q15_35, ...)
    
    # Case 7 and 8: e0 and 45q15 or 35q15 are known; must find both 5q0 and k
    if (my_case %in% c("C7", "C8")) {
        k    <- q0_5 <- 0
        iter <- crit <- 1
        
        while (crit > tol & iter <= maxit) {
            k.old  <- k
            q0_5.old <- q0_5
            # Get new 5q0 from e0 assuming k (case 9 from MortalityEstimate::wilmothLT)
            fun.q0_5 = function(q0_5) { 
                lthat.logquad(cf, x, q0_5, k, 
                              q0_1 = NULL, q15_45 =  NULL, q15_35 = NULL, 
                              radix, Sex = Sex)$lt$ex[1] - e0 
            }
            root <- uniroot(f = fun.q0_5, interval = c(1e-4, 0.8))$root
            q0_5 <- lthat.logquad(
                cf, x, q0_5 = root, k, 
                q0_1 = NULL, q15_45 =  NULL, q15_35 = NULL,
                radix, Sex = Sex
            )$values$q0_5 
            # Get k from 45q15 or 35q15 assuming 5q0
            if (my_case == "C7") tmp = lt_model_lq(fitted_logquad = fitted_logquad, q0_5 = q0_5, q15_45 = q15_45, ...)
            if (my_case == "C8") tmp = lt_model_lq(fitted_logquad = fitted_logquad, q0_5 = q0_5, q15_35 = q15_35, ...)
            k  <- tmp$values$k
            crit <- sum(abs(c(k, q0_5) - c(k.old, q0_5.old)))
            iter <- iter + 1
        }
        if (iter > maxit) {
            warning("number of iterations reached maximum without convergence", 
                    call. = FALSE)
        }
    }
    
    # Return life table plus values of the 6 possible inputs
    out = list(lt = tmp$lt, 
               values = tmp$values)
    out = structure(class = "lt_model_lq", out)
    return(out)
}


#' Estimated life table using the log-quadratic model
#' 
#' @param coefs Estimated coefficients
#' @inheritParams lt_model_lq
#' @keywords internal
#' @export
lthat.logquad <- function(coefs, 
                          x, 
                          q0_5, 
                          k, 
                          radix,
                          Sex, # needed to pass over to lt_id_morq_a
                          ...) {
    
    # Sex <- Sex
    h     <- log(q0_5)
    mx    <- with(as.list(coefs), exp(ax + bx*h + cx*h^2 + vx*k))
    # estimate ax
    age_int <- c(diff(x), Inf)
    ax    <- lt_id_morq_a(
        nMx = mx, 
        Age = x,
        AgeInt = age_int,
        axmethod = "un",
        Sex = Sex,
        region = "w"
    )
    # qx from mx and estimated ax
    qx <- lt_id_ma_q(nMx = mx, nax = ax, AgeInt = age_int, IMR = NA)
    # Force 4q1 (and thus 4m1) to be consistent with 1q0 and 5q0
    qx[2] <- 1 - (1 - q0_5)/(1 - qx[1])
    mx[2] <- lt_id_qa_m(nqx = qx, nax = ax)[2]
    names(mx) = names(qx) <- rownames(coefs)
    
    LT     <- lt_abridged(Age = x, nMx = mx, lx0 = radix)
    e0     <- LT$ex[1]
    q0_1   <- LT$nqx[1]
    q15_45 <- 1 - LT[LT$Age == 60, "lx"] / LT[LT$Age == 15, "lx"]
    q15_35 <- 1 - LT[LT$Age == 50, "lx"] / LT[LT$Age == 15, "lx"]
    values <- data.frame(k, q0_1, q0_5, q15_35, q15_45, e0, row.names = "")
    
    # Exit
    out <- list(lt = LT, values = values)
    return(out)
}


#' Function that determines the case/problem we have to solve
#' It also performs some checks
#' @inheritParams lt_model_lq
#' @keywords internal
find.my.case <- function(q0_5, q0_1, q15_45, q15_35, e0) {
    # Test that at least of one of 1q0 and 5q0 is null
    input   <- as.list(environment())
    my_case <- unlist(lapply(input, is.null))
    
    if (sum(my_case[c(1, 2)]) == 0) {
        stop("cannot have both 'q0_1' and 'q0_5' as inputs", call. = FALSE)
    }
    
    # Test that at least of one of 45q15 and 35q15 is null
    if (sum(my_case[c(3, 4)]) == 0) {
        stop("cannot have both 'q15_45' and 'q15_35' as inputs", call. = FALSE)
    }
    
    # Test that exactly two inputs are non-null
    if (sum(my_case) != 3) {
        stop("must have exactly two inputs", call. = FALSE)
    }
    
    # There are 8 cases:  "5 choose 2" = 10, but we disallow two cases 
    # (1q0 and 5q0, or 45q15 and 35q15)

    if (all(my_case == c(F,T,T,T,F))) case = "C1"
    if (all(my_case == c(F,T,F,T,T))) case = "C2" 
    if (all(my_case == c(F,T,T,F,T))) case = "C3"

    if (all(my_case == c(T,F,T,T,F))) case = "C4"   
    if (all(my_case == c(T,F,F,T,T))) case = "C5"
    if (all(my_case == c(T,F,T,F,T))) case = "C6" 

    if (all(my_case == c(T,T,F,T,F))) case = "C7" 
    if (all(my_case == c(T,T,T,F,F))) case = "C8" 
    return(case)
}
