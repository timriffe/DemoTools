# Author: ...
###############################################################################

#' Lee-Carter method with limited data.
#' 
#' @description Given a data frame with dates, sex and mortality data by age (rates, conditionated probabilities of death 
#' or survival function), this function interpolate/extrapolate life tables 
#' using the method for limited data suggested by  Li et. al (2004) (at least three observed years). 
#'
#' @details Based on spreedsheet "Li_2018_Limited_Lee-Carter-v4.xlsm" from UN. Useful for abridged or single ages, and allows
#' output in both formats also.
#' The main options are the use of non-divergent method for sex coherency (Li & Lee, 2005),
#' and the possibility of fitting `"k"` to replicate `"e_0"` at some given dates. 
#'
#' @note Draft Version
#'
#' @param data data.frame with cols: Date, Sex, Age, nMx (opt), nqx (opt), lx (opt)  
#' @param dates_out numeric. Vector of decimal years to interpolate or extrapolate.
#' @param Single logical. Wheter or not the lifetable output is by single ages.
#' @param dates_e0 numeric. Vector of decimal years where `"e_0"` should be fitted when apply method.
#' @param e0_Males numeric. Vector of life expectancy by year to be fitted. Same length than `"dates_e0"`.
#' @param e0_Females numeric. Vector of life expectancy by year to be fitted. Same length than `"dates_e0"`.
#' @param prev_divergence logical. Whether or not prevent divergence and sex crossover. Default `FALSE.`
#' @param OAG logical. Whether or not the last element of `nMx` (or `nqx` or `lx`) is an open age group. Default `TRUE.`
#' @param extrapLaw character. If extrapolating, which parametric mortality law should be invoked?
#' @param OAnew integer. Inferior limit of the open age group when extrapolate with `extrapLaw`.
#' @param OAnew logical.
#' @param verbose logical. Default `FALSE`.
#' @param ... Other arguments to be passed on to the \code{\link[lt_abridged]{lt_abridged}} function.
#' @seealso
#' \code{\link[DemoTools]{lt_abridged}}
#' @export
# TR: you can use markdown for this sort of thing, just getting used to it
#' @return Lifetable in a data.frame with columns
#' \itemize{
#'   \item{Date}{numeric. Dates included in dates_out},
#'   \item{Sex}{character. Male `"m"` or female `"f"`},
#'   \item{Age}{integer. Lower bound of abridged age class},
#'   \item{AgeInt}{integer. Age class widths.}
#'   \item{nMx}{numeric. Age-specific central death rates.}
#'   \item{nAx}{numeric. Average time spent in interval by those deceased in interval. }
#'   \item{nqx}{numeric. Age-specific conditional death probabilities.}
#'   \item{lx}{numeric. Lifetable survivorship}
#'   \item{ndx}{numeric. Lifetable deaths distribution.}
#'   \item{nLx}{numeric. Lifetable exposure.}
#'   \item{Sx}{numeric. Survivor ratios in uniform 5-year age groups.}
#'   \item{Tx}{numeric. Lifetable total years left to live above age x.}
#'   \item{ex}{numeric. Age-specific remaining life expectancy.}
#' }
#' @references
#' \insertRef{Li2005}{DemoTools}
#' \insertRef{Li2004}{DemoTools}
#' 
#' @examples
#' # Get data for Sweden (loaded from Human Mortality Database)
#' Age <- c(0,1,seq(5,100,5))
#' dates_in <- as.Date(c("1990-07-01", "2000-07-01", "2010-07-01"))
#' data <- data.frame(Date= c(rep(sort(rep(dates_in, length(Age))),2)),
#'                    Age = rep(Age,  2 * length(dates_in)),
#'                    Sex = c(rep("m", length(Age) * length(dates_in)),
#'                            rep("f", length(Age) * length(dates_in))),
#'                    nMx = c(0.006807, 0.000333, 0.000185, 0.000165, 0.000624, 
#'                            0.000859, 0.000955, 0.001081, 0.001375, 0.002098, 0.003206, 0.004996, 
#'                            0.009028, 0.014406, 0.023247, 0.039929, 0.067969, 0.110134, 0.180245, 
#'                            0.265287, 0.402896, 0.5, 0.004052, 0.000131, 9.8e-05, 0.000139, 
#'                            0.000473, 0.000761, 0.000733, 0.000759, 0.000984, 0.001526, 0.0026, 
#'                            0.004038, 0.00628, 0.010852, 0.018433, 0.031859, 0.0543, 0.093116, 
#'                            0.1638, 0.264152, 0.423773, 0.567665, 0.002712, 0.000159, 5e-05, 
#'                            9.9e-05, 0.00034, 0.000661, 0.000726, 0.000674, 0.000703, 0.001146, 
#'                            0.001939, 0.003303, 0.005405, 0.009007, 0.014105, 0.023629, 0.043337, 
#'                            0.077892, 0.141126, 0.241389, 0.379467, 0.678112, 0.005418, 0.000218, 
#'                            0.000131, 0.000165, 0.000268, 0.000325, 0.000433, 0.000536, 0.000776, 
#'                            0.001246, 0.002037, 0.002875, 0.005079, 0.007333, 0.012102, 0.020741, 
#'                            0.038245, 0.07144, 0.128579, 0.220927, 0.347505, 0.50933, 0.002837, 
#'                            0.000116, 9.3e-05, 0.000129, 0.000232, 0.000248, 0.000256, 0.000281, 
#'                            0.000655, 0.000936, 0.0017, 0.002723, 0.004497, 0.006826, 0.010464, 
#'                            0.017542, 0.031884, 0.061459, 0.115478, 0.204387, 0.335873, 0.484215, 
#'                            0.002402, 0.00014, 7.7e-05, 8e-05, 0.000196, 0.000248, 0.000257, 
#'                            0.000327, 0.000428, 0.000673, 0.001203, 0.002056, 0.003418, 0.005943, 
#'                            0.009141, 0.015038, 0.027323, 0.051479, 0.103552, 0.192314, 0.311358, 
#'                            0.479488))
#' dates_out <- as.Date(paste0(seq(1948,2018,5),"-07-01"))
#' 
#' # with abrev ages
#' lc_lim_data <- interp_lc_lim(data = data, dates_out = dates_out, OAG = F)
#' 
#' \dontrun{
#'   lc_lim_data %>% ggplot(aes(Age,nMx,col=factor(Date))) +
#'     geom_step() + scale_y_log10() + facet_wrap(~Sex)
#' }
#' 
#' # with simple ages as output
#' lc_lim_data_single <- interp_lc_lim(data = data, dates_out = dates_out, OAG = F,
#'                                     Single = TRUE)
#' 
#' \dontrun{
#'   lc_lim_data_single %>% ggplot(aes(Age,nMx,col=factor(Date))) +
#'     geom_step() + scale_y_log10() + facet_wrap(~Sex)
#' }
#' 
#' # Avoiding cross-over between sex
#' lc_lim_nondiv <- interp_lc_lim(data = data, dates_out = dates_out, OAG = F,
#'                                prev_divergence = TRUE)
#' \dontrun{
#'   lc_lim_nondiv %>% ggplot(aes(Age,nMx,col=factor(Date))) + 
#'     geom_step() + scale_y_log10() + facet_wrap(~Sex)
#' }
#' 
#' # Using information about e0 for past years
#' dates_e0 <- as.Date(paste0(seq(1960,2015,5),"-07-01"))   
#' e0_Males <- c(71.24, 71.74, 72.23, 72.17, 72.78, 73.78, 74.81, 76.18, 77.38, 78.42, 79.52, 80.32)
#' e0_Females <- c(74.88, 76.08, 77.21, 77.95, 78.85, 79.69, 80.4, 81.44, 82.02, 82.75, 83.47, 84.02)
#' lc_lim_fite0 <- interp_lc_lim(data = data, dates_out = dates_out, OAG = F,
#'                               dates_e0 = dates_e0,
#'                               e0_Males = e0_Males, 
#'                               e0_Females = e0_Females)
#' \dontrun{                           
#'   ggplot() + 
#'     geom_point(data = data.frame(Sex = c(rep("m",length(e0_Males)), rep("f",length(e0_Males))),
#'                                  ex = c(e0_Males, e0_Females),
#'                                  Date = rep(dec.date(dates_e0),2)),
#'                aes(Date,ex,col=factor(Sex)))+
#'     geom_line(data = lc_lim_fite0[lc_lim_fite0$Age==0,], aes(Date,ex,col=factor(Sex)))
#' }
#' 
#' # smooth and/or extend open age group, in this case input is for 80+
#' lc_lim_extOAg <- interp_lc_lim(data = data[data$Age<=80,], 
#'                                dates_out = dates_out, OAnew = 100,
#'                                OAG = F, extrapLaw = "kannisto")
#' \dontrun{ 
#'   ggplot() + 
#'     geom_step(data = lc_lim_extOAg, aes(Age,nMx,col=factor(Date))) +
#'     scale_y_log10() + facet_wrap(~Sex)
#' }
#' #End

interp_lc_lim <- function(data = NULL, # with cols: Date, Sex, Age, nMx (opt), nqx (opt), lx (opt)  
                          dates_out = dates_in, 
                          Single = FALSE,
                          dates_e0 = NULL,
                          e0_Males = NULL, 
                          e0_Females = NULL, 
                          prev_divergence = FALSE, 
                          OAnew = max(data$Age), 
                          OAG = TRUE, 
                          extrapLaw = NULL,
                          verbose = TRUE,
                          ...){
  
  # TR: capture args to filter out optional lifetable args to pass on in a named way?
  ExtraArgs <- as.list(match.call())
  
  # get always Mx -----------------------------------------------------------
  
  dates_in  <- unique(data$Date) %>% dec.date()
  dates_out <- dec.date(dates_out)
  
  # Two tries for dates_e0, otherwise we error
  # IW: the option for replicate e0 is optional
  if(!is.null(e0_Males)){
    if (is.null(dates_e0)){
      if (length(e0_Males) == length(dates_in)){
        dates_e0 <- dates_in
        if (verbose){
          cat("\ndates_e0 not specified, assuming:\n",paste(dates_in,collapse = ", "),"\n" )
        }
      }
    }
    if (is.null(dates_e0)){
      if (length(e0_Males) == length(dates_out)){
        dates_e0 <- dates_out
        if (verbose){
          cat("\ndates_e0 not specified, assuming:\n",paste(dates_out, collapse = ", "),"\n" )
        }
      }
    }
    if (is.null(dates_e0)){
      stop("\nSorry we can't guess the argument dates_e0, you'll need to specify it\n")
    }
  }
  
  
  # TR: note, we may as well offer to pass in args to lt_abridged(), which offers lots of control
  # over extrapolation.
  
  # axmethod = "pas", a0rule = "ak", Sex = "m", 
  # IMR = NA, region = "w", mod = TRUE, SRB = 1.05, OAG = TRUE, 
  # OAnew = max(Age), extrapLaw = "kannisto", extrapFrom = max(Age), 
  # extrapFit = Age[Age >= 60 & ifelse(OAG, Age < max(Age), TRUE)]
  # I added a capture of extra args entering the function all at the top.
  # instead of ... pass things in by name, above a list of default values for
  # lt_abridged(), most of which are the same for lt_single*()
  
  if (!any(names(data)%in%c("nMx", "nqx", "lx"))){
    stop("\nSorry we need some column called nMx, nqx or lx\n")
  }
  
  # cases for smooth older ages by default. And warn it (not sure if warning or message). Better cat?
  # IW: maybe you want more conditions, go ahead.
  # Extrap
  if(is.null(extrapLaw)){
    Ageext <- sort(unique(data$Age))
    extrapFrom <- max(Ageext)
    OAnew      <- 100
    if(max(Ageext) < 90){
      extrapLaw  <- "gompertz"
      if (verbose) message("A Gompertz function was fitted for older ages.")
      extrapFit = Ageext[Ageext >= 30 & ifelse(OAG, Ageext < max(Ageext), TRUE)]
    }else{
      extrapLaw  <- "kannisto"
      if (verbose) message("A Kannisto function was fitted for older ages.")
    }  
  }
  
  # here written out, not elegant. data.table() is preferred, 
  # but this is also fine
  . <- NULL
  data <- split(data, list(data$Sex, data$Date)) %>% 
    lapply(function(X){
      types     <- c("nMx","nqx","lx")
      lt_ambiguous_arg_types <- c("m","q","l")
      this_type <- lt_ambiguous_arg_types[types %in% colnames(X)]
      this_col  <- types[types %in% colnames(X)]
      this_sex  <- unique(X[["Sex"]])
      this_date <- unique(X[["Date"]])
      # TR: other args not passed in are scoped one level up
      out <- lt_ambiguous(x = X[[this_col]], 
                   Age = X[["Age"]], 
                   type = this_type,
                   Sex = this_sex,
                   Single = Single,
                   OAnew = OAnew, 
                   extrapLaw = extrapLaw,
                   ... = ...)
      out$Sex <- this_sex
      out$Date <- this_date
      out
    }) %>% 
    do.call("rbind", .)
    
  datadt <-
    data %>% 
    as.data.table()  
  
  # avoids 'no visible binding' warning
  .Sex <- NULL
  
  nMxf <-
    datadt %>% 
    subset(Sex == "f") %>% 
    data.table::dcast(Age ~ Date, value.var = "nMx") %>% 
    .[order(Age)]
  Age  <- nMxf[["Age"]]

  nMxf <- nMxf[, -1] %>% as.matrix()
  rownames(nMxf) <- Age
  
  nMxm <-
    datadt %>% 
    subset(Sex == "m") %>% 
    data.table::dcast(Age ~ Date, value.var = "nMx") %>% 
    .[order(Age)]
  nMxm <- nMxm[, -1] %>% as.matrix()
  rownames(nMxm) <- Age
  
  # LC at unequal intervals ---------------------------------------------------------
  
  #Age = sort(unique(data$Age)) # defined above for rownames
  ndates_in  <- length(dates_in)
  ndates_out <- length(dates_out)
  nAge       <- length(Age)
  
  # IW: make this modular
  # males
  lc_estimate_m <- lc_estimate(nMxm, dates_in, dates_out)
  axm <- lc_estimate_m[[1]]
  bxm <- lc_estimate_m[[2]]
  ktm <- lc_estimate_m[[3]]
  k0m <- lc_estimate_m[[4]]
  # females
  lc_estimate_f <- lc_estimate(nMxf, dates_in, dates_out)
  axf <- lc_estimate_m[[1]]
  bxf <- lc_estimate_m[[2]]
  ktf <- lc_estimate_m[[3]]
  k0f <- lc_estimate_m[[4]]
  
  # ask if prevent divergence and replicate target e0 ---------------------------------------------------------
  
  if (is.null(dates_e0)){ # not rep e0
    
    # basic
    nMxm_hat <- exp(axm + bxm %*% t(ktm))
    nMxf_hat <- exp(axf + bxf %*% t(ktf))
    
    # avoid divergence extrapolating to 1950: same bx and kt
    if (prev_divergence){
      kt = (ktm + ktf) * .5 # equal size male and female
      # error in vba code line 335
      # for (j in 1:ndates_out){
      #   for (i in 1:nAge){
      #     bxm[i] = (bxm[i] + bxf[i]) * .5
      #   }
      # }
      bx = (bxm + bxf) * .5
      k0 = (k0m + k0f) * .5
      nMxm_hat_div <- nMxm[,1] * exp(bx %*% t(kt-k0))
      nMxf_hat_div <- nMxf[,1] * exp(bx %*% t(kt-k0))
      
      # only for those years before min(dates_in). vba code explicit on that. 
      # IW: why not for dates_out>max(dates_in) also?
      dates_extrap <- dates_out < min(dates_in)
      nMxm_hat[,dates_extrap] <- nMxm_hat_div[,dates_extrap]
      nMxf_hat[,dates_extrap] <- nMxf_hat_div[,dates_extrap]
    }
  } else { # fit e0 at each target year
    stopifnot(length(e0_Males) == length(dates_e0))
    stopifnot(length(e0_Females) == length(dates_e0))
    # stepwise linear intra/extrapolation to target years. 
    # IW: Use interp(). Accepts matrix, so have to rbind and only get 1st row
    # source("R/AGEINT.R")
    e0m <- interp(rbind(e0_Males, 
                        e0_Males), 
                  dates_e0,
                  dates_out,
                  extrap = TRUE)[1, ]
    
    e0f <- interp(rbind(e0_Females, 
                        e0_Females),
                  dates_e0,
                  dates_out,
                  extrap = TRUE)[1, ]
    
    # avoid divergence: same bx but not kt.
    if (prev_divergence){
      bxm <- bxf <- (bxm + bxf) * .5
    }
    
    # Optimize kt for each LC extrap/interp and sex
    ktm_star = ktf_star = c()
    for (j in 1:ndates_out){
      ktm_star[j] <- optimize(f = lc_lim_kt_min,
                              interval = c(-20, 20),
                              ax = axm,
                              bx = bxm,
                              age = Age,
                              sex = "m",
                              e0_target = e0m[j])$minimum
      ktf_star[j] <- optimize(f = lc_lim_kt_min, # TR: add ...
                              interval = c(-20, 20),
                              ax = axf,
                              bx = bxf,
                              age = Age,
                              sex = "f",
                              e0_target = e0f[j])$minimum
    }
    
    # get rates with optim k.
    nMxm_hat <- exp(axm + bxm %*% t(ktm_star))
    nMxf_hat <- exp(axf + bxf %*% t(ktf_star))
  }
  
  # life tables output ------------------------------------------------------------
  
  colnames(nMxm_hat) <- dates_out
  colnames(nMxf_hat) <- dates_out
  . = NULL
  # TR: can use ... to pass in optional args. 
  Males_out <-
    lapply(colnames(nMxm_hat), function(x,MX,Age) {
      mx <- MX[, x]
      LT <- lt_ambiguous(x = mx, 
                            Age = Age, 
                            Sex = "m", 
                            Single = Single,
                            axmethod = "un", 
                            a0rule = "ak",
                            ... = ...) # quida de esto
      LT$Sex  <- "m"
      LT$Date <- as.numeric(x)
      LT
    }, MX = nMxm_hat, Age = Age) %>% 
    do.call("rbind", .)
  
  Females_out <-
    lapply(colnames(nMxf_hat), function(x,MX,Age) {
      mx <- MX[, x]
      LT <- lt_ambiguous(x = mx, 
                            Age = Age, 
                            Sex = "f", 
                            Single = Single,
                            axmethod = "un", 
                            a0rule = "ak",
                            ... = ...)
      LT$Sex  <- "f"
      LT$Date <- as.numeric(x)
      LT
    }, MX = nMxf_hat, Age = Age) %>% 
    do.call("rbind", .)
  
  out <- rbind(Males_out, Females_out)
  return(out)
  
}

# additional functions -------------------------------------------------------------------
# TR: new little helper to simplify some lines of code here and there.
interp_lc_lim_abk_m <- function(k,ax,bx){
  exp(ax + bx * k)
}
# optimize k for fitting e0
# TR: recommend lc_lim_kt_min(), following name conventions elsewhere in package
lc_lim_kt_min <- function(k,
                          ax,
                          bx,
                          age,
                          sex,
                          e0_target,
                          ...){
  Mx_hat <- interp_lc_lim_abk_m(k, ax, bx)
  e0 <- lt_ambiguous(x = Mx_hat, 
                        Age = age, 
                        Sex = sex,
                         ...)$ex[1]
  return(((e0-e0_target)/e0_target)^2)
}

# estimate lc functions for limited data
# IW: if wants normal LC (and data suitable for that), this is the fun to modify
lc_estimate <- function(M, dates_in, dates_out){
      ndates_in  <- length(dates_in)
      ax         <- rowSums(log(M))/ndates_in
      M_svd      <- svd(log(M)-ax)
      bx         <- M_svd$u[, 1]/sum(M_svd$u[, 1])
      kto        <- M_svd$d[1] * M_svd$v[, 1] * sum(M_svd$u[, 1])
      c          <- 0
      c[2]       <- (kto[ndates_in] - kto[1])/(dates_in[ndates_in] - dates_in[1])
      c[1]       <- kto[1] - c[2] * dates_in[1]
      kt         <- c[1] + c[2] * dates_out
      k0         <- c[1] + c[2] * dates_in[1]
      return(list(ax,bx,kt,k0))
}

# get lt for abrev/single ages and m/l/q input
# TR: this looks sufficiently ambiguous, love it, haha.
# I'm removing the else {} clauses, since this function should be pretty hermetic.
lt_ambiguous <- function(x = NULL, 
                         type = "m", # accepts"m", "q", or "l"
                         Age = NULL, 
                         Sex = NULL, 
                         Single = FALSE,
                         ...){
  if (type == "l"){
    x = lt_id_l_q(x)
    type = "q"
  }
  
  # a final catch
  out <- NULL
  # Abridged input lt
  if (is_abridged(Age)){
    
    # If we have nMx
    if (type == "m" & Single){
      out <- lt_abridged2single(nMx = x, Age = Age, Sex = Sex, ...)
    }
    if (type == "m" & !Single){
      out <- out <- lt_abridged(nMx = x, Age = Age, Sex = Sex, ...)  
    }
    # If we have nMx
    if (type == "q" & Single){
      out <- lt_abridged2single(nqx = x, Age = Age, Sex = Sex, ...)
    }
    if (type == "q" & !Single){
      out <- out <- lt_abridged(nqx = x, Age = Age, Sex = Sex, ...)  
    }
  }

  if (is_single(Age)){
    if (type == "m" & Single){
      out <- lt_single_mx(nMx = x, Age = Age, Sex = Sex, ...)
    }
    if (type == "m" & !Single){
      out <- lt_single_mx(nMx = x, Age = Age, Sex = Sex, ...)
      out <- lt_single2abridged(lx = out$lx,nLx = out$Lx, ex = out$ex) 
    }
    if (type == "q" & Single){
      out <- lt_single_qx(nqx = x, Age = Age, Sex = Sex, ...)
    }
    if (type == "q" & !Single){
      out <- lt_single_qx(qx = x, Age = Age, Sex = Sex, ...)
      out <- lt_single2abridged(lx = out$lx,nLx = out$Lx, ex = out$ex) 
    }
  }
    
  if (is.null(out)){
    # a final catch
    stop("please check function arguments")
  }  
  return(out)
}
