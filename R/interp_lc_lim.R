#' Lee-Carter method with limited data.
#' 
#' @description Given a data frame with dates, sex and mortality data by age (rates, conditionated probabilities of death 
#' or survival function), this function interpolate/extrapolate life tables 
#' using the method for limited data suggested by  Li et. al (2004) (at least three observed years). 
#'
#' @details Based on spreedsheet "Li_2018_Limited_Lee-Carter-v4.xlsm" from UN. 
#' Useful for abridged or single ages, and allows output in both formats also.
#' One option is the use of non-divergent method for sex coherency (Li & Lee, 2005).
#' The other is the possibility of fitting `"k"` to replicate `"e_0"` at some given dates. 
#'
#' @note Draft Version
#'
#' @param input data.frame with cols: Date, Sex, Age, nMx (opt), nqx (opt), lx (opt)  
#' @param dates_out numeric. Vector of decimal years to interpolate or extrapolate.
#' @param Single logical. Wheter or not the lifetable output is by single ages.
#' @param dates_e0 numeric. Vector of decimal years where `"e_0"` should be fitted when apply method.
#' @param e0_Males numeric. Vector of life expectancy by year to be fitted. Same length than `"dates_e0"`.
#' @param e0_Females numeric. Vector of life expectancy by year to be fitted. Same length than `"dates_e0"`.
#' @param prev_divergence logical. Whether or not prevent divergence and sex crossover. Default `FALSE.`
#' @param OAG logical. Whether or not the last element of `nMx` (or `nqx` or `lx`) is an open age group. Default `TRUE.`
#' @param verbose logical. Default `FALSE`.
#' @param SVD logical. Use Singular Value Decomposition for estimate b and k or Maximum Likelihood Estimation. Default `FALSE` for Maximum Likelihood Estimation.
#' @param ... Other arguments to be passed on to the \code{\link[DemoTools]{lt_abridged}} function.
#' @seealso
#' \code{\link[DemoTools]{lt_abridged}}
#' @export
# TR: you can use markdown for this sort of thing, just getting used to it
#' @return List with:
#' \itemize{
#'   \item Interpolated/extrapolated lifetables in a data.frame with columns:
#'   * `Date` numeric. Dates included in dates_out,
#'   * `Sex` character. Male `"m"` or female `"f"`,
#'   * `Age` integer. Lower bound of abridged age class,
#'   * `AgeInt`` integer. Age class widths.
#'   * `nMx` numeric. Age-specific central death rates.
#'   * `nAx` numeric. Average time spent in interval by those deceased in interval. 
#'   * `nqx` numeric. Age-specific conditional death probabilities.
#'   * `lx` numeric. Lifetable survivorship
#'   * `ndx` numeric. Lifetable deaths distribution.
#'   * `nLx` numeric. Lifetable exposure.
#'   * `Sx` numeric. Survivor ratios in uniform 5-year age groups.
#'   * `Tx` numeric. Lifetable total years left to live above age x.
#'   * `ex` numeric. Age-specific remaining life expectancy.
#'   \item List with estimated Lee-Carter parameters for each sex:
#'   * `kt` numeric time vector. Time trend in mortality level.
#'   * `ax` numeric age vector. Average time of `log(m_{x,t})`.
#'   * `bx` numeric age vector. Pattern of change in response to `kt`.
#' 
#' @references
#' \insertRef{Li2005}{DemoTools}
#' \insertRef{Li2004}{DemoTools}
#' 
#' @examples
#' # mortality rates from Sweden, for specific dates
#' data("mA_swe")
#' 
#' # needs mortality rates in this dates: 
#' dates_out <- as.Date(paste0(seq(1948,2018,5),"-07-01"))
#' 
#' # apply LC with limited data to extrap/interpolate
#' lc_lim_data <- interp_lc_lim(input = mA_swe, dates_out = dates_out, OAG = FALSE)
#' 
#' \dontrun{
#' lc_lim_data %>% ggplot(aes(Age,nMx,col=factor(round(Date,1)))) +
#'   geom_step() + scale_color_viridis_d() + 
#'   scale_y_log10() + theme_classic() + facet_wrap(~Sex)
#' }
#' 
#' # with simple ages as output
#' lc_lim_data_single <- interp_lc_lim(input = mA_swe, dates_out = dates_out, OAG = FALSE,
#'                                     Single = TRUE)
#' 
#' \dontrun{
#' lc_lim_data_single %>% ggplot(aes(Age,nMx,col=factor(round(Date,1)))) +
#'   geom_step() + scale_color_viridis_d() + 
#'   scale_y_log10() + theme_classic() + facet_wrap(~Sex)
#' }
#' 
#' # Avoiding cross-over between sex.
#' lc_lim_nondiv <- interp_lc_lim(input = mA_swe, dates_out = dates_out, OAG = FALSE,
#'                                prev_divergence = TRUE)
#' \dontrun{
#' lc_lim_nondiv %>% ggplot(aes(Age,nMx,col=factor(round(Date,1)))) +
#'   geom_step() + scale_color_viridis_d() + 
#'   scale_y_log10() + theme_classic() + facet_wrap(~Sex)
#' }
#' 
#' # Fitting information about e0 in Sweden for past years.
#' data("e0_swe")
#' lc_lim_fite0 <- interp_lc_lim(input = mA_swe, dates_out = dates_out, OAG = FALSE,
#'                               dates_e0 = unique(e0_swe$Date),
#'                               e0_Males = e0_swe$e0[e0_swe$Sex=="m"], 
#'                               e0_Females = e0_swe$e0[e0_swe$Sex=="f"])
#' \dontrun{                               
#' ggplot() + 
#'   geom_point(data = e0_swe, aes(Date,e0,col=factor(Sex)))+
#'   geom_line(data = lc_lim_fite0[lc_lim_fite0$Age==0,], aes(Date,ex,col=factor(Sex)))+
#'   labs(color = "Sex")+
#'   theme_classic()
#' }
#' 
#' # smooth and/or extend open age group, in this case input is for 80+, and chosen law is Makeham.
#' lc_lim_extOAg <- interp_lc_lim(input = mA_swe[mA_swe$Age<=80,], dates_out = dates_out,
#'                                OAG = FALSE,
#'                                OAnew=100,
#'                                extrapLaw = "makeham")
#' \dontrun{
#' ggplot() + 
#'   geom_step(data = lc_lim_extOAg, aes(Age,nMx,col=factor(round(Date,1)))) +
#'   scale_y_log10() + scale_color_viridis_d() + theme_classic() + facet_wrap(~Sex)
#'   }
#' #End

interp_lc_lim <- function(input = NULL,
                          dates_out = dates_in, 
                          Single = FALSE,
                          dates_e0 = NULL,
                          e0_Males = NULL, 
                          e0_Females = NULL, 
                          prev_divergence = FALSE, 
                          OAG = TRUE,
                          verbose = TRUE,
                          SVD = FALSE,
                          ...){
  
  # TR: ExtraArgs has a problem in that it won't capture & pass NULL defaults
  # mget(names(formals()),sys.frame(sys.nframe()))
  # IW: this captures everything: dots & NULLS
  ExtraArgs = c(as.list(environment()), list(...))
  ExtraArgs = ExtraArgs[! names(ExtraArgs) %in% c("input","Single")]

  # dates
  dates_in  <- unique(input$Date) %>% dec.date()
  dates_out <- dec.date(dates_out)
  
  # take care if lc with limited data is suitable
  if (length(dates_in)<3){
    stop("\nYou need three observed dates at least.")
  }
  if (all(diff(dates_in)==1)){
    stop("\nYou have single-year-interval data and probably should use basic Lee-Carter method.")
  }
  
  # Two tries for dates_e0, otherwise we error
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
  
  if (!any(names(input)%in%c("nMx", "nqx", "lx"))){
    stop("\nSorry we need some column called nMx, nqx or lx\n")
  }
  
  # get always Mx -----------------------------------------------------------
  
  # TR data.table() is preferred, 
  # and maybe this function shouldn't be anonymous, but rather called inside
  # inputdt[, new_function(.SD,...),by=list(Sex,Date)] # or similar.
  . <- NULL
  inputdt <- split(input, list(input$Sex, input$Date)) %>% 
              lapply(
                function(X) do.call(lt_smooth_ambiguous,
                                    c(list(input=X), ExtraArgs))) %>% 
              do.call("rbind", .)%>% 
              as.data.table() 

  # avoids 'no visible binding' warning
  Sex <- NULL
  
  nMxf <-
    inputdt %>% 
    subset(Sex == "f") %>% 
    data.table::dcast(Age ~ Date, value.var = "nMx") %>% 
    .[order(Age)]
  Age  <- nMxf[["Age"]]

  nMxf <- nMxf[, -1] %>% as.matrix()
  rownames(nMxf) <- Age
  
  nMxm <-
    inputdt %>% 
    subset(Sex == "m") %>% 
    data.table::dcast(Age ~ Date, value.var = "nMx") %>% 
    .[order(Age)]
  nMxm <- nMxm[, -1] %>% as.matrix()
  rownames(nMxm) <- Age
  
  # LC at unequal intervals ---------------------------------------------------------
  
  #Age = sort(unique(input$Age)) # defined above for rownames
  ndates_in  <- length(dates_in)
  ndates_out <- length(dates_out)
  nAge       <- length(Age)
  
  # males
  lc_estimate_m <- interp_lc_lim_estimate(nMxm, dates_in, dates_out, SVD)
  axm <- lc_estimate_m[[1]]
  bxm <- lc_estimate_m[[2]]
  ktm <- lc_estimate_m[[3]]
  k0m <- lc_estimate_m[[4]]
  # females
  lc_estimate_f <- interp_lc_lim_estimate(nMxf, dates_in, dates_out, SVD)
  axf <- lc_estimate_f[[1]]
  bxf <- lc_estimate_f[[2]]
  ktf <- lc_estimate_f[[3]]
  k0f <- lc_estimate_f[[4]]
  
  # ask if prevent divergence and replicate target e0 ---------------------------------------------------------
  
  if (is.null(dates_e0)){ # not rep e0
    
    # basic
    nMxm_hat <- exp(axm + bxm %*% t(ktm))
    nMxf_hat <- exp(axf + bxf %*% t(ktf))
    
    # avoid divergence extrapolating
    if (prev_divergence){
      kt = (ktm + ktf) * .5 # equal size male and female
      bx = (bxm + bxf) * .5 # # error in vba code line 335. A parameter controls that only for reproducing purpose
      k0 = (k0m + k0f) * .5
      
      # apply common factor to rates with already specific factor (formula 6 in Li (2005)), 
      # not like Li´s paper way in this case. IW: we can improve this if UN wants.
      nMxm_hat_div <- nMxm[,1] * exp(bx %*% t(kt-k0))
      nMxf_hat_div <- nMxf[,1] * exp(bx %*% t(kt-k0))
      
      # only for those years before min(dates_in). UN code explicit on that. 
      # IW: why not for dates_out>max(dates_in) also? Ask UN.
      dates_extrap <- dates_out < min(dates_in)
      nMxm_hat[,dates_extrap] <- nMxm_hat_div[,dates_extrap]
      nMxf_hat[,dates_extrap] <- nMxf_hat_div[,dates_extrap]
    }
  } else { # fit e0 at each target year
    stopifnot(length(e0_Males) == length(dates_e0))
    stopifnot(length(e0_Females) == length(dates_e0))
    # stepwise linear intra/extrapolation to target years. 
    # IW: Use interp(). Accepts matrix, so have to rbind and only get 1st row
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
      ktm_star[j] <- optimize(f = interp_lc_lim_kt_min,
                              interval = c(-20, 20),
                              ax = axm,
                              bx = bxm,
                              age = Age,
                              sex = "m",
                              e0_target = e0m[j],
                              ...)$minimum
      ktf_star[j] <- optimize(f = interp_lc_lim_kt_min, # TR: add ...
                              interval = c(-20, 20),
                              ax = axf,
                              bx = bxf,
                              age = Age,
                              sex = "f",
                              e0_target = e0f[j],
                              ...)$minimum
    }
    
    # get rates with optim k.
    nMxm_hat <- exp(axm + bxm %*% t(ktm_star))
    nMxf_hat <- exp(axf + bxf %*% t(ktf_star))
    ktm <- ktm_star
    ktf <- ktf_star
  }
  
  # life tables output ------------------------------------------------------------
  
  colnames(nMxm_hat) <- dates_out
  colnames(nMxf_hat) <- dates_out
  . = NULL
  
  Males_out <-
    lapply(colnames(nMxm_hat), function(x,MX,Age) {

      mx <- MX[, x]
      LT <- lt_ambiguous(x = mx,
                            type = "m",
                            Age = Age, 
                            Sex = "m",
                            Single = Single,
                            ...)
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
                            ...)
      LT$Sex  <- "f"
      LT$Date <- as.numeric(x)
      LT
    }, MX = nMxf_hat, Age = Age) %>% 
    do.call("rbind", .)
  lt_hat <- rbind(Males_out, Females_out)
  
  # for output
  lc_params <- list(ax = data.frame(Male = axm, Female = axf),
                    bx = data.frame(Male = bxm, Female = bxf),
                    kt = data.frame(Male = ktm, Female = ktf))
  return(list(lt_hat = lt_out,
              lc_params = lc_params)
         )
}

#'  Optimize k
#' @description Optimize estimated k from LC with limited data model, 
#' for fitting given e_0 at same dates
#' @details Given LC parameters at some date, change a bit k for replicate already know e_0 values. 
#' This is useful to give some sort of flexibility, and not follow strictly linear model implied in LC model,
#' but taking advantage of estimated structure (ax) and change by age (bx) for some trustable period.
#' @param k numeric. k parameter from LC model.
#' @param ax numeric. Vector (same length of age) of parameters from LC model.
#' @param bx numeric. Vector (same length of age) of parameters from LC model.
#' @param age numeric. 
#' @param sex numeric.
#' @param e0_target numeric.
#' @param ... Other arguments to be passed on to the \code{\link[DemoTools]{lt_abridged}} function.
#' @export
interp_lc_lim_kt_min <- function(k,
                          ax,
                          bx,
                          age,
                          sex,
                          e0_target,
                          ...){
  Mx_hat <- as.numeric(interp_lc_lim_abk_m(k, ax, bx))
  e0 <- lt_ambiguous(x = Mx_hat, 
                        Age = age, 
                        Sex = sex,
                         ...)$ex[1]
  return(((e0-e0_target)/e0_target)^2)
}

#' wrapper fun for `"interp_lc_lim_estimate"` function
#' @description wrapper fun to estimate rates from LC parameters
#' @inheritParams interp_lc_lim_kt_min 
#' @export
interp_lc_lim_abk_m <- function(k,ax,bx){
  exp(ax + bx * k)
}

# estimate LC for limited data
#' Estimate LC with limited data params
#' @description Estimate LC with limited data from a matrix of rates (age by dates).
#' @details SVD for ax and bx. Fit a simmple linear model for k and interp/extrapolate for objective dates.
#' @param M numeric. Matrix with many rows as ages and columns as dates_in.
#' @param dates_in numeric. Vector of dates with input rates.
#' @param dates_out numeric. Vector of dates for estimate a set of rates.
#' @param SVD logical. Use Singular Value Decomposition for estimate b and k or Maximum Likelihood Estimation. Default `FALSE` for Maximum Likelihood Estimation.
#' @references
#' \insertRef{Li2004}{DemoTools}
#' @export
interp_lc_lim_estimate <- function(M, dates_in, dates_out, SVD = F){
      ndates_in  <- length(dates_in)
      ax         <- rowSums(log(M))/ndates_in
      if(SVD==TRUE){
        # Singular Value Decomposition
        M_svd      <- svd(log(M)-ax)
        bx         <- M_svd$u[, 1]/sum(M_svd$u[, 1])
        kto        <- M_svd$d[1] * M_svd$v[, 1] * sum(M_svd$u[, 1])  
      }else{
        # likelihood method
        kto  <- colSums(log(M))-sum(ax) # because sum(bx)==1
        bx   <- (log(M) - ax)%*%kto/sum(kto^2)
      }
      c          <- 0
      c[2]       <- (kto[ndates_in] - kto[1])/(dates_in[ndates_in] - dates_in[1])
      c[1]       <- kto[1] - c[2] * dates_in[1]
      # explanation ratio
      R = 1- sum((log(M) - (ax + bx %*% t(kto)))^2)/sum((log(M)-ax)^2)
      # extrapolated k
      kt         <- c[1] + c[2] * dates_out
      # initial k (useful for avoiding divegence case)
      k0         <- c[1] + c[2] * dates_in[1]
      return(list(ax=ax,bx=bx,kt=kt,k0=k0,R=R))
}

# smooth rule previous to solve ambiguous
#' Smooth and apply lt_ambiguous
#' @description Considering different mortality input for each sex/year data, 
#' smooth olders with gompertz or kannisto in case no law was specified, 
#' and return a data.frame with standard LT. 
#' @details Gompertz is chosen if last age is less than 90. Else Kannisto. 
#' @param input data.frame. with cols: Date, Sex, Age, nMx (opt), nqx (opt), lx (opt)  
#' @param ... Other arguments to be passed on to the \code{\link[DemoTools]{lt_abridged}} function.
#' @export
lt_smooth_ambiguous <- function(input, ...){
  
  ExtraArgs = c(as.list(environment()), list(...))
  ExtraArgs = ExtraArgs[! names(ExtraArgs) %in% c("input")]

  # get only cols with values (could entry nMx for some sex/year and nqx or lx for other)
  X <- input[!sapply(input, function(x) all(is.na(x)))]
  
  types     <- c("nMx","nqx","lx")
  this_type <- types[types %in% colnames(X)]
  this_sex  <- unique(X[["Sex"]])
  this_date <- unique(X[["Date"]])
  
  # cases for smooth older ages by default
  if(!"extrapLaw" %in% names(ExtraArgs)){  
    Ageext <- sort(unique(X$Age))
    this_extrapFrom <- max(Ageext)
    this_OAnew = 100
    if(this_extrapFrom < 90){
      this_extrapLaw  <- "gompertz"
      if (ExtraArgs$verbose) cat(paste0("A Gompertz function was fitted for older ages for sex ",
                                  this_sex, " and date ",this_date,".\n"))
      # TR: changed this. 30 could be sort of low in some situations.
      this_extrapFit = Ageext[Ageext >= (this_extrapFrom - 30) & ifelse(ExtraArgs$OAG, Ageext < max(Ageext), TRUE)]
    }else{
      this_extrapLaw  <- "kannisto"
      if (ExtraArgs$verbose) cat(paste0("A Kannisto function was fitted for older ages for sex ",
                                  this_sex, " and date ",this_date,".\n"))
      this_extrapFit = Ageext[Ageext >= 60 & ifelse(ExtraArgs$OAG, Ageext < max(Ageext), TRUE)]
    }
    # TR: other args not passed in are scoped one level up
    thisExtraArgs <- ExtraArgs[!names(ExtraArgs) %in% 
                                 c("extrapLaw","extrapFit","extrapFrom","OAnew","Single")]
    out <- lt_ambiguous(x = X[[this_type]], 
                        Age = X[["Age"]], 
                        type = this_type,
                        Sex = this_sex,
                        extrapLaw = this_extrapLaw,
                        extrapFit = this_extrapFit,
                        extrapFrom = this_extrapFrom, 
                        OAnew = this_OAnew,
                        ...)
  }else{
    out <- lt_ambiguous(x = X[[this_type]], 
                        Age = X[["Age"]], 
                        type = this_type,
                        Sex = this_sex,
                        ...)  
  }
  
  out$Sex <- this_sex
  out$Date <- this_date
  out
}
