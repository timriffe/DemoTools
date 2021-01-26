# Author: ...
###############################################################################

#' Lee-Carter method with limited data.
#' 
#' @description Given a matrix of age (rows) by year (cols) for each sex, this function interpolate/extrapolate life tables 
#' using the method for limited data suggested by  Li et. al (2004) (at least three observed years). 
#'
#' @details Based on spreedsheet "Li_2018_Limited_Lee-Carter-v4.xlsm" from UN. Only useful for abridged ages.
#' The main options are the use of non-divergent method for sex coherency (Li & Lee, 2005),
#' and the possibility of fitting `"k"` to replicate `"e_0"` at some given years. 
#' Coale-Guo method allows to smooth rates at older ages.
#'
#' @note Draft Version
#'
#' @param data data.frame with cols: Date, Sex, Age, nMx (opt), nqx (opt), lx (opt)  
#' @param dates_out numeric. Vector of decimal years to interpolate or extrapolate.
#' @param ... Arguments passed to `\link{lt_abridged}`.
#' @param dates_e0 numeric. Vector of decimal years where `"e_0"` should be fitted when apply method.
#' @param e0_Males numeric. Vector of life expectancy by year to be fitted. Same length than `"dates_e0"`.
#' @param e0_Females numeric. Vector of life expectancy by year to be fitted. Same length than `"dates_e0"`.
#' @param prev_divergence logical. Whether or not prevent divergence and sex crossover. Default `FALSE.`
#' @param OAG logical. Whether or not the last element of `nMx` (or `nqx` or `lx`) is an open age group. Default `TRUE.`
#' @param extrapLaw character. If extrapolating, which parametric mortality law should be invoked? Options include  `"Coale-Guo"` or other functions in `"lt_abridged"` function.
#' @export
# TR: you can use markdown for this sort of thing, just getting used to it
#' @return Lifetable in a data.frame with columns
#' \itemize{
#'   \item{Year}{integer. Years included in dates_out},
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
# TR: Example commented out because we don't know dimension of first data objects, as Age and dates_in
# were not previously defined.
#' @examples
#' # Using Lee-Carter method for Sweden, assuming only is available mortality rates
#' for years 1980, 1990, 2000 and 2010.
#' # LC with limited data
#' dates_in <- seq(1980,2010,10)
#' dates_out <- seq(1948,2018,5)
#' data <- readHMDweb("SWE", "Mx_5x1", user, pass, fixup = TRUE)%>% 
#'   select(Date = Year, Age, Male, Female) %>% 
#'   gather(Sex,nMx,-Date,-Age) %>% 
#'   mutate(Sex=ifelse(Sex=="Male","m","f")) %>% 
#'   filter(Date %in% dates_in, Age<=100)
#' 
#' # abr ages
#' lc_lim_data <- interp_lc_lim(data = data, dates_out = dates_out)
#' 
#' lc_lim_data %>% ggplot(aes(Age,nMx,col=factor(Date))) +
#'   geom_step() + scale_y_log10() + facet_wrap(~Sex)
#' 
#' # simple ages --- takes a loooong time, not performanced yet
#' lc_lim_data <- interp_lc_lim(data = data, dates_out = dates_out, 
#'                              Single = T, OAnew = 100)
#' 
#' lc_lim_data %>% ggplot(aes(Age,nMx,col=factor(Date))) +
#'   geom_step() + scale_y_log10() + facet_wrap(~Sex)
#' #End

interp_lc_lim <- function(data = NULL, # with cols: Date, Sex, Age, nMx (opt), nqx (opt), lx (opt)  
                          dates_out = dates_in, 
                          Single = F,
                          dates_e0 = NULL,
                          e0_Males = NULL, 
                          e0_Females = NULL, 
                          prev_divergence = FALSE, 
                          OAnew = max(Age), 
                          OAG = TRUE, 
                          extrapLaw = NULL,
                          verbose = TRUE,
                          ...){
  
  # TR: capture args to filter out optional lifetable args to pass on in a named way?
  ExtraArgs <- as.list(match.call())
    
  # get always Mx -----------------------------------------------------------
  
  dates_in <- unique(data$Date)
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
  
  # define kind of input
  types <- c(nMx = NA_real_, nqx = NA_real_, lx = NA_real_)
  data <- data %>% 
              add_column(!!!types[setdiff(names(types), names(data))]) %>%
              mutate(type = ifelse(!is.na(nMx),"m",
                                   ifelse(!is.na(nqx), "q", "l")),
                     Value = ifelse(!is.na(nMx), nMx,
                                   ifelse(!is.na(nqx), nqx, lx)))
  
  # build asked lt with Single arg and get rates. 
  # IW: tried show_query lazing data with dtplyr but not wking yet
  # 
  data <- data %>% 
    arrange(Date, Sex, Age) %>% 
    group_by(Date, Sex) %>%
    summarise(lt_mx_ambiguous(x = Value, 
                              Age = Age, 
                              type = unique(type),
                              Single = Single,
                              OAnew = OAnew, 
                              extrapLaw = extrapLaw, ...)[c("Age","nMx")]) %>% 
    ungroup()

  # need a matrix by sex. Maybe with split in one sentence. sorry the spread,will change w pivot
  nMxm <- data %>% 
              filter(Sex == "m") %>% 
              spread(Date,nMx) %>% 
              select(-Age, -Sex) %>% 
              as.matrix()
  
  nMxf <- data %>% 
              filter(Sex == "f") %>% 
              spread(Date,nMx) %>% 
              select(-Age, -Sex) %>% 
              as.matrix()
  
  # LC at unequal intervals ---------------------------------------------------------

  Age = sort(unique(data$Age))
  ndates_in  <- length(dates_in)
  ndates_out <- length(dates_out)
  nAge       <- length(Age)

  # males
  # lee carter using svd better. That´s what paper suggests
  axm      <- rowSums(log(nMxm))/ndates_in 
  # ktom   <- colSums(log(nMxm))-sum(axm)
  # bxm    <- rowSums(sweep(log(nMxm) - axm, MARGIN = 2, ktom, `*`))/sum(ktom^2)
  nMxm_svd <- svd(log(nMxm)-axm)
  bxm      <- nMxm_svd$u[,1]/sum(nMxm_svd$u[,1])
  ktom     <- nMxm_svd$d[1] * nMxm_svd$v[,1] * sum(nMxm_svd$u[,1])
  cm       <- 0 # parameters
  cm[2]    <- d1m <- (ktom[ndates_in] - ktom[1])/(dates_in[ndates_in] - dates_in[1]) # slope
  cm[1]    <- ktom[1] - cm[2] * dates_in[1] # intercept
  ktm      <- cm[1] + cm[2] * dates_out # interpolate/extrapolate
  k0m      <- cm[1] + cm[2] * dates_in[1] # first year
  
  # females. IW: not happy with this sex duplication lines, maybe a lee-carter function?
  axf      <- rowSums(log(nMxf))/ndates_in
  nMxf_svd <- svd(log(nMxf)-axf)
  bxf      <- nMxf_svd$u[, 1]/sum(nMxf_svd$u[, 1])
  ktof     <- nMxf_svd$d[1] * nMxf_svd$v[, 1] * sum(nMxf_svd$u[, 1])
  cf       <- 0
  cf[2]    <- d1f <- (ktof[ndates_in] - ktof[1])/(dates_in[ndates_in] - dates_in[1])
  cf[1]    <- ktof[1] - cf[2] * dates_in[1]
  ktf      <- cf[1] + cf[2] * dates_out
  k0f      <- cf[1] + cf[2] * dates_in[1]
  
  # ask if prevent divergence and replicate target e0 ---------------------------------------------------------
  
  if (is.null(dates_e0)){ # not rep e0
    
    # basic
    nMxm_hat <- exp(axm + sweep(matrix(bxm, nAge, length(dates_out)),MARGIN=2, ktm,`*`))
    nMxf_hat <- exp(axf + sweep(matrix(bxf, nAge, length(dates_out)),MARGIN=2, ktf,`*`))
    
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
      nMxm_hat_div <- nMxm[,1] * exp(sweep(matrix(bx,nAge,length(dates_out)),MARGIN=2,kt-k0,`*`))
      nMxf_hat_div <- nMxf[,1] * exp(sweep(matrix(bx,nAge,length(dates_out)),MARGIN=2,kt-k0,`*`))
      
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
      ktf_star[j] <- optimize(f = lc_lim_kt_min,
                              interval = c(-20, 20),
                              ax = axf,
                              bx = bxf,
                              age = Age,
                              sex = "f",
                              e0_target = e0f[j])$minimum
    }
    
    # get rates with optim k.
    # TR: this expression could look less R-esoteric. So, we want to multiply k into columns 
    nMxm_hat <- exp(axm + sweep(matrix(bxm, nAge, length(dates_out)), MARGIN = 2, ktm_star, `*`))
    nMxf_hat <- exp(axf + sweep(matrix(bxf, nAge, length(dates_out)), MARGIN = 2, ktf_star, `*`))
  }
  
  # life tables output ------------------------------------------------------------
  
 
  colnames(nMxm_hat) <- dates_out
  colnames(nMxf_hat) <- dates_out
  . = NULL
  # TR: can use ... to pass in optional args. 
  Males_out <-
    lapply(colnames(nMxm_hat), function(x,MX,Age) {
        mx <- MX[, x]
        LT <- lt_mx_ambiguous(x = mx, 
                              Age = Age, 
                              Sex = "m", 
                              axmethod = "un", 
                              a0rule = "ak")
        LT$Sex  <- "m"
        LT$Date <- as.numeric(x)
        LT
      }, MX = nMxm_hat, Age = Age) %>% 
    do.call("rbind", .)
  
  Females_out <-
    lapply(colnames(nMxf_hat), function(x,MX,Age) {
      mx <- MX[, x]
      LT <- lt_mx_ambiguous(x = mx, 
                            Age = Age, 
                            Sex = "f", 
                            axmethod = "un", 
                            a0rule = "ak")
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
                          e0_target){
  Mx_hat <- interp_lc_lim_abk_m(k, ax, bx)
  e0 <- lt_mx_ambiguous(x = Mx_hat, 
                        Age = age, 
                        Sex = sex)$ex[1]
  return(((e0-e0_target)/e0_target)^2)
}

# coale-guo extension nMx from 80 - INCOMPLETE!!!
# TR: can this generalize to arbitrary age classes;
# and also remove position indexing?
# coale_guo_kisker <- function(m, x, xout){
#   #m = nMxf[,1]; x=Age; xout=seq(80,110,5); xextr=seq(80,100,5)
#   
#   if(is_abridged(x)){
#     x_80 = which(x==80)
#     r = .2 * log(m[x_80]/m[x_80-1])
#     m80 = (m[x_80] + m[x_80-1])/2
#     mA = 0.66 + m80 # max at 110
#     s = (log(mA/m80) - 30 * r)/900 
#     m = c(m[1:18], m[18] * exp(r*(xout-80)+s*(xout-80)^2))
#     m = m[c(rep(T,x_80-1), xout %in% xextr)]
#   } else {
#     k_80 = log(m[x_80]/m[x_80-1])
#     obj_fun <- function(s,x,m){
#       k_start = log(m[2]/m[1])
#       m_hat = m[1] * exp(k_start + (x-x[1]) * s)
#       sum(((m-m_hat)/m)^2)
#     }
#     coale_guo_kisker(m,x, xout){
#       m = Mx[80:101,"Male"]
#       mmm = optimize(obj_fun,interval = c(0,1000),x=79:100,m=m)
#     }
#     }
#   return(m)
# }


# get lt for abrev/single ages and m/l/q input
lt_mx_ambiguous <- function(x = NULL, 
                            type = "m", # accepts "q" or "l"
                            Age = NULL, 
                            Sex = NULL, 
                            Single = F,
                            ...){
  if (type == "l"){
    x = lt_id_l_q(x)
    type = "q"
  }
  if (is_abridged(Age)){
    if (type == "m"){
      if(Single == T){
        out <- lt_abridged2single(nMx = x, Age = Age, Sex = Sex, ...)
      }else{
        out <- lt_abridged(nMx = x, Age = Age, Sex = Sex, ...)  
      }
    } else {
      if(Single == T){
        out <- lt_abridged2single(nqx = x, Age = Age, Sex = Sex, ...)
      }else{
        out <- lt_abridged(nqx = x, Age = Age, Sex = Sex, ...)   
      } 
    }
  } else {
    if (type == "m"){
      out <- lt_single_mx(nMx = x, Age = Age, Sex = Sex, ...)
      if(Single != T){
        out <- lt_single2abridged(lx = out$lx, ...)  
      }
    } else {
      out <- lt_single_qx(nqx = x, Age = Age, Sex = Sex, ...)
      if(Single != T){
        out <- lt_single2abridged(lx = out$lx, ...)  
      }
    }
  }
  return(out)
}
