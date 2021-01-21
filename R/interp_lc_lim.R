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
#' @param Males numeric. Matrix of abridged age classes by year.
#' @param Females numeric. Matrix of abridged age classes by year. 
#' @param type character. Either rates `"m"`, conditionated probabilities `"q"` or survival function `"l"` included in matrices
#' @param dates_in numeric. Vector of observed decimal years for input matrix. Same number of columns as input matrices.
#' @param Age numeric. Vector with inferior limit of abridged age classes with data. Same number of rows than input matrix.
#' @param Sex character. Either male `"m"`, female` "f"`.
#' @param dates_out numeric. Vector of decimal years to interpolate or extrapolate.
#' @param ... Arguments passed to `\link{lt_abridged}`.
#' @param dates_e0 numeric. Vector of decimal years where `"e_0"` should be fitted when apply method.
#' @param e0_Males numeric. Vector of life expectancy by year to be fitted. Same length than `"dates_e0"`.
#' @param e0_Females numeric. Vector of life expectancy by year to be fitted. Same length than `"dates_e0"`.
#' @param prev_divergence logical. Whether or not prevent divergence and sex crossover. Default `FALSE.`
#' @param OAG logical. Whether or not the last element of `nMx` (or `nqx` or `lx`) is an open age group. Default `TRUE.`
#' @param extrapLaw character. If extrapolating, which parametric mortality law should be invoked? Options include
#'   `"Coale-Guo"` or other functions in `"\link{lt_abridged}"` function.
#' @inheritParams lt_a_closeout #??? IW: what put here
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
#' @examples
#' 


interp_lc_lim <- function(Males = NULL, 
                          Females = NULL, 
                          type = "m", 
                          dates_in = as.numeric(colnames(Males)), 
                          Age = as.numeric(rownames(Males)), 
                          dates_out = NULL, 
                                     # TR: is there a default? If not, then default can be dates_in, right? 
                                     # IW: It would be the same that user already has, but building lifetables
                          # In which case it's a cheap smoother. IW: I don´t understand this
                          dates_e0 = NULL,
                          e0_Males = NULL, 
                          e0_Females = NULL, 
                          prev_divergence = FALSE, 
                          ...){
  
  # a basic check
  stopifnot(ncol(Males) ==  length(dates_in) & ncol(Females) == length(dates_in))
  
  # same fun for abrev or simple --------------------------------------------
  
  # probably a better name here
  lt_mx <- function(x = NULL, type = "m", # accepts "q" or "l"
                    Age = NULL, Sex = NULL, ...){
    if(type=="l"){
      x = lt_id_l_d(x)/x[1]
      type = "q"
    }
    if(is_abridged(Age)){
      if(type=="m"){
        lt_abridged(nMx = x, Age = Age, Sex = Sex, ...)  
      }else{
        lt_abridged(nqx = x, Age = Age, Sex = Sex, ...)  
      }
    }else{
      if(type=="m"){
        lt_single_mx(nMx = x, Age = Age, Sex = Sex, ...) 
      }else{
        lt_single_qx(nqx = x, Age = Age, Sex = Sex, ...) 
      }
    }
  } 
  
  # get always Mx and smooth it -----------------------------------------------------------
  
  # TR: note, we may as well offer to pass in args to lt_abridged(), which offers lots of control
  # over extrapolation.
  # IW: added ...
  
  # TR: Perhaps coale-guo should just become an extra option for lt_rule_m_extrapolate()?
  # IW: decided to be included in ... (?)
  
  nMxm <- apply(Males, 2, function(x) {
    lt_mx(x, Age = Age, type = type, Sex = "m", ...)$nMx})
  nMxf <- apply(Females, 2, function(x) {
    lt_mx(x, Age = Age, type = type, Sex = "f", ...)$nMx}) 

  # LC at unequal intervals ---------------------------------------------------------
  
  ndates_in = length(dates_in)
  ndates_out = length(dates_out)
  nAge = length(Age)
  
  # males
  # lee carter using svd better maybe? That´s what paper suggests
  # mx_svd <- svd(log(nMxm)-axm)
  # bx     <- mx_svd$u[, 1]/sum(mx_svd$u[, 1])
  # kt     <- mx_svd$d[1] * mx_svd$v[, 1] * sum(mx_svd$u[, 1])
  axm   <- rowSums(log(nMxm))/ndates_in 
  ktom  <- colSums(log(nMxm))-sum(axm)
  bxm   <- rowSums(sweep(log(nMxm) - axm, MARGIN = 2, ktom, `*`))/sum(ktom^2)
  cm    <- 0
  cm[2] <- d1m <- (ktom[ndates_in] - ktom[1])/(dates_in[ndates_in] - dates_in[1]) # slope
  cm[1] <- ktom[1] - cm[2] * dates_in[1] # intercept
  ktm   <- cm[1] + cm[2] * dates_out # interpolate/extrapolate
  k0m   <- cm[1] + cm[2] * dates_in[1] # first year
  
  # females. Not happy with this sex duplication lines
  axf   <- rowSums(log(nMxf))/ndates_in
  ktof  <- colSums(log(nMxf))-sum(axf)
  bxf   <- rowSums(sweep(log(nMxf) - axf, MARGIN = 2, ktof, `*`))/sum(ktof^2)
  cf    <- 0
  cf[2] <- d1f <- (ktof[ndates_in] - ktof[1])/(dates_in[ndates_in] - dates_in[1])
  cf[1] <- ktof[1] - cf[2] * dates_in[1]
  ktf   <- cf[1] + cf[2] * dates_out
  k0f   <- cf[1] + cf[2] * dates_in[1]
  
  # ask if prevent divergence and replicate target e0 ---------------------------------------------------------
  
  if(is.null(dates_e0)){ # not rep e0
    
    # basic
    nMxm_hat <- exp(axm + sweep(matrix(bxm, nAge, length(dates_out)),MARGIN=2, ktm,`*`))
    nMxf_hat <- exp(axf + sweep(matrix(bxf, nAge, length(dates_out)),MARGIN=2, ktf,`*`))
    
    # avoid divergence extrapolating to 1950: same bx and kt
    if (prev_divergence){
      kt = (ktm + ktf) * .5 # equal size male and female
      # possible error in line 335, should be bxm = (bxm + bxf) * .5
      for (j in 1:ndates_out){
        for (i in 1:nAge){
          bxm[i] = (bxm[i] + bxf[i]) * .5
        }
      }
      k0 = (k0m + k0f) * .5 # same start point
      nMxm_hat_div <- nMxm[,1] * exp(sweep(matrix(bxm,nAge,length(dates_out)),MARGIN=2,kt-k0,`*`))
      # TR: maybe use bxf here? IW: already both sex unified. vba code rares. For sure improve this
      nMxf_hat_div <- nMxf[,1] * exp(sweep(matrix(bxm,nAge,length(dates_out)),MARGIN=2,kt-k0,`*`))
      # only for those years extra to 1950
      Years_extrap <- dates_out<min(dates_in)
      nMxm_hat[,Years_extrap] <- nMxm_hat_div[,Years_extrap]
      nMxf_hat[,Years_extrap] <- nMxf_hat_div[,Years_extrap]
      
    }
  } else { # fit e0 at each target year
    
    # stepwise linear intra/extrapolation to target years
    # TR: can you just use DemoTools::interp() here?
    # it's build for Lexis surfaces. There's also stats::approx()
    # IW: interp modified. Accepts matrix, so have to rbind and only get 1st row
    e0m = interp(rbind(e0_Males,e0_Males),dates_e0,dates_out,extrap = T)[1,]
    e0f = interp(rbind(e0_Females,dates_e0),dates_out,extrap = T)[1,]
    
    # avoid divergence extrapolating to 1950 (?): same bx but not kt
    if(prev_divergence){
      bxm = bxf = (bxm + bxf) * .5
    }
    
    # males. Optimize kt for each LC extrap/interp
    ktm_star = c()
    for (j in 1:ndates_out){
      d1m = d1f # error here in the LC female output
      if (d1m >= 0){
        k1 = k0m
      } else {
        k1 = k0m + d1m * (j+1)
      }
      k0 = k0m + d1m * j
      k0 = ifelse(abs(k1) > .1, abs(k1), .1)
      ktm_star[j] <- optimize(f = lc_lim_kt_min,
                              interval = c(-20, 20),
                              ax = axm,
                              bx = bxm,
                              age = Age,
                              sex = "m",
                              e0_target=e0m[j])$minimum
    }
    
    # females
    ktf_star = c()
    for(j in 1:ndates_out){
      #j=1
      d1f = 0 # error here in the LC female output
      if(d1f>=0){
        k1 = k0f
      }else{
        k1 = k0f + d1f * (j + 1)
      }
      k0 = k0f + d1f * j
      k0 = ifelse(abs(k1) > .1, abs(k1), .1)
      ktf_star[j] <- optimize(f = lc_lim_kt_min,
                              interval = c(-20, 20),
                              ax = axf,
                              bx = bxf,
                              age = Age,
                              sex = "f",
                              e0_target = e0f[j])$minimum
    }
    
    # get rates with optim k
    nMxm_hat <- exp(axm + sweep(matrix(bxm,nAge,length(dates_out)),MARGIN=2,ktm_star,`*`))
    nMxf_hat <- exp(axf + sweep(matrix(bxf,nAge,length(dates_out)),MARGIN=2,ktf_star,`*`))
  }
  
  # life tables output ------------------------------------------------------------
  # TR: can use ... to pass in optional args.
  out <- rbind(
    do.call(rbind, apply(nMxm_hat, 2, function(x) {
      lt_mx(x, Age = Age, Sex = "m", ...)})),
    do.call(rbind, apply(nMxf_hat, 2, function(x) {
      lt_mx(x, Age = Age, Sex = "f", ...)}))) %>% 
    mutate(Year = rep(sort(rep(dates_out,nAge)), 2), 
           Sex = c(rep("m", nAge * ndates_out),
                   rep("f", nAge * ndates_out))) %>% 
    select(Year, Sex, everything())
  return(out)
  
}

# additional functions -------------------------------------------------------------------

# optimize k for fitting e0

# TR: recommend lc_lim_kt_min(), following name conventions elsewhere in package
lc_lim_kt_min <- function(k,ax,bx,age,sex,e0_target){
  Mx_hat <- exp(ax + bx * k)
  e0 <- lt_mx(nMx = Mx_hat, Age = age, Sex = sex)$ex[1]
  return(((e0-e0_target)/e0_target)^2)
}

# coale-guo extension nMx from 80 (improve it)
# TR: can this generalize to arbitrary age classes;
# and also remove position indexing?
coale_guo <- function(m,age=Age){
  r = .2 * log(m[18]/m[17])
  m80 = (m[18]+m[17])/2
  mA = 0.66 + m80
  s = (log(mA/m80) - 30 * r)/900
  h = (19:22-18) * 5
  m = c(m[1:18], m[18] * exp(r*h+s*h^2))
  return(m)
}




