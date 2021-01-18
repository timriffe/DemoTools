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
#' @param dates_e0 numeric. Vector of decimal years where `"e_0"` should be fitted when apply method.
#' @param e0_Males numeric. Vector of life expectancy by year to be fitted. Same length than `"dates_e0"`.
#' @param e0_Females numeric. Vector of life expectancy by year to be fitted. Same length than `"dates_e0"`.
#' @param prev_divergence logical. Whether or not prevent divergence and sex crossover. Default `FALSE.`
#' @param OAG logical. Whether or not the last element of `nMx` (or `nqx` or `lx`) is an open age group. Default `TRUE.`
#' @param extrapLaw character. If extrapolating, which parametric mortality law should be invoked? Options include
#'   `"Coale-Guo"` or other functions in `"lt_abridged"` function.
#' @inheritParams lt_a_closeout #???
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
#' Date
#' Data taken from UN spreadsheet
#' # Years with observed rates
#' Years <- c(1980.671, 1991.668, 1996.586, 2000.586, 2010.71)
#' # Age groups (inferior limit of classes)
#' Age <- c(0,1,seq(5,100,5))
#' # Years to intra/extrapolate, following midpoints for WPP19
#' dates_out <- 1948 + 1:14*5
#' Age (row) by Year (col) arrays, by sex
#' Males <- matrix(c(0.058872, 0.002992, 0.000697, 0.000703, 0.001424, 
#'                   0.002102, 0.002519, 0.003109, 0.004072, 0.005968, 0.00739, 0.010927, 
#'                   0.013366, 0.018798, 0.028653, 0.037304, 0.052714, 0.059629, 0.070922, 
#'                   0.093256, 0.135563, 0.217859, 0.030825, 0.001351, 0.000517, 0.000528, 
#'                   0.001434, 0.002436, 0.002954, 0.003341, 0.003971, 0.004966, 0.006267, 
#'                   0.009662, 0.012983, 0.019135, 0.024503, 0.032664, 0.047827, 0.057952, 
#'                   0.073104, 0.099948, 0.148105, 0.237862, 0.026198, 0.001109, 0.000457, 
#'                   0.000555, 0.001571, 0.002404, 0.003012, 0.003674, 0.004129, 0.005016, 
#'                   0.006223, 0.008328, 0.012217, 0.017762, 0.025149, 0.032561, 0.042365, 
#'                   0.057642, 0.080202, 0.116701, 0.177582, 0.282593, 0.018484, 0.000883, 
#'                   0.000382, 0.000455, 0.001646, 0.002304, 0.002467, 0.003097, 0.003724, 
#'                   0.004507, 0.005908, 0.00794, 0.010738, 0.016865, 0.022493, 0.032624, 
#'                   0.040211, 0.051478, 0.068234, 0.09696, 0.147703, 0.241212, 0.01295, 
#'                   0.00063, 0.000332, 0.000433, 0.001641, 0.002581, 0.002578, 0.002547, 
#'                   0.00289, 0.004012, 0.005381, 0.007316, 0.009889, 0.013273, 0.018334, 
#'                   0.028212, 0.03749, 0.052073, 0.073922, 0.109615, 0.169785, 0.274699),22,5)
#' Females <- matrix(c(0.045269, 0.002704, 0.000507, 0.00046, 0.000734, 
#'                     0.000895, 0.001126, 0.001495, 0.002197, 0.003143, 0.003983, 0.005939, 
#'                     0.007469, 0.01166, 0.018486, 0.026548, 0.042649, 0.050858, 0.063509, 
#'                     0.086965, 0.130587, 0.215029, 0.023838, 0.001154, 0.000358, 0.000318, 
#'                     0.000502, 0.000698, 0.000918, 0.001144, 0.001572, 0.002207, 0.003151, 
#'                     0.005038, 0.007183, 0.011023, 0.014718, 0.022267, 0.035953, 0.048153, 
#'                     0.066424, 0.097196, 0.150869, 0.248412, 0.020248, 0.000933, 0.00031, 
#'                     0.000339, 0.000525, 0.000652, 0.000901, 0.001251, 0.001599, 0.00223, 
#'                     0.00313, 0.004514, 0.007125, 0.01058, 0.015764, 0.021294, 0.032344, 
#'                     0.049166, 0.07543, 0.117877, 0.18764, 0.304247, 0.014603, 0.000768, 
#'                     0.000271, 0.000287, 0.000487, 0.000565, 0.000715, 0.001059, 0.001481, 
#'                     0.002049, 0.002936, 0.004201, 0.006039, 0.009984, 0.013853, 0.021179, 
#'                     0.02809, 0.042159, 0.064247, 0.100939, 0.163497, 0.273028, 0.010488, 
#'                     0.000521, 0.00025, 0.00029, 0.000453, 0.000581, 0.000725, 0.000901, 
#'                     0.001171, 0.001816, 0.002734, 0.003782, 0.005293, 0.007575, 0.011174, 
#'                     0.018559, 0.026524, 0.041711, 0.066135, 0.106604, 0.174691, 0.291021),22,5)

#' # Some visual review of input data - some inconsistencies by age and year
#' plot(Age, Males[,5],t="o",log="y",main="nMx input Males")
#' lines(Age,Males[,1],t="o",log="y",col=2)
#' plot(Age, Males[,5],t="o",log="y",main="nMx input Females")
#' lines(Age,Males[,1],t="o",log="y",col=2)
#' 
#' # LC with limited data
#' lc_lim_data <- interp_lc_lim(Males, Females, "m", dates_in = Years, Age, dates_out)
#' 
#' # LC with  limited data and non-divergence btw sex
#' lc_lim_nondiv <- interp_lc_lim(Males, Females, "m", dates_in = Years, Age, dates_out,
#'                         prev_divergence = T)
#'                         
#' # LC with limited data. rReplicate e0  (see performance issue)
#' lc_lim_fite0 <- interp_lc_lim(Males, Females, "m", Years, Age, dates_out,
#'                        dates_e0 = c(1988,1995,2001,2012),
#'                        e0_Males = 77:80, 
#'                        e0_Females = 77:80)
#'                        
#' # LC with limited data. Smoothing older ages rates
#' lc_lim_CoaleGuo <- interp_lc_lim(Males,Females, type = "m", Years, Age, dates_out,
#'                           extrapLaw = "Coale-Guo")

interp_lc_lim <- function(Males = NULL, 
                     Females = NULL, 
                     type = "m", 
                     dates_in = as.numeric(colnames(Males)), 
                     Age = as.numeric(rownames(Males)), 
                     dates_out, # TR: is there a default? If not, then default can be dates_in, right? 
                                # In which case it's a cheap smoother
                     dates_e0 = NULL,
                     e0_Males = NULL, 
                     e0_Females = NULL, 
                     prev_divergence = FALSE, 
                     OAnew = max(Age), 
                     OAG = TRUE, 
                     extrapLaw = NULL){

  # get always Mx -----------------------------------------------------------
  
  # TR: note, we may as well offer to pass in args to lt_abridged(), which offers lots of control
  # over extrapolation.
  
  # TR Perhaps coale-guo should just become an extra option for lt_rule_m_extrapolate()?
  
  
  if (type == "m"){
    nMxm = Males
    nMxf = Females
    } else {
      if (type == "l"){
        # TR: also lt_abridged() accepts lx as an input arg. So this could loo like the below
        Males = apply(Males, 2, function(x){lt_id_l_d(x)/x})
        Females = apply(Females, 2, function(x){lt_id_l_d(x)/x})
      } # is q
      nMxm <- apply(Males, 2, function(x) {
                    lt_abridged(nqx = x, Age = Age, Sex = "m", axmethod = "un")$nMx})
      nMxf <- apply(Females, 2, function(x) {
                    lt_abridged(nqx = x, Age = Age, Sex = "f", axmethod = "un")$nMx})
    }
  
  # smooth to 100+ from 80-------------------------------------------------------
  if(!is.null(extrapLaw)){
    if(extrapLaw == "Coale-Guo"){ # always from 80
      nMxm = apply(Males,2, coale_guo, age=Age)
      nMxf = apply(Females,2, coale_guo, age=Age)
    }else{
      # call lt_rule_m_extrapolate
    }
    Age = c(0,1,seq(5,100,5))
  }
  
  # LC at unequal intervals ---------------------------------------------------------
  
  nYears = length(dates_in)
  nYears_Target = length(dates_out)
  nAge = length(Age)
  
  # males
  # lee carter using svd better maybe? That´s what paper suggests
    # mx_svd <- svd(log(nMxm)-axm)
    # bx     <- mx_svd$u[, 1]/sum(mx_svd$u[, 1])
    # kt     <- mx_svd$d[1] * mx_svd$v[, 1] * sum(mx_svd$u[, 1])
  axm   <- rowSums(log(nMxm))/nYears
  ktom  <- colSums(log(nMxm))-sum(axm)
  bxm   <- rowSums(sweep(log(nMxm) - axm, MARGIN = 2, ktom, `*`))/sum(ktom^2)
  cm    <- 0
  cm[2] <- d1m <- (ktom[nYears] - ktom[1])/(dates_in[nYears] - dates_in[1]) # slope
  cm[1] <- ktom[1] - cm[2] * dates_in[1] # intercept
  ktm   <- cm[1] + cm[2] * dates_out # interpolate/extrapolate
  k0m   <- cm[1] + cm[2] * dates_in[1] # first year
  
  # females. Not happy with this sex duplication lines
  axf   <- rowSums(log(nMxf))/nYears
  ktof  <- colSums(log(nMxf))-sum(axf)
  bxf   <- rowSums(sweep(log(nMxf) - axf, MARGIN = 2, ktof, `*`))/sum(ktof^2)
  cf    <- 0
  cf[2] <- d1f <- (ktof[nYears] - ktof[1])/(dates_in[nYears] - dates_in[1])
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
      for (j in 1:nYears_Target){
        for (i in 1:nAge){
          bxm[i] = (bxm[i] + bxf[i]) * .5
        }
      }
      k0 = (k0m + k0f) * .5 # same start point
      nMxm_hat_div <- nMxm[,1] * exp(sweep(matrix(bxm,nAge,length(dates_out)),MARGIN=2,kt-k0,`*`))
      # TR: maybe use bxf here?
      nMxf_hat_div <- nMxf[,1] * exp(sweep(matrix(bxm,nAge,length(dates_out)),MARGIN=2,kt-k0,`*`))
      # only for those years extra to 1950
      Years_extrap <- dates_out<min(Years)
      nMxm_hat[,Years_extrap] <- nMxm_hat_div[,Years_extrap]
      nMxf_hat[,Years_extrap] <- nMxf_hat_div[,Years_extrap]
  
    }
  } else { # fit e0 at each target year
    
    # stepwise linear intra/extrapolation to target years
    # TR: can you just use DemoTools::interp() here?
    # it's build for Lexis surfaces. There's also stats::approx()
    e0m = interp_linear(dates_e0,e0_Males,dates_out)
    e0f = interp_linear(dates_e0,e0_Females,dates_out)
    
    # avoid divergence extrapolating to 1950 (?): same bx but not kt
    if(prev_divergence){
      bxm = bxf = (bxm + bxf) * .5
    }
    
    # males. Optimize kt for each LC extrap/interp
    ktm_star = c()
    for (j in 1:nYears_Target){
      d1m = d1f # error here in the LC female output
      if (d1m >= 0){
        k1 = k0m
      } else {
        k1 = k0m + d1m * (j+1)
      }
      k0 = k0m + d1m * j
      k0 = ifelse(abs(k1) > .1, abs(k1), .1)
      ktm_star[j] <- optimize(f = kt_obj_fun,
                              interval = c(-20, 20),
                              ax = axm,
                              bx = bxm,
                              age = Age,
                              sex = "m",
                              e0_target=e0m[j])$minimum
    }
    
    # females
    ktf_star = c()
    for(j in 1:nYears_Target){
      #j=1
      d1f = 0 # error here in the LC female output
      if(d1f>=0){
        k1 = k0f
      }else{
        k1 = k0f + d1f * (j + 1)
      }
      k0 = k0f + d1f * j
      k0 = ifelse(abs(k1) > .1, abs(k1), .1)
      ktf_star[j] <- optimize(f = kt_obj_fun,
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
      lt_abridged(nMx = x, Age = Age, Sex = "m", axmethod = "un",a0rule = "ak")})),
    do.call(rbind, apply(nMxf_hat, 2, function(x) {
      lt_abridged(nMx = x, Age = Age, Sex = "f", axmethod = "un",a0rule = "ak")}))) %>% 
    mutate(Year = rep(sort(rep(dates_out,nAge)), 2), 
           Sex = c(rep("m", nAge * nYears_Target),
                   rep("f", nAge * nYears_Target))) %>% 
    select(Year, Sex, everything())
  return(out)

}

# additional functions -------------------------------------------------------------------

# interpolation between pairwise points
# TR: we have interp(), for example. Does this do the same?
interp_linear <- function(x,y,x_star){
  out = c()
  for(i in 1:length(x_star)){
    xt=x_star[i]
    # in case is out or in observed range 
    if(between(xt,min(x),max(x))){
      x1 = max(x[x<=xt])
      x2 = min(x[x>=xt])
    }else{
      nodes_interp = which(abs(xt-x) %in% sort(abs(xt-x))[1:2])
      x1 = x[nodes_interp[1]]
      x2 = x[nodes_interp[2]]
    }
    # avoid dividing by zero if target years matches observed ones 
    out[i] = y[x==x1] + 
            ifelse(x[x==x2] - x[x==x1]==0,0,(y[x==x2] - y[x==x1])/(x[x==x2] - x[x==x1])) * 
            (xt-x[x==x1])
  }
  return(out)
}

# optimize k for fitting e0

# TR: recommend lc_lim_kt_min(), following name conventions elsewhere in package
kt_obj_fun <- function(k,ax,bx,age,sex,e0_target){
  Mx_hat <- exp(ax + bx * k)
  e0 <- lt_abridged(nMx = Mx_hat, Age = age, Sex = sex, axmethod = "un", a0rule = "cd")$ex[1]
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

