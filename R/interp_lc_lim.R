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
#' @examples
#' 
#'# Using Lee-Carter method for Sweden, assuming only is available mortality rates 
#'for years 1980, 1990, 2000 and 2010
#'# Got from HMD, period rates:
#'Males <- matrix(
#'  c(0.008098, 0.000375, 0.000252, 0.000227, 0.000631, 
#'    0.000996, 0.001241, 0.001317, 0.001759, 0.002725, 0.004439, 0.006827, 
#'    0.010767, 0.017827, 0.028302, 0.046967, 0.075899, 0.123596, 0.191041, 
#'   0.300943, 0.395159, 0.770878, 0.006807, 0.000333, 0.000185, 0.000165, 
#'    0.000624, 0.000859, 0.000955, 0.001081, 0.001375, 0.002098, 0.003206, 
#'    0.004996, 0.009028, 0.014406, 0.023247, 0.039929, 0.067969, 0.110134, 
#'    0.180245, 0.265287, 0.402896, 0.5, 0.004052, 0.000131, 9.8e-05, 
#'    0.000139, 0.000473, 0.000761, 0.000733, 0.000759, 0.000984, 0.001526, 
#'    0.0026, 0.004038, 0.00628, 0.010852, 0.018433, 0.031859, 0.0543, 
#'    0.093116, 0.1638, 0.264152, 0.423773, 0.567665, 0.002712, 0.000159, 
#'    5e-05, 9.9e-05, 0.00034, 0.000661, 0.000726, 0.000674, 0.000703, 
#'    0.001146, 0.001939, 0.003303, 0.005405, 0.009007, 0.014105, 0.023629, 
#'    0.043337, 0.077892, 0.141126, 0.241389, 0.379467, 0.678112),
#'  nrow = length(Age),length(dates_in))
#'Females <- matrix(
#'  c(0.0057, 0.000346, 0.000161, 0.000197, 0.000284, 0.000379, 
#'   0.000485, 0.000651, 0.000949, 0.001666, 0.002325, 0.003553, 0.00563, 
#'    0.00864, 0.013978, 0.024696, 0.045976, 0.084165, 0.143814, 0.232647, 
#'    0.354088, 0.452577, 0.005418, 0.000218, 0.000131, 0.000165, 0.000268, 
#'    0.000325, 0.000433, 0.000536, 0.000776, 0.001246, 0.002037, 0.002875, 
#'    0.005079, 0.007333, 0.012102, 0.020741, 0.038245, 0.07144, 0.128579, 
#'   0.220927, 0.347505, 0.50933, 0.002837, 0.000116, 9.3e-05, 0.000129, 
#'    0.000232, 0.000248, 0.000256, 0.000281, 0.000655, 0.000936, 0.0017, 
#'    0.002723, 0.004497, 0.006826, 0.010464, 0.017542, 0.031884, 0.061459, 
#'    0.115478, 0.204387, 0.335873, 0.484215, 0.002402, 0.00014, 7.7e-05, 
#'    8e-05, 0.000196, 0.000248, 0.000257, 0.000327, 0.000428, 0.000673, 
#'    0.001203, 0.002056, 0.003418, 0.005943, 0.009141, 0.015038, 0.027323, 
#'    0.051479, 0.103552, 0.192314, 0.311358, 0.479488),
#'  nrow = length(Age),length(dates_in))
#'  
#'dates_in <- seq(1980,2010,10)
#'
#'# Estimate rates from mid SXX but using LC with limited data method (4 observed points and unequal intervals)
#'dates_out <- seq(1948,2018,5)
#'
#'# LC with limited data
#'lc_lim_data <- interp_lc_lim(Males, Females, type = "m", 
#'                             dates_in=dates_in, Age=Age, dates_out=dates_out)
#'\dontrun{
#'lc_lim_data %>% ggplot(aes(Age,nMx,col=factor(Year))) + 
#'                geom_step() + scale_y_log10() + facet_wrap(~Sex)
#'}
#'
#'# Avoiding cross-over between sex
#'lc_lim_nondiv <- interp_lc_lim(Males, Females, type = "m", 
#'                               dates_in=dates_in, Age=Age, dates_out=dates_out,
#'                               prev_divergence = T)
#'\dontrun{
#'lc_lim_nondiv %>% ggplot(aes(Age,nMx,col=factor(Year))) + 
#'                  geom_step() + scale_y_log10() + facet_wrap(~Sex)
#'}
#'                        
#'# Using information about e0 in some years, replicate an estimate of e0 in dates_out 
#'# improve performance next step maybe
#'dates_e0 <- seq(1960,2015,5)
#'e0_Males <- readHMDweb("SWE", "E0per", user, pass, fixup = TRUE) %>% 
#'  filter(Year %in% dates_e0) %>% pull(Male)
#'e0_Females <- readHMDweb("SWE", "E0per", user, pass, fixup = TRUE) %>% 
#'  filter(Year %in% dates_e0) %>% pull(Female)
#'lc_lim_fite0 <- interp_lc_lim(Males, Females, type = "m", 
#'                              dates_in=dates_in, Age=Age, dates_out=dates_out,
#'                               dates_e0 = dates_e0,
#'                               e0_Males = e0_Males, 
#'                               e0_Females = e0_Females)
#'\dontrun{                           
#'ggplot() + 
#'  geom_line(data = lc_lim_fite0 %>% filter(Age==0), aes(Year,ex,col=factor(Sex))) + 
#'  geom_point(data = data.frame(Sex = c(rep("m",length(e0_Males)), rep("f",length(e0_Males))),
#'                               ex = c(e0_Males, e0_Females),
#'                              Year = rep(dates_e0,2)),
#'             aes(Year,ex,col=factor(Sex)))
#'}
#'
#'# smooth and/or extend open age group, in this case input is for 80+, and smoothing with Kannisto law
#'# Coale-Guo is coming ;)...
#'lc_lim_extOAg <- interp_lc_lim(Males[1:18,], Females[1:18,],
#'                               type = "m", 
#'                               dates_in=dates_in, Age=Age[1:18], dates_out=dates_out,
#'                               OAG = F, OAnew = 100, extrapLaw = "kannisto")
#'\dontrun{ 
#'lc_lim_extOAg %>% ggplot(aes(Age,nMx,col=factor(Year))) + 
#'  geom_step() + scale_y_log10() + facet_wrap(~Sex)
#'}
#'#End

interp_lc_lim <- function(Males = NULL, 
                          Females = NULL, 
                          type = "m", 
                          dates_in = as.numeric(colnames(Males)), 
                          Age = as.numeric(rownames(Males)), 
                          dates_out = NULL, 
                                     # TR: is there a default? If not, then default can be dates_in, right? 
                                     # IW: It would be the same that user already has, but building lifetables
                          dates_e0 = NULL,
                          e0_Males = NULL, 
                          e0_Females = NULL, 
                          prev_divergence = FALSE, 
                          OAnew = max(Age), 
                          OAG = TRUE, 
                          extrapLaw = NULL,
                          ...){
  
  # TR: capture args to filter out optional lifetable args to pass on in a named way?
  ExtraArgs <- as.list(match.call())
    
  # get always Mx -----------------------------------------------------------
  
  # TR: note, we may as well offer to pass in args to lt_abridged(), which offers lots of control
  # over extrapolation.
  
  # axmethod = "pas", a0rule = "ak", Sex = "m", 
  # IMR = NA, region = "w", mod = TRUE, SRB = 1.05, OAG = TRUE, 
  # OAnew = max(Age), extrapLaw = "kannisto", extrapFrom = max(Age), 
  # extrapFit = Age[Age >= 60 & ifelse(OAG, Age < max(Age), TRUE)]
  # I added a capture of extra args entering the function all at the top.
  # instead of ... pass things in by name, above a list of default values for
  # lt_abridged(), most of which are the same for lt_single*()
  nMxm <- apply(Males, 2, function(x) {
    lt_mx_ambiguous(x, 
                    Age = Age, 
                    type = type, 
                    Sex = "m", 
                    OAnew = OAnew, 
                    extrapLaw = extrapLaw)$nMx}) # IW: add this last 2 args to smooth,
                                                 # but change it as you suggested lines before
  nMxf <- apply(Females, 2, function(x) {
    lt_mx_ambiguous(x, 
                    Age = Age, 
                    type = type, 
                    Sex = "f", 
                    OAnew = OAnew, 
                    extrapLaw=extrapLaw)$nMx}) 
  
  # extend age in case was asked, dont add another if abr/simple
  # TR: Um, this is a bit laborius, right?
  # Age <- lt_mx_ambiguous(x = Females[,1], 
  #                        Age = Age, 
  #                        type = type, 
  #                        Sex = "f", 
  #                        OAnew = OAnew, 
  #                        extrapLaw = extrapLaw)$Age
  if (is_abridged(Age)){
    Age <- inferAgeAbr(nMxf[,1])
  } else {
    Age <- 1:nrow(nMxf)-1
  }
  # LC at unequal intervals ---------------------------------------------------------

  ndates_in  <- length(dates_in)
  ndates_out <- length(dates_out)
  nAge       <- length(Age)
  
  # a basic check
  stopifnot(ncol(Males) ==  ndates_in & ncol(Females) == ndates_in)

  # males
  # lee carter using svd better. That´s what paper suggests
  axm   <- rowSums(log(nMxm))/ndates_in 
  # ktom  <- colSums(log(nMxm))-sum(axm)
  # bxm   <- rowSums(sweep(log(nMxm) - axm, MARGIN = 2, ktom, `*`))/sum(ktom^2)
  nMxm_svd <- svd(log(nMxm)-axm)
  bxm     <- nMxm_svd$u[,1]/sum(nMxm_svd$u[,1])
  ktom    <- nMxm_svd$d[1] * nMxm_svd$v[,1] * sum(nMxm_svd$u[,1])
  cm    <- 0 # parameters
  cm[2] <- d1m <- (ktom[ndates_in] - ktom[1])/(dates_in[ndates_in] - dates_in[1]) # slope
  cm[1] <- ktom[1] - cm[2] * dates_in[1] # intercept
  ktm   <- cm[1] + cm[2] * dates_out # interpolate/extrapolate
  k0m   <- cm[1] + cm[2] * dates_in[1] # first year
  
  # females. IW: not happy with this sex duplication lines, maybe a lee-carter function?
  axf   <- rowSums(log(nMxf))/ndates_in
  nMxf_svd <- svd(log(nMxf)-axf)
  bxf     <- nMxf_svd$u[, 1]/sum(nMxf_svd$u[, 1])
  ktof    <- nMxf_svd$d[1] * nMxf_svd$v[, 1] * sum(nMxf_svd$u[, 1])
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
    
    # stepwise linear intra/extrapolation to target years. 
    # IW: Use interp(). Accepts matrix, so have to rbind and only get 1st row
    # source("R/AGEINT.R")
    e0m = interp(rbind(e0_Males,e0_Males),dates_e0,dates_out,extrap = TRUE)[1,]
    e0f = interp(rbind(e0_Females,e0_Females),dates_e0,dates_out,extrap = TRUE)[1,]

    # avoid divergence: same bx but not kt.
    if(prev_divergence){
      bxm = bxf = (bxm + bxf) * .5
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
    
    # get rates with optim k
    nMxm_hat <- exp(axm + sweep(matrix(bxm,nAge,length(dates_out)),MARGIN=2,ktm_star,`*`))
    nMxf_hat <- exp(axf + sweep(matrix(bxf,nAge,length(dates_out)),MARGIN=2,ktf_star,`*`))
  }
  
  # life tables output ------------------------------------------------------------
  
  # TR: can use ... to pass in optional args.
  colnames(nMxm_hat) <- dates_out
  colnames(nMxf_hat) <- dates_out
  . = NULL
  
  Males_out <-
    lapply(colnames(nMxm_hat), function(x,MX,Age) {
        mx <- MX[, x]
        LT <- lt_mx_ambiguous(x = mx, 
                              Age = Age, 
                              Sex = "m", 
                              axmethod = "un", 
                              a0rule = "ak")
        LT$Sex  <- "m"
        LT$Year <- as.numeric(x)
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
      LT$Year <- as.numeric(x)
      LT
    }, MX = nMxf_hat, Age = Age) %>% 
    do.call("rbind", .)

    out <- rbind(Males_out, Females_out)
  return(out)
  
}

# additional functions -------------------------------------------------------------------

# optimize k for fitting e0
# TR: recommend lc_lim_kt_min(), following name conventions elsewhere in package
lc_lim_kt_min <- function(k,ax,bx,age,sex,e0_target){
  Mx_hat <- exp(ax + bx * k)
  e0 <- lt_mx_ambiguous(x = Mx_hat, Age = age, Sex = sex)$ex[1]
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
lt_mx_ambiguous <- function(x = NULL, type = "m", # accepts "q" or "l"
                  Age = NULL, Sex = NULL, ...){
  if(type=="l"){
    x = lt_id_l_q(x)
    type = "q"
  }
  if(is_abridged(Age)){
    if (type=="m"){
      lt_abridged(nMx = x, Age = Age, Sex = Sex, ...)  
    } else {
      lt_abridged(nqx = x, Age = Age, Sex = Sex, ...)  
    }
  } else {
    if (type == "m"){
      lt_single_mx(nMx = x, Age = Age, Sex = Sex, ...) 
    } else {
      lt_single_qx(nqx = x, Age = Age, Sex = Sex, ...) 
    }
  }
} 

