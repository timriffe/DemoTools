#' estimate adult mortality using two censuses.
#' @description 
#' @param c1 numeric vector. Population counts in five-age groups, from first census with an exact reference date.
#' @param c2 numeric vector. Population counts in five-age groups, from second census with an exact reference date.
#' @param date1 Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}. Reference date for first census.
#' @param date2 Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}. Reference date for second census.
#' @param age integer vector. Lower bound of age groups from first census. Last age is assumed OAG.
#' @param age integer vector. Lower bound of age groups from second census. Last age is assumed OAG.
#' @param sex character string. Either `"m"` for males, `"f"` for females (default). 
#' @param mlt_family character string. Model life table family: `"Chilean"`, `"East"`, `"Far_East_Asian"`, `"General"`, `"Latin"`, `"North"`, `"South"`, `"South_Asian"` or `"West"` (default).
#' @param method character. Options include:
#' * `"match"` use survival ratios for finding best mlt level (default).
#' * `"bproj"` back-projection from second census (UN, 1983).
#' * `"fproj"` forward-projection from first census (UN, 1983).
#' * `"var-r"` variable-r method (Preston & Benet, 1986).
#' * `"feeney"` smoothing survival ratios into survival function, implementation by Feeney (UN, 2002).
#' * `logit` logit model changing level (Preston 2001, section 11.5.2)
#' @param ages_fit integer vector. Ages to be considered when calculating the median of implied level by age. By default 10 to 70, but depends data and method.  
#' @param e0_accept integer vector. Range acceptable when calculating the median of implied level by age (avoid non-possible extrapolations). By deafult between 20 and 100.
#' @param q01_q05 numeric. Both values in case final life table be computed considering child mortality input values.
#' @param span_pre_smooth numeric. smooth is applied to survival ratios by age in case is not null. Value is `span` parameter from `loess` function (between 0 and 1). Default `FALSE`.
#' @param mlt_e0_logit_feeney numeric. Level as reference in `logit` and `feeney` methods.
#' @param mlt_input_data data.frame. By default is used data from  \code{Morcast} package. But specific model can be included (be sure same columns are included with same names), ant that `mlt_family` be included as category in `type` column.
#' @param extrapLaw character. Which mortality law is chosen for extension to 100 age in lt output. See `Demotools::lt_abridged` for mode details
#' @param HIV_prev numeric. Estimates of population-based estimate of HIV prevalence by age. If some value is assigned, then will be assumed AIDS variation of the method based on Spectrum patterns.
#' @param HIV_art numeric. Estimate of proportion using antiretroviral therapy (ART).
#' @return A list with 45q15, 35q15, the life table as result of the method, the general implicit level, the age-specific implicit level, and nSx from input data.
#' @export
#' @examples
#' \dontrun{
#' # PANAMA 1960-1970 from UN Manual X (Table 174). Get only level.
#' intercensal_survival(
#'     c1 = c(90071,76598,63635,54431,45634,37818,32179,28724,23974,20618,15068,11999,10283,6737,5242,6756),
#'     c2 = c(114017,106944,85253,73381,63010,50924,40885,36115,29409,25360,21775,17632,13004,10061,6690,9873),
#'     date1 = "1960/12/11",
#'     date2 = "1970/05/10",
#'     age1 = seq(0,75,5),
#'     age2 = seq(0,75,5),
#'     sex = "f",
#'     mlt_family = "West",
#'     method = "match",
#'     ages_fit = seq(0,75,5))$selected_level
#' }

# general function ---------------------------------------------------
intercensal_survival <- function(c1,
                                 c2,
                                 date1,
                                 date2,
                                 age1 = seq(0, length(c1) * 5, 5),
                                 age2 = seq(0, length(c2) * 5, 5),
                                 sex = "f",
                                 mlt_family = "CD_West", 
                                 method = "match",
                                 ages_fit = seq(10, 70, 5),
                                 e0_accept = c(20, 100),
                                 q01_q05 = NULL,
                                 mlt_e0_logit_feeney = NULL,
                                 span_pre_smooth = NULL,
                                 mlt_input_data = NULL,
                                 HIV_prev = NULL,
                                 HIV_art = NULL,
                                 extrapLaw = "makeham",
                                 verbose = TRUE){
  
  # check family and method
  mlt_families <- c("CD_East", "CD_North", "CD_South", "CD_West", "UN_Chilean", 
                    "UN_Far_Eastern", "UN_General", "UN_Latin_American", "UN_South_Asian")
  # which families consider
  if(!is.null(mlt_family)){
    mlt_family <- match.arg(mlt_family, mlt_families, several.ok = TRUE)
  }else{
    mlt_family <- mlt_families
  }
  
  # is HIV
  if((is.null(HIV_prev) + is.null(HIV_art)) == 0) is.HIV = TRUE
  if((is.null(HIV_prev) + is.null(HIV_art)) == 2) is.HIV = FALSE
  if((is.null(HIV_prev) + is.null(HIV_art)) == 1) stop("for hiv you need prevalence and art")
  
  # check input method
  method <- match.arg(method, c("match", "bproj", "fproj", "var-r", "feeney", "logit"))
  
  # manage ages: fit to more wider group and minor OAG.
  # Always assume last age is an OAG
  if(length(age1)!=length(c1) | length(age2)!=length(c2)) stop("Not same interval between pop and age")
  age1_int <- max(unique(diff(age1)))
  age2_int <- max(unique(diff(age2)))
  age_int <- max(age1_int, age2_int)  
  OAG <- min(max(age1), max(age2))
  c1 <- data.frame(c1=c1, age1=pmin(trunc(age1/age_int)*age_int, OAG)) %>% dplyr::group_by(age1) %>% dplyr::summarise(c1 = sum(c1)) %>% dplyr::pull(c1)
  c2 <- data.frame(c2=c2, age2=pmin(trunc(age2/age_int)*age_int, OAG)) %>% dplyr::group_by(age2) %>% dplyr::summarise(c2 = sum(c2)) %>% dplyr::pull(c2)
  age <- seq(0, OAG, age_int)
  ages <- length(age)
  
  # fitting ages, important for some methods
  if(is.null(ages_fit)) ages_fit <- age
  
  # find a close date for date1, so difference is a multiple of age_int
  if(!is.numeric(date1)) date1 <- round(DemoTools::dec.date(date1),2)
  if(!is.numeric(date2)) date2 <- round(DemoTools::dec.date(date2),2)
  interc_t <- date2-date1
  if((interc_t)>30) stop("More than 30 years in intercensal period. Too much")
  interc_length <- (interc_t)/age_int
  interc_round <- round(interc_length)
  date1_adj <- date1 + (interc_length-interc_round)*age_int
  
  # translate first census using general growth rate (UN Manual)
  # IW: using age specific growth rate, probably better: c1 <- DemoTools::interpolatePop(c1, c2, date1, date2, DesiredDate = date1_adj, method = "exponential")
  r_hat <- log(sum(c2)/sum(c1))/(date2-date1)
  c1_noadj <- c1
  date1_noadj <- date1
  c1 <- c1 * exp(r_hat * (date1_adj-date1))
  date1 <- date1_adj
  interc_t <- date2-date1
  interc_length <- (interc_t)/age_int
  
  # get survival ratios and cum pop, preparing for methods
  surv_data <- data.frame(age = age, age_int = age_int, c1_noadj = c1_noadj, c1 = c1, c2 = c2) %>%  
    dplyr::mutate(c1_cum = purrr::accumulate(c1, `+`, .dir = "backward"), 
                  c2_cum = purrr::accumulate(c2, `+`, .dir = "backward"),
                  nSx = dplyr::case_when(age < (OAG - interc_t)  ~ dplyr::lead(c2, interc_length)/c1,
                                         T ~ dplyr::lead(c2_cum, interc_length)/c1_cum),
                  n = interc_t,
                  age_int = ifelse(age == OAG, NA, age_int))
  
  # give a message on how many survival ratios greater than 1 you have, no matter ages_fit
  ages_nSx_greater1 <- surv_data$age[surv_data$nSx>1 & !is.na(surv_data$nSx)]
  if(length(ages_nSx_greater1)>0 & verbose) message(
    paste0("survival ratios greater than 1 for group ages starting on ", paste(ages_nSx_greater1,collapse=", "))
  )
  
  # pre-smooth on nSx
  if(!is.null(span_pre_smooth)){
    surv_data$nSx[!is.na(surv_data$nSx)] <- loess(nSx ~ age, data = surv_data, span = span_pre_smooth) %>% 
      predict(age = surv_data$age[!is.na(surv_data$nSx)])
  }
  
  # select family and standardize for same age_int, OAG and risk interval
  if(!is.null(mlt_input_data)){
    if(!all(colnames(mlt_input_data) %in% colnames(MortCast::MLTlookup))) stop("Be sure to have same col names and cathegories than MortCast::MLTlookup")
    mlt_data <- mlt_data_input %>% mutate(type = "user", e0 = "user")
    mlt_e0_logit_feeney <- "user"
  }else{
    if(!is.HIV){
      this_sex <- ifelse(sex == "f", 2, 1)
      mlt_this_family_all_ages <- MortCast::MLTlookup %>% filter(type %in% mlt_family, sex == this_sex)
      # e0 should be rounded to proximate available level
      mlt_e0_logit_feeney <- round(mlt_e0_logit_feeney/2.5,0)*2.5
    }else{
      if(is.null(q01_q05[2]) | is.null(HIV_art)) stop("needs 5q0 and/or HIV ART")
      # Spectrum model metadata
      load("R/modsr-vr-dhs-spectrum-25-04.RData")
      this_sex <- ifelse(sex == "f", "female", "male")
      mlt_this_family_all_ages <- lapply(seq(.1,.9,.05), function(x){
        hiv_svd_comp_x <- predictNQX(this_sex, 
                                     cm = q01_q05[2], 
                                     am = x, 
                                     hiv = HIV_prev, 
                                     art = HIV_art, 
                                     adult = "q45") %>% pull()
        lx_hiv_svd_comp_x <- lt_id_q_l(expit(hiv_svd_comp_x))
        lt_abridged(lx = lx_hiv_svd_comp_x[c(0,1,seq(5,100,5))+1], Age = c(0,1,seq(5,100,5))) %>% 
          mutate(type = "HIVSpectrum", sex = ifelse(sex== "f", 2, 1)) %>% 
          group_by(type) %>% 
          mutate(e0 = ex[Age == 0])}) %>% 
        bind_rows() %>% 
        select(type, e0, sex, age = Age, mx = nMx, lx, Lx = nLx)
      mlt_family <- "HIVSpectrum"
      if(method %in% c("feeney", "logit")){
        actual_levels <- unique(mlt_data$e0)
        closer_level <- which(abs(actual_levels-mlt_e0_logit_feeney)==min(abs(actual_levels-mlt_e0_logit_feeney)))
        mlt_e0_logit_feeney <- actual_levels[closer_level]
      }
    }
  }
  
  # standarize mlt to same age structure than data
  this_sex = ifelse(sex == "f", 2, 1)
  this_mlt_family <- mlt_this_family_all_ages %>% 
    dplyr::mutate(age_desired = as.integer(pmin(trunc(age/age_int)*age_int, OAG))) %>% 
    dplyr::group_by(type, e0, age_desired) %>% 
    dplyr::summarise(mxn = sum(mx*Lx)/sum(Lx),
                     lx = max(lx),
                     Lxn = sum(Lx), .groups = 'drop') %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(type, e0) %>% 
    dplyr::mutate(Tx = purrr::accumulate(Lxn, `+`, .dir = "backward"),
                  ex = Tx/lx,
                  nSx = dplyr::case_when(age_desired < (OAG - interc_t)  ~ dplyr::lead(Lxn, interc_length)/Lxn,
                                         T ~ dplyr::lead(Tx, interc_length)/Tx),
                  n = interc_t) %>% 
    dplyr::select(type, e0, age = age_desired, mxn, lx, Lxn, Tx, ex, nSx, n) %>% 
    dplyr::ungroup()
  
  # Find best level for each nSx (Table 171, B.2 in MX)
  if(method == "match" | method == "logit"){
    
    # get closest level for each nSx
    obs_data <- data.frame(age = age, nSx = surv_data$nSx) %>% dplyr::filter(nSx<1)
    mlt_data <- this_mlt_family %>% dplyr::select(type, age, e0, nSx)
    mlt_closest <- interp_level_mlt(obs_data, mlt_data, var_interp = "e0", var_ref = "nSx")
  }
  
  # variable-r 
  if(method == "var-r"){
    
    # get closest level for each nSx
    obs_data <- lt_ManualX_variable_r(age = age, 
                                      Nx1 = pmax(c1_noadj,1), # to avoid 0 pop in some age (age-specific rate goes Inf) 
                                      Nx2 = pmax(c2,1), # to avoid 0 pop in some age (age-specific rate goes Inf) 
                                      ts = date2-date1_noadj, radix = 100000, full_lt = F) %>% dplyr::select(age, ex)
    mlt_data <- this_mlt_family %>% dplyr::select(type, age, e0, ex)
    mlt_closest <- interp_level_mlt(obs_data, mlt_data, var_interp = "e0", var_ref = "ex")
  }
  
  # projection (forw and backw)
  if(method %in% c("bproj", "fproj")){
    
    # do proj and bproj for each level
    surv_data_method <- data.frame(age = age, c1 = c1, c2 = c2) %>% 
      dplyr::mutate(c1_cum = purrr::accumulate(c1, `+`, .dir = "backward"), 
                    c2_cum = purrr::accumulate(c2, `+`, .dir = "backward")) %>% 
      dplyr::inner_join(this_mlt_family %>% select(type, age, e0, nSx_mlt = nSx), by = "age") %>% 
      dplyr::arrange(type, e0, age) %>% 
      dplyr::group_by(type, e0) %>% 
      dplyr::mutate(c1_bproj = lead(c2, interc_length, default = NA)/ nSx_mlt,
                    c2_fproj = case_when(age < OAG ~ lag(c1 * nSx_mlt, interc_length, default = NA),
                                         T ~ lag(c1_cum * nSx_mlt, interc_length, default = NA)),
                    c1_bproj_cum = ifelse(is.na(c1_bproj), NA, purrr::accumulate(c1_bproj, sum, .dir = "backward", na.rm = T)),
                    c2_fproj_cum = ifelse(is.na(c2_fproj), NA, purrr::accumulate(c2_fproj, sum, .dir = "backward", na.rm = T))) %>% ungroup()
    
    # get best level for each direction
    if(method == "bproj"){
      obs_data <- surv_data_method %>% dplyr::select(age, c1_cum) %>% distinct()
      mlt_data <- surv_data_method %>% dplyr::select(age, type, e0, c1_cum = c1_bproj_cum)
      mlt_closest <- interp_level_mlt(obs_data, mlt_data, "e0", "c1_cum")
    }
    if(method == "fproj"){
      obs_data <- surv_data_method %>% dplyr::select(age, c2_cum) %>% distinct()
      mlt_data <- surv_data_method %>% dplyr::select(age, type, e0, c2_cum = c2_fproj_cum)
      mlt_closest <- interp_level_mlt(obs_data, mlt_data, "e0", "c2_cum") %>% 
        dplyr::mutate(e0_interp = lead(e0_interp, interc_length))
    }
  }
  
  # Feeney UN´s manual (2002)
  if(method == "feeney"){
    
    # always apply 5 or 10 interval variants, except greater than 10, where using r-variable variant
    if(!interc_t %in% c(5, 10)) {
      interc_t = date2 - date1_noadj
      surv_data$c1 = surv_data$c1_noadj
      if(verbose) message("r-variable variant was applied for feeney method")
    }
    # get e(OAG-5). If multiple mlt are selected, take the mean
    if(!is.null(mlt_e0_logit_feeney)){
      TOAGless5_l2.5 <-  mlt_this_family_all_ages %>% 
        dplyr::filter(e0 == trunc(mlt_e0_logit_feeney/2.5)*2.5) %>% 
        dplyr::group_by(type) %>% 
        dplyr::mutate(Tx = purrr::accumulate(Lx, `+`, .dir = "backward")) %>% 
        dplyr::summarise(TOAGless5_l2.5 = sum(Tx[age==OAG-5])/sum(lx[age==1]*2.5/4+lx[age==5]*1.5/4)) %>% 
        dplyr::pull()
      mlt_data <- this_mlt_family %>%  
        dplyr::select(type, age, e0, ex)
    }else{
      trunc_age <- OAG - 5
      mlt_data <- this_mlt_family %>% 
        dplyr::filter(age<trunc_age) %>% 
        dplyr::group_by(type, e0) %>% 
        dplyr::mutate(Tx = purrr::accumulate(Lxn, `+`, .dir = "backward"), ex_temp = Tx/lx) %>%  
        dplyr::select(type, age, e0, ex_temp)
      TOAGless5_l2.5 <- NULL
    }
    
    # get ex
    surv_data_method_lx <- intercensal_surv_feeney(c1 = surv_data$c1, 
                                                   c2 = surv_data$c2, 
                                                   age = surv_data$age, 
                                                   age_int = surv_data$age_int,
                                                   interc_t = interc_t,
                                                   TOAGless5_l2.5 = TOAGless5_l2.5)
    names(mlt_data)[ncol(mlt_data)] <- "ex"
    names(surv_data_method_lx)[ncol(surv_data_method_lx)] <- "ex"
    # interpolate
    obs_data <- surv_data_method_lx %>% 
      dplyr::select(age, ex) %>% 
      dplyr::filter(!is.na(ex))
    mlt_closest <- interp_level_mlt(obs_data, mlt_data, "e0", "ex")
  }
  
  # given selected method, get best general level for all ages
  mlt_level <- mlt_closest %>%
    dplyr::filter(age %in% ages_fit) %>%
    mutate(e0_interp = ifelse(e0_interp < e0_accept[1], e0_accept[1],
                              ifelse(e0_interp > e0_accept[2], e0_accept[2],
                                     e0_interp))) %>% 
    dplyr::group_by(type) %>% 
    dplyr::summarise(ages_fitting = paste(range(ages_fit),collapse = "-"),
                     mean_e0 = mean(e0_interp, na.rm = T),
                     median_e0 = median(e0_interp, na.rm = T),
                     coef_var = sd(e0_interp, na.rm = T)/mean_e0,
                     p25 = quantile(e0_interp, probs = .25, na.rm = T),
                     p75 = quantile(e0_interp, probs = .75, na.rm = T))
  
  # get mlt best fit interpolating lx
  lx_mlt_level <- interp_level_mlt(
    mlt_level %>% select(type, e0 = median_e0) %>% expand_grid(age = c(0,1,seq(5,100,5))),
    mlt_this_family_all_ages %>% 
      dplyr::filter(type %in% mlt_family, sex == this_sex) %>% 
      dplyr::select(type, age, e0, lx), 
    "lx", "e0")
  
  # logit method using match results
  if(method == "logit"){
    match_mlt_e0 <- mlt_level %>% dplyr::mutate(mlt_e0 = round(median_e0/2.5)*2.5)
    lx_mlt_level <- map_df(mlt_family, function(X){
      if(is.null(mlt_e0_logit_feeney) | is.HIV){
        # message("logit e0 input is taken from match method level result") 
        mlt_e0 <-  match_mlt_e0$mlt_e0[match_mlt_e0$type==X]
      }else{
        mlt_e0 <- round(mlt_e0_logit_feeney/2.5)*2.5
      }
      this_mlt_family_X <- mlt_this_family_all_ages %>% dplyr::filter(type == X, e0 == mlt_e0) 
      if(is.null(q01_q05)) {
        # message("logit 5p0 input is taken from match method level result") 
        q0_5 <- 1 - this_mlt_family_X$lx[this_mlt_family_X$age==5]/this_mlt_family_X$lx[this_mlt_family_X$age==0]
      }else{
        q0_5 <- q01_q05[2]
      }
      
      # apply method
      e5_hat <- intercensal_surv_Preston_logit_var_r(c1, c2, date1, date2, age, age, 
                                                     sex = sex,
                                                     mlt_family = X,
                                                     mlt_e0 = mlt_e0,
                                                     ages_fit = ages_fit,
                                                     q0_5 = q0_5)$e5
      
      # limit range
      range_e5_family <- range(this_mlt_family$ex[this_mlt_family$age==5])
      if(e5_hat > range_e5_family[2] | e5_hat<range_e5_family[1]){
        e5_hat <- min(max(e5_hat, range_e5_family[2]),range_e5_family[2])
        if(verbose) message("value of e5 out of range, replaced by some of the limits. Review")
      }
      
      # interp considering age 1
      lx_logit_rr <- interp_level_mlt(data.frame(age =  c(0,1,seq(5,100,5)), e5 = e5_hat), 
                                      mlt_this_family_all_ages %>%
                                        dplyr::filter(type == X) %>%
                                        group_by(type, e0) %>% 
                                        dplyr::mutate(Tx = purrr::accumulate(Lx, `+`, .dir = "backward"), ex = Tx/lx) %>%  
                                        dplyr::mutate(e5 = ex[age==5]) %>% 
                                        dplyr::ungroup() %>% 
                                        dplyr::select(type, age, e5, lx), 
                                      "lx", "e5")
      lx_logit_rr$lx_interp <- pmax(lx_logit_rr$lx_interp, 0)
      return(lx_logit_rr)
    })
  }  
  
  # if both infant and child input are ok, impute first ages and recalculate lx
  if(length(q01_q05) == 2){
    lx_mlt_level <- lx_mlt_level %>% 
      dplyr::mutate(lx_interp = case_when(age == 1 ~ 1e5 * (1-q01_q05[1]),
                                          age == 5 ~ 1e5 * (1-q01_q05[2]),
                                          T ~ lx_interp))
  }
  
  # build lt with OAG from last fit age
  lt_selected <- lx_mlt_level %>% 
    dplyr::select(type, age, lx_interp) %>% 
    dplyr::filter(age<=100) %>% 
    split(list(.$type)) %>% 
    purrr::map_df(function(X){
      X <- X[X$lx_interp>0,]
      DemoTools::lt_abridged(lx = X$lx_interp, Age = X$age, Sex = sex, extrapLaw = extrapLaw) %>% 
        dplyr::mutate(type = unique(X$type), .before = 1)
    })
  
  # list of results to return
  mlt_out <- list(Adult = lt_selected %>% 
                    dplyr::group_by(type) %>% 
                    dplyr::summarise(q15_45 = 1-lx[Age==60]/lx[Age==15],
                                     q15_35 = 1-lx[Age==50]/lx[Age==15]),
                  lt_fit = lt_selected,
                  selected_level = mlt_level,
                  rank_match = mlt_closest,
                  surv_data = surv_data)
  if(!is.null(span_pre_smooth)) mlt_out$plotplot_smooth_nSx <- plot(plot_smooth_nSx)
  return(mlt_out)
}


# Preston logit r ---------------------------------------------------------

# based on Preston (2001, section 11.5.2)
intercensal_surv_Preston_logit_var_r <- function(c1,
                                                 c2,
                                                 date1,
                                                 date2,
                                                 age1 = 1:length(c1) - 1,
                                                 age2 = 1:length(c2) - 1,
                                                 sex = "f",
                                                 mlt_family = "CD_West",
                                                 mlt_e0 = 60,
                                                 ages_fit = NULL, 
                                                 q0_5 = NULL,
                                                 verbose = FALSE){
  
  # manage ages: fit to more wider group
  if(is.null(q0_5)) stop("some value in q0_5 is needed.")
  age1_int <- max(unique(diff(age1)))
  age2_int <- max(unique(diff(age2)))
  age_int <- max(age1_int, age2_int)  
  c1 <- data.frame(c1=c1, age1=trunc(age1/age_int)*age_int) %>% group_by(age1) %>% summarise(c1 = sum(c1)) %>% pull(c1)
  c2 <- data.frame(c2=c2, age2=trunc(age2/age_int)*age_int) %>% group_by(age2) %>% summarise(c2 = sum(c2)) %>% pull(c2)
  age <- seq(0, min(max(age1), max(age2)), age_int)
  OAG <- max(age)
  ages <- length(age)
  if(!is.numeric(date1)) date1 <- round(DemoTools::dec.date(date1),2)
  if(!is.numeric(date2)) date2 <- round(DemoTools::dec.date(date2),2)
  # if((date2-date1)>15) message("more than 15 years in intercensal period")
  if(is.null(ages_fit)) ages_fit <- age 
  if(!is.numeric(date1)) date1 <- round(DemoTools::dec.date(date1),2)
  if(!is.numeric(date2)) date2 <- round(DemoTools::dec.date(date2),2)
  interc_t <- date2-date1
  
  # check family and method
  mlt_families <- c("CD_East", "CD_North", "CD_South", "CD_West", "UN_Chilean", 
                    "UN_Far_Eastern", "UN_General", "UN_Latin_American", "UN_South_Asian")
  mlt_family <- match.arg(mlt_family, mlt_families)
  
  # mlt data
  this_sex <- ifelse(sex == "f", 2, 1)
  mlt_data <- MortCast::MLTlookup %>% 
    dplyr::filter(type == mlt_family, e0 == mlt_e0, sex == this_sex) %>% 
    dplyr::mutate(lx = lx,
                  age_desired = as.integer(pmin(trunc(age/age_int)*age_int, OAG))) %>% 
    dplyr::group_by(type, e0, age_desired) %>% 
    dplyr::summarise(mxn = sum(mx*Lx)/sum(Lx),
                     lx = max(lx),
                     Lxn = sum(Lx), .groups = 'drop') %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(type, e0) %>% 
    dplyr::mutate(Tx = purrr::accumulate(Lxn, `+`, .dir = "backward"),
                  ex = Tx/lx) %>% 
    dplyr::select(type, e0, age = age_desired, mxn, lx, Lxn, Tx, ex) %>% 
    dplyr::filter(age %in% age[age!=OAG]) %>% 
    dplyr::mutate(p5a = lx/lx[age==5], q5a = 1-p5a, `q5a/p5a`=q5a/p5a)
  
  # surv data
  surv_data <- data.frame(age = age, age_int = age_int, c1 = c1, c2 = c2) %>%
    dplyr::filter(age!=OAG) %>% 
    dplyr::mutate(r = log(c2/c1)/(date2-date1),
                  r_cum = purrr::accumulate(r, `+`, .dir = "forward"),
                  r_cum = dplyr::lag(r_cum, default = 0),
                  c1.5 = (c2-c1)/(r*interc_t),
                  c_a = (dplyr::lag(c1.5)+c1.5)/(interc_t * sum(c1.5)),
                  r_cum_5 = r_cum * 5,
                  ratio_rc = exp(-r_cum_5)*(1-q0_5)/c_a)
  
  # regression
  regression_data <- inner_join(surv_data, mlt_data, by = "age") %>% dplyr::filter(age >= 5)
  model <- lm(ratio_rc ~ `q5a/p5a`, data = regression_data, subset = age %in% ages_fit)
  plot_obs_fit <- ggplot2::ggplot(regression_data %>% dplyr::filter(age %in% ages_fit), aes(ratio_rc, `q5a/p5a`)) + 
    ggplot2::geom_point() + 
    ggplot2::stat_smooth(method = "lm")
  
  # new level
  coeff <- model$coefficients
  K <- coeff[2]/coeff[1]
  
  # don´t accept and alpha<.1 and alpha>3
  if(K > exp(3) | K<exp(.1)){
    K <- min(max(exp(.1), K), exp(3))
    if(verbose) message("Revise age fit selection. results were truncated for alpha on limits .1 or 3")
  }
  regression_data <- regression_data %>% 
    dplyr::mutate(p5a_hat = 1/(1+`q5a/p5a`*K), nLx = (p5a_hat + lead(p5a_hat))/2*5)   
  
  # find ex for last age
  last_age <- max(surv_data$age)
  obs_data_last_age <- data.frame(age  = last_age, p5last_age = regression_data$p5a_hat[regression_data$age==last_age])
  mlt_data_last_age <- MortCast::MLTlookup %>% 
    dplyr::filter(type == mlt_family, sex == this_sex) %>% 
    dplyr::group_by(type, e0) %>% 
    dplyr::summarise(elast_age = sum(Lx[age>=last_age])/lx[age==last_age],
                     p5last_age = lx[age==last_age]/lx[age==5], .groups = 'drop') %>% 
    dplyr::mutate(age = last_age) %>% 
    dplyr::select(-e0)
  mlt_closest <- interp_level_mlt(obs_data_last_age, mlt_data_last_age, 
                                  var_interp = "elast_age", var_ref = "p5last_age")
  
  # compute e5 and return
  e5 <- sum(regression_data$nLx[regression_data$age<last_age]) + regression_data$p5a_hat[regression_data$age==last_age] * mlt_closest$elast_age_interp
  return(list(e5 = e5, model = summary(model), plot = plot_obs_fit))
}


# level interpolation for any lt function ------------------------------------

# TR: if this is just picking out two patterns for different e0 levels.
# Can this be done with less dplyr verbage? Can it make use of DemoTools interp functions
# that already exist? or even approx() from base?
interp_level_mlt <- function(obs_data, mlt_data, var_interp = "e0", var_ref = "nSx"){
  
  # match var names
  if(!all(colnames(obs_data) %in% colnames(mlt_data))) stop("names in obs_data are not in mlt_data")
  match.arg(var_interp, colnames(mlt_data))
  
  # rename
  var_reference_index_obs <- which(colnames(obs_data)==var_ref)
  colnames(obs_data)[var_reference_index_obs] <- "ref"
  var_reference_index_mlt <- which(colnames(mlt_data)==var_ref)
  colnames(mlt_data)[var_reference_index_mlt] <- "var_ref"
  var_interp_index_mlt <- which(colnames(mlt_data)==var_interp)
  colnames(mlt_data)[var_interp_index_mlt] <- "var_int"
  
  # join and find closest
  data <- dplyr::inner_join(obs_data, mlt_data, 
                            by = if(!"type" %in% colnames(obs_data)) "age" else {c("type", "age")}) %>%
    dplyr::mutate(diff = abs(var_ref/ref-1)) %>% 
    dplyr::arrange(type, age, var_int) %>% 
    dplyr::group_by(type, age) %>% 
    dplyr::arrange(diff) %>% 
    dplyr::slice(1:2) %>%
    dplyr::mutate(levels = c("left", "right")) %>% 
    tidyr::pivot_wider(id_cols = -diff, names_from = levels, values_from = c(var_int, var_ref)) %>% 
    dplyr::mutate(diff = min(abs(ref - var_ref_left), abs(ref - var_ref_right)),
                  interp = var_int_left + (ref - var_ref_left)/(var_ref_right - var_ref_left) * (var_int_right - var_int_left)) %>% 
    dplyr::select(type, age, ref, var_ref_left, var_ref_right, var_int_left, var_int_right, interp) %>% 
    ungroup()
  
  # rename
  colnames(data) <- c("type", "age", paste0(var_ref,"_obs"), 
                      paste0(var_ref,"_left"), paste0(var_ref,"_right"), 
                      paste0(var_interp,"_left"), paste0(var_interp,"_right"), 
                      paste0(var_interp,"_interp"))
  return(data)
}


# Feeney´s method (UN, 2002) ---------------------------------------------

intercensal_surv_feeney <- function(c1,
                                    c2,
                                    interc_t = 10,
                                    age = NULL,
                                    age_int = NULL,
                                    TOAGless5_l2.5 = NULL){
  
  # arrange data
  census_data <- data.frame(c1, c2, age, age_int, age_mp = age + age_int/2)
  
  # remove OAG because DemoTools function already assume that
  age_OAG <- age[age_int %in% c(0, -1, NA)]
  surv_data <- census_data[census_data$age!=age_OAG,]
  
  # apply method
  if(interc_t == 5){
    # TR: These functions are currently in CENSUR.R
    # Possibly we should be grouping functionality differently now.
    surv_data$lx_l2.5 <- surv5(surv_data$c1, surv_data$c2)
  }
  if(interc_t == 10){
    surv_data$lx_l2.5 <- surv10(surv_data$c1, surv_data$c2)
  }
  if(!interc_t %in% c(5,10)){
    surv_data$lx_l2.5 <- survN(surv_data$c1, surv_data$c2, interval = interc_t)
  }
  
  # complete Feeney development
  surv_data <- surv_data %>% 
    mutate(lx_l2.5 = (lx_l2.5+lag(lx_l2.5))/2,
           Lx5_l2.5 = (lx_l2.5 + lead(lx_l2.5))*age_int/2)
  
  # depending TOAGless5_l2.5 is given, can return full e_x or temporary until OAG-5
  if(is.null(TOAGless5_l2.5)){
    surv_data <- surv_data %>% 
      # TR use lt_id_L_T()
      mutate(Tx_l2.5 = purrr::accumulate(Lx5_l2.5, sum, .dir = "backward", na.rm=T),
             Tx_l2.5 = ifelse(age == 0, NA, Tx_l2.5),
             ex_temp = Tx_l2.5/lx_l2.5)
  }else{
    surv_data <- surv_data %>% 
      mutate(Lx5_l2.5 = ifelse(age == last(age), TOAGless5_l2.5, Lx5_l2.5),
             Tx_l2.5 = purrr::accumulate(Lx5_l2.5, sum, .dir = "backward", na.rm=T),
             Tx_l2.5 = ifelse(age == 0, NA, Tx_l2.5),
             ex = Tx_l2.5/lx_l2.5)
  }
  return(surv_data[!is.na(surv_data$Lx5_l2.5),])
}

# variable r method --------------------------------------------------------------

# TODO: document
intercensal_surv_var_r <- function(c1,
                                   c2,
                                   date1,
                                   date2,
                                   age1 = 1:length(c1) - 1,
                                   age2 = 1:length(c2) - 1,
                                   sex = "f"){
  # manage ages: fit to more wider group
  age1_int <- max(unique(diff(age1)))
  age2_int <- max(unique(diff(age2)))
  age_int <- max(age1_int, age2_int)  
  # TR: check for DemoTools function that does this: groupAges()
  c1 <- data.frame(c1=c1, age1=trunc(age1/age_int)*age_int) %>% group_by(age1) %>% summarise(c1 = sum(c1)) %>% pull(c1)
  c2 <- data.frame(c2=c2, age2=trunc(age2/age_int)*age_int) %>% group_by(age2) %>% summarise(c2 = sum(c2)) %>% pull(c2)
  age <- seq(0, min(max(age1), max(age2)), age_int)
  OAG <- max(age)
  ages <- length(age)
  # TR: lubridate seems lightweight, we might as well replace all instances of dec.date with the lubridate version
  # and deprecate our own version.
  if(!is.numeric(date1)) date1 <- round(dec.date(date1),2)
  if(!is.numeric(date2)) date2 <- round(dec.date(date2),2)
  # if((date2-date1)>15) stop("more than 15 years in intercensal period")
  
  # apply method
  lt_ManualX_variable_r(Age = age, Nx1 = c1, Nx2 = c2, ts = date2-date1, 
                        radix = 100000, full_lt = T)
}

# TODO: make this legible, document, test.
# from Manual X, r-var method
# from Michael Lachanski (mikelach@sas.upenn.edu)
# be sure to add Preston reference, or the exact manual X reference, or something
#' @author author
lt_ManualX_variable_r <- function(age, Nx1, Nx2, ts, radix = 1000, full_lt = T){
  DT <- data.table(age, Nx1, Nx2, ts, radix, key = "age")
  
  # set rho parameters - taken from Table 185 on page 219 in Manual X.
  
  if (max(age) >= 85)     { a = 0.006; b =  2.68;  c = 0.006 }
  else if (max(age) >= 80){ a = 0.025; b =  4.30;  c = 0.029 }
  else if (max(age) >= 75){ a = 0.053; b =  6.40;  c = 0.063 }
  else if (max(age) >= 70){ a = 0.086; b =  8.77;  c = 0.102 }
  else if (max(age) >= 65){ a = 0.119; b = 11.22;  c = 0.141 }
  else if (max(age) >= 60){ a = 0.150; b = 13.66;  c = 0.176 }
  else if (max(age) >= 55){ a = 0.179; b = 16.02;  c = 0.207 }
  else if (max(age) >= 50){ a = 0.205; b = 18.28;  c = 0.235 }
  else if (max(age)  < 50){ a = 0.229; b = 20.43;  c = 0.258 }
  
  DT_lt <-
    DT %>%
    # easier to just keep track of two n's follow Preston et al. (2001), p.186 - 187
    .[ ,  `:=`
       (
         # Step 1: adjustment for net intercensal migration and territorial coverage
         # Step 2: calculation of age-specific intercensal growth rates
         # This function computes the intercensal age-specific growth rate.
         n_r_x  = log(Nx2/Nx1)/ts,
         # Step 3: calculation of average intercensal age distribution.
         nNx_mid = (Nx1 + Nx2)/2,
         # The above function takes the mean of exactly two scalars.
         n_before = age - shift(age, type = "lag", n = 1, fill = 0),
         n_after = shift(age, type = "lead", n = 1, fill = 0) - age
       )
    ] %>%
    .[age >= 10, `:=` (N1_10_plus = sum(Nx1), N2_10_plus = sum(Nx2))] %>%
    .[age >= 45, `:=` (N1_45_plus = sum(Nx1), N2_45_plus = sum(Nx2))] %>%
    # make variables used for Step 4
    .[ , `:=` (nNx_mid_10plus = 0.5*(N1_10_plus + N2_10_plus),
               nNx_mid_45plus = 0.5*(N1_45_plus + N2_45_plus),
               r10plus = log(N2_10_plus/N1_10_plus)/ts)
    ] %>%
    # Preston deals with the open interval by assuming that there is an age or duration
    # equidistant from the last two durations. At the final duration, everyone dies.
    .[nrow(DT) , n_after  := n_before]  %>%
    # Step 4: cumulation of age-specific growth rates from age 5 upward.
    .[ , Rx := 0.5*n_after*n_r_x + cumsum(n_before*shift(n_r_x, fill= 0))]
  
  if(full_lt == F){
    DT_lt[nrow(DT),  n_r_x := NA_real_] %>%
      .[ ,
         Rx := 0.5*n_before*replace_na(n_r_x, 0) +
           cumsum(shift(n_before, fill = 0)*shift(n_r_x, fill= 0))]   %>%
      # only in Preston (1983) does he do this.
      .[1, Rx :=NA_real_]
  }
  DT_lt %>%
    # Accumulation function: Rx is what Preston calls it in Manual X
    .[ , rho := a + b*r10plus + c*log(nNx_mid_45plus/nNx_mid_10plus)] %>%
    .[nrow(DT) , Rx := Rx  + rho] %>%
    # Step 5: reduction of age distribution to a stationary form.
    .[ ,  Lx_star := nNx_mid*exp(Rx)] %>%
    # Step 6: calculation of the expectation of life.
    .[               , lx := 1/(n_before + n_after)*(Lx_star + shift(Lx_star))] %>%
    .[               , temp_lx := shift(lx, n = 1, type = "lead")] %>%
    .[1              , lx := fifelse(radix > temp_lx, radix, temp_lx)] %>%
    .[1              , lx := fifelse(full_lt == F, NA_real_, lx)] %>%
    .[nrow(DT)       , lx := fifelse(full_lt == F, NA_real_, lx)] %>%
    .[  , temp_lx := NULL] %>%
    # DemoTools only cares about some columns of lifetable based on example
    # https://timriffe.github.io/DemoTools/articles/lifetables_with_demotools.html
    .[         , ndx := lx - shift(lx, n = 1, type = "lead")] %>%
    .[nrow(DT) , ndx := lx ]   %>%
    .[         , Tx := revcumsumSkipNA(Lx_star)] %>%
    .[         , Sx := lx/lx[1]] %>%
    .[         , ex := Tx/lx]
  if(full_lt == F){DT_lt[age > 50 , ex := NA_real_] }
  
  
  # Once the life expectancy figures have been calculated, usually for x ranging
  # from 10 to 50, the levels they imply in a model life-table system can be found,
  # and a final estimate of mortality can be obtained by averaging the most reliable
  # estimates of mortality level (those left after discarding any clearly unsuitable
  # values). In practice, the mortality estimates for values of x up to age 30 or
  # so are reasonably consistent, but after age 30 or 35 they often show progressively
  # lower mortality as age x increases. The best estimate of overall mortality may
  # therefore be an average of the levels associated with ex or x ranging from 10
  # to 30, though this conclusion implies that the results will not be a useful
  # basis for the selection of an age pattern of mortality, nor will they be good
  # indicators of the necessity of adjustment when errors in the growth rates arise
  # because of changes in enumeration completeness.
  DT_lt %>% `[`
  return(DT_lt)
}
