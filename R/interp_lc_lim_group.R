#' Lee-Carter method with limited data for groups.

# tests:
# Tests againts spreadsheet
# test output againts `interp_lc_lim` function.
# testing args from main fun
# mixing input: single/abr with output single/abr, and mixing input nMx and lx
# passing lt arguments
# text/messages/warnings. Specially the case when no `id` is given

#' 
#' @description Given a data frame with groups (country, region, sex), dates, sex and mortality data by age (rates, conditionated probabilities of death 
#' or survival function), this function interpolate/extrapolate life tables 
#' using the method for limited data suggested by  Li et. al (2004) (at least three observed years). 
#'
#' @details Based on spreedsheet "Li_2018_Limited_Lee-Carter-v4.xlsm" from UN. 
#' Useful for abridged or single ages, and allows output in both formats also.
#' One option is the use of non-divergent method for sex coherency (Li & Lee, 2005).
#' The other is the possibility of fitting `"k"` to replicate `"e_0"` at some given dates.
#' `id` column in `input` argument works for separate between groups. In case only one population/sex is given, 
#' is recommended to give some group name to `id`, if not the function will try to infer the case. 
#'
#' @note Draft Version
#'
#' @param input data.frame. Columns: id, Date, Sex, Age, nMx (opt), nqx (opt), lx (opt). 
#' The first column (id) cn be a numeric index or charcter vector identifying each group.  
#' @param dates_out numeric. Vector of decimal years to interpolate or extrapolate.
#' @param Single logical. Wheter or not the lifetable output is by single ages.
#' @param input_e0 data.frame with cols: id, Date, Sex and `"e_0"`. This should be fitted when apply method.
#' @param prev_divergence logical. Whether or not prevent divergence and sex crossover between groups. Default `FALSE.`
#' @param weights list. For `prev_divergence` option. A double for each element of a list with names as `id` columns. Should sum up to 1. Default: same weight for each group.
#' @param OAG logical. Whether or not the last element of `nMx` (or `nqx` or `lx`) is an open age group. Default `TRUE.`
#' @param verbose logical. Default `FALSE`.
#' @param SVD logical. Use Singular Value Decomposition for estimate b and k or Maximum Likelihood Estimation. Default `FALSE` for Maximum Likelihood Estimation.
#' @param ... Other arguments to be passed on to the \code{\link[DemoTools]{lt_abridged}} function.
#' @seealso
#' \code{\link[DemoTools]{lt_abridged}}
#' @export
#' @importFrom data.table rbindlist
#' @importFrom data.table setDT
#' @importFrom data.table uniqueN
#' @return List with:
#' \itemize{
#'   \item Lifetable in a data.frame with columns:
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
#'   \item List with parameters estimated for each group:
#'   * `kt` numeric time vector. Time trend in mortality level.
#'   * `ax` numeric age vector. Average time of `log(m_{x,t})`.
#'   * `bx` numeric age vector. Pattern of change in response to `kt`.
#' }
#' @references
#' \insertRef{Li2005}{DemoTools}
#' \insertRef{Li2004}{DemoTools}
#'
#' @examples
#' # mortality rates from Sweden, for specific dates. Each sex a group.
#' mA_swe$id = c(rep("A",nrow(mA_swe)/2),
#'              rep("B",nrow(mA_swe)/2))
#' 
#' # needs mortality rates in this dates: 
#' dates_out <- as.Date(paste0(seq(1948,2018,5),"-07-01"))
#' 
#' # apply LC with limited data to extrap/interpolate
#' lc_lim_data <- interp_lc_lim_group(input = mA_swe, dates_out = dates_out)
#' 
#' \dontrun{
#' lc_lim_data[["lt_hat"]] %>% ggplot(aes(Age,nMx,col=factor(round(Date,1)))) +
#'   geom_step() + scale_color_viridis_d() + 
#'   scale_y_log10() + theme_classic() + facet_wrap(~Sex)
#' }
#' 
#' # avoid cross-over between groups
#' lc_lim_data <- interp_lc_lim_group(input = mA_swe, dates_out = dates_out,
#'                                     prev_divergence = TRUE, weights=list(A=.4,B=.6))
#' 
#' \dontrun{
#' lc_lim_data[["lt_hat"]] %>% ggplot(aes(Age,nMx,col=factor(round(Date,1)))) +
#'   geom_step() + scale_color_viridis_d() + 
#'   scale_y_log10() + theme_classic() + facet_wrap(~id)
#' }

# fun ---------------------------------------------------------------------

interp_lc_lim_group <- function(input = NULL,
                                dates_out = NULL, 
                                Single = FALSE,
                                input_e0 = NULL, 
                                prev_divergence = FALSE,
                                weights = NULL,
                                OAG = TRUE,
                                verbose = TRUE,
                                SVD = FALSE,
                                ...){
  
  # just make up empty placeholders for stuff used inside data.table
  .    <- NULL
  Date <- NULL
  id   <- NULL
  
  ExtraArgs = c(as.list(environment()), list(...))
  ExtraArgs = ExtraArgs[! names(ExtraArgs) %in% c("input","Single")]
  dates_out <- dec.date(dates_out)
  ndates_out <- length(dates_out)
  
  # enough obs dates
  min_dates_in <- min(setDT(input)[, .(count = uniqueN(Date)), 
                                               by = id][,2])
  if (min_dates_in<3){
    stop("\nYou need three observed dates at least.")
  }
  
  # is data limited?
  diff_dates_in <- unique(setDT(input)[, .(count = diff.Date(Date)), 
                                                   by = id][,2])
  if (all(diff_dates_in %in% c(0,1))){
    stop("\nYou have single-year-interval data and probably should use basic Lee-Carter method.")
  }
  
  # check input
  if (!any(names(input)%in%c("nMx", "nqx", "lx"))){
    stop("\nSorry we need some column called nMx, nqx or lx\n")
  }
  
  # TR: I commented this out. 
  # 1) id is used earlier than this, so there's still an error if it's missing
  # 2)
  # # you gave no id - save it
  # if (!"id" %in% colnames(input)){
  #   # but two sex
  #   cases <- aggregate(Age~Date+Sex,input,FUN=length)
  #   if(!any(cases$Age)==cases$Age[1]){
  #     input$id = ifelse(Sex=="f",1,2)
  #   }
  # }
  
  ngroups <- length(unique(input$id))
  groups  <- unique(input$id)
  # three objects, with number of elements as groups
  nMx         <- list()
  nMx_hat     <- list()
  lc_estimate <- list()
  
  # get always Mx -----------------------------------------------------------
  . <- NULL
  for(i in groups){
    input_id <- as.data.frame(input[input$id == i,])
      # inputdt <- split(input_id, list(input_id$Date)) %>% 
      # lapply(
      #   function(X) do.call(lt_smooth_ambiguous,
      #                       c(list(input=X), ExtraArgs))) %>% 
      # do.call("rbind", .)%>% 
      # as.data.table()%>% 
      # data.table::dcast(Age ~ Date, value.var = "nMx") %>% 
      # .[order(Age)]
      # 
      # 
      input_id <- as.data.frame(input[input$id == i,])
      inputdt <- split(input_id, list(input_id$Date)) %>% 
        lapply(., function(X) {
          Age <- X$Age
          Sex_i <- unique(X$Sex)
          
          
          types     <- c("nMx","nqx","lx")
          this_type <- types[types %in% colnames(X)]
          if (length(this_type) > 1){
            ind <- X[,this_type] %>% 
              as.matrix() %>% 
              is.na() %>% 
              colSums() %>% 
              which.min()
            this_type <- this_type[ind]
          }
      
          LT <- lt_ambiguous(nMx_or_nqx_or_lx = X[[this_type]],
                             type = this_type,
                             Age = Age, 
                             Sex = Sex_i,
                             Single = Single,
                             ...)
         
          LT$Date <- unique(X$Date)
          LT
        }) %>% 
        do.call("rbind", .)%>% 
        as.data.table()%>% 
        data.table::dcast(Age ~ Date, value.var = "nMx") %>% 
        .[order(Age)]
      
      
    Age  <- inputdt[["Age"]]
    nMx[[i]] <- inputdt[, -1] %>% as.matrix()
    rownames(nMx[[i]]) <- Age
    nAge       <- length(Age)
    
    # LC and estimate
    dates_in  <- unique(input_id$Date) %>% dec.date()
    lc_estimate[[i]] <- interp_lc_lim_estimate(nMx[[i]], dates_in, dates_out, SVD)
    nMx_hat[[i]] <- exp(lc_estimate[[i]][["ax"]] + lc_estimate[[i]][["bx"]] %*% t(lc_estimate[[i]][["kt"]]))
  }
  
  # options -----------------------------------------------------------------
  
  # prevent divergence/cross-over
  if (prev_divergence){
    if(ngroups==1) stop("No subgroups no divergence.")
    # weigths
    if(is.null(weights)){
      weights <- list()
      for(i in 1:ngroups) weights[[i]] = 1/ngroups
    }else{
      if(sum(unlist(weights))!=1) stop("Weights do not sum up to 1.")
    }
    # weighted mean of parameters
    bx_div <- rep(0,nrow(lc_estimate[[1]][["bx"]]))
    kt_div <- rep(0,length(lc_estimate[[1]][["kt"]]))
    k0_div <- 0
    for(i in 1:ngroups){
      # i =1
      bx_div = bx_div + lc_estimate[[i]][["bx"]] * weights[[i]]
      kt_div = kt_div + lc_estimate[[i]][["kt"]] * weights[[i]]
      k0_div = k0_div + lc_estimate[[i]][["k0"]] * weights[[i]]
    }
    for(i in 1:ngroups){
    lc_estimate[[i]][["bx"]] <- bx_div
    lc_estimate[[i]][["kt"]] <- kt_div
    lc_estimate[[i]][["k0"]] <- k0_div
    }
  }
  
  # fit e_0 and/or prev_divergence
  for(i in groups){
    ax_i <- lc_estimate[[i]][["ax"]]
    bx_i <- lc_estimate[[i]][["bx"]]
    kt_i <- lc_estimate[[i]][["kt"]]
    k0_i <- lc_estimate[[i]][["k0"]]
    e0 = input_e0[input_e0$id==i,"e0"]
    dates_in  <- unique(input$Date[input$id == i]) %>% dec.date()
    
    if (!is.null(e0)){
      dates_e0 = input_e0[input_e0$id==i,"Date"]
      ndates_e0 = length(dates_e0)
      Sex_e0 = unique(input_e0[input_e0$id==i,"Sex"])
      e0_star <- interp(rbind(e0, e0), 
                        dates_e0,dates_out,
                        extrap = TRUE)[1, ]
      kt_star = c()
      for (j in 1:ndates_out){
        kt_star[j] <- optimize(f = interp_lc_lim_kt_min,
                                interval = c(-20, 20),
                                ax = ax_i,
                                bx = bx_i,
                                age = Age,
                                sex = Sex_e0,
                                e0_target = e0_star[j],
                                ...)$minimum
      }
      nMx_hat[[i]] <- exp(ax_i + bx_i %*% t(kt_star))
    }else{
      if(prev_divergence){
        nMx_hat_div <- nMx[[i]][,1] * exp(bx_i %*% t(kt_i-k0_i))
        dates_extrap <- dates_out < min(dates_in)
        nMx_hat[[i]][,dates_extrap] <- nMx_hat_div[,dates_extrap]
      }
    }
   
    # return lt
    colnames(nMx_hat[[i]]) <- dates_out
    Sex_i = unique(input$Sex[input$id == i])
    out <-
      lapply(colnames(nMx_hat[[i]]), function(xx,MX,Age) {
        mx <- MX[, xx]
        LT <- lt_ambiguous(nMx_or_nqx_or_lx = mx,
                           type = "m",
                           Age = Age, 
                           Sex = Sex_i,
                           Single = Single,
                           ...)
        LT$Sex  <- Sex_i
        LT$Date <- as.numeric(xx)
        LT
      }, MX = nMx_hat[[i]], Age = Age) %>% 
      rbindlist()
    nMx_hat[[i]] <- out
  }

  return(list(
            lt_hat = rbindlist(nMx_hat, idcol = "id"),
            lc_params = lc_estimate #IW: must bind

            ))
}

