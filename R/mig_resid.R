
# TODO 
# This is a high priority
# -[ ] make sure mig_resid_cohort() handles dimensions properly (named indexing; no waste dims)
# -[ ] make sure mig_resid_time() handles dimensions properly
# -[ ] check dims of incoming arguments.
# -[ ] new args years_pop, years_asfr, years_sr, years_srb (to be fed to checker)
# -[ ] write a dimension checker + trimming mig_resid_dim_check()
# -[ ] make this checker/trimmer the first step in mig_resid*()

# This can come next
# -[ ] make new package data. usethis::use_data(pop_m_mat)
# -[ ] document new package data in data.R following other examples

# Then this
# -[ ] write wrapper function, mig_resid() with an argumet 'method' 
#      with options "cohort", "stock" or "time", and all other args the same.

# Then this
# -[ ] unit tests

# Then this
# -[ ] sanity checks: do estimated migration patterns actually look reasonable in 
#      periods/places that are known to be strong in or out migration places.



#' Estimate net migration using residual methods: stock change,
#' time even flow and cohort even flow
#'
#' @details
#'
#' 1. The stock method (\code{mig_resid_stock}) is the difference in stocks that
#' survive between t and t+5, and the first age group is based on the difference
#' with the surviving births by sex.  It provides net migrants by lexis cohort
#' parallelograms, and basically such info gets used as end-period migration
#' since the migrants don't get exposed to mortality within the period.
#'
#' 2. The time even flow (\code{mig_resid_time}) method uses the result from
#' the first option, but splits it back into lexis period squares and assumes
#' that half of the net migrants get exposed to the mortality risk during this
#' period. Such info can get used as evenly distributed migration by period,
#' but the assumptions lead to zig-zag age patterns that are highly implausible.
#'
#' 3. The cohort even flow (\code{mig_resid_cohort}) method provides the most
#' meaningful pattern of net migration by age consistent by cohort and assumes
#' an evenly distribution within the 5-year period, and half of the migrants
#' get exposed both fertility and mortality within this period.
#'
#' @param pop_m_mat A \code{numeric} matrix with population counts. Rows should
#' be ages and columns should be years. Only five year age groups are supported.
#' See examples.
#'
#' @param pop_f_mat A \code{numeric} matrix with population counts. Rows should
#' be ages and columns should be years. Only five year age groups are supported.
#' See examples.
#'
#' @param sr_m_mat A \code{numeric} matrix with survival rates for males. Rows
#' should be ages and columns should be years. ** This matrix should have
#' one column less than \code{pop_m_mat} and \code{pop_f_mat}. For example,
#' if the last year in these matrices is 2050, then the last year in
#' \code{sr_m_mat} should be 2045. **
#'
#' @param sr_f_mat A \code{numeric} matrix with survival rates for females. Rows
#' should be ages and columns should be years. ** This matrix should have
#' one column less than \code{pop_m_mat} and \code{pop_f_mat}. For example,
#' if the last year in these matrices is 2050, then the last year in
#' \code{sr_f_mat} should be 2045. **.
#'
#' @param asfr_mat A \code{numeric} matrix with age specific fertility rates.
#' Rows should be ages and columns should be years. ** This matrix should have
#' one column less than \code{pop_m_mat} and \code{pop_f_mat}. For example,
#' if the last year in these matrices is 2050, then the last year in
#' \code{asfr_mat} should be 2045**. This row will usually have fewer age groups
#' (rows) than in the population matrices or survival matrices, so the user
#' needs to supply the specific ages in the \code{ages_fertility} argument.
#'
#' @param srb_vec A \code{numeric} vector of sex ratios at birth for every year.
#' The years should be the same as the years in \code{sr_m_mat},
#' \code{sr_f_mat}, and \code{asfr_mat}.
#'
#' @param ages A \code{numeric} vector of ages used in the rows in
#' \code{pop_m_mat}, \code{pop_f_mat}, \code{sr_m_mat}, \code{sr_f_mat}.
#'
#' @param ages_fertility A \code{numeric} vector of ages used in the rows in
#' \code{asfr_mat}.
#'
#' @return A list with two matrices. One is for males (called `mig_m`) and the
#' other for females (called `mig_f`). Both matrices contain net migration
#' estimates by age/period using one of the three methods.
#'
#' @examples
#'
#' ################ Stock change method #####################
#'
#' mig_res <-
#'  mig_resid_stock(
#'    pop_m_mat = pop_m_mat,
#'     pop_f_mat = pop_f_mat,
#'    sr_m_mat = sr_m_mat,
#'    sr_f_mat = sr_f_mat,
#'    asfr_mat = asfr_mat,
#'    srb_vec = srb_vec,
#'     ages = ages,
#'    ages_fertility = ages_fertility
#'  )
#'
#' # Net migration for males using stock change method
#' mig_res$mig_m
#'
#' # Net migration for females using stock change method
#' mig_res$mig_f
#'
#'
#' ################ cohort even flow method  #####################
#'
#' # We reuse the same data from before
#'
#' mig_res <-
#'   mig_resid_cohort(
#'     pop_m_mat = pop_m_mat,
#'     pop_f_mat = pop_f_mat,
#'     sr_m_mat = sr_m_mat,
#'     sr_f_mat = sr_f_mat,
#'     asfr_mat = asfr_mat,
#'     srb_vec = srb_vec,
#'     ages = ages,
#'     ages_fertility = ages_fertility
#'   )
#'
#' # Net migration for males using the cohort even flow method
#' mig_res$mig_m
#'
#' # Net migration for females using the cohort even flow method
#' mig_res$mig_f
#'
#' ################ time even flow method  #####################
#'
#' # We reuse the same data from before
#'
#' mig_res <-
#'   mig_resid_time(
#'     pop_m_mat = pop_m_mat,
#'     pop_f_mat = pop_f_mat,
#'     sr_m_mat = sr_m_mat,
#'     sr_f_mat = sr_f_mat,
#'     asfr_mat = asfr_mat,
#'     srb_vec = srb_vec,
#'     ages = ages,
#'     ages_fertility = ages_fertility
#'   )
#'
#' # Net migration for males using the time even flow method
#' mig_res$mig_m
#'
#' # Net migration for females using the time even flow method
#' mig_res$mig_f
#'
#' @export
mig_resid_stock <- function(pop_m_mat,
                            pop_f_mat,
                            sr_m_mat,
                            sr_f_mat,
                            asfr_mat,
                            srb_vec,
                            ages = NULL,
                            ages_fertility = NULL,
                            years_pop = NULL,
                            years_sr = NULL,
                            years_asfr = NULL,
                            years_srb = NULL) {
  # this arg list can feed into the checker
  args_list <- as.list(match.call())

  mig_resid_dim_checker(args_list)
  
  pop_m_mat <- args_list$pop_m_mat
  pop_f_mat <- args_list$pop_f_mat
  sr_m_mat  <- args_list$sr_m_mat
  sr_f_mat  <- args_list$sr_f_mat
  asfr_mat  <- args_list$asfr_mat
  srb_vec   <- args_list$srb_vec



  stopifnot(
    is.matrix(pop_m_mat),
    is.matrix(pop_f_mat),
    is.matrix(sr_m_mat),
    is.matrix(sr_f_mat),
    is.matrix(asfr_mat),
    is.numeric(srb_vec),
    is.numeric(ages),
    is.numeric(ages_fertility)
  )

  

# <<<<<<< HEAD
#   # Check in dimensions are ok - still working on this
#   if(ncol(asfr_mat) == ncol(pop_f_mat) -1 & nrow(sr_f_mat) == nrow(pop_f_mat) -1){
#     print("matrix dimensions are correct")
#   }
#     else {
#     print("check matrix dimensions")
#   }
#   
#   #if there are extra years, drop it - still thinking the best way to deal with it
#   if(ncols(asfr_mat) != ncols(sr_f_mat)){
#     asfr_mat <- asfr_mat[, colnames(sr_f_mat)]
#     sr_f_mat <- sr_f_mat[, colnames(asfr_mat)]
#   }
#   else {
#     asfr_mat
#     sr_f_mat
#   }
# =======
# >>>>>>> 362ae9857574b05c519c8de40548461d6b9070dd



  # Migration net of only survivors
  net_mig_m <- migresid_net_surv(pop_m_mat, sr_m_mat)
  net_mig_f <- migresid_net_surv(pop_f_mat, sr_f_mat)

 # fertility_index <- which(ages %in% ages_fertility)

  # Returns all births for all years
  age_interval <- unique(diff(ages))
  all_births <- migresid_births(
    pop_f_mat,
    asfr_mat,
   # fertility_index,
    age_interval
  )

  # With all_births already calculated, separate between
  # female/male births with the sex ratio at birth
  byrs     <- names(all_births)
  births_m <- all_births * (srb_vec[byrs] / (1 + srb_vec[byrs]))
  births_f <- all_births * (1 / (1 + srb_vec[byrs]))

  net_mig_m <- migresid_net_surv_first_ageg(
    net_mig_m,
    pop_m_mat,
    births_m,
    sr_m_mat
  )

  net_mig_f <- migresid_net_surv_first_ageg(
    net_mig_f,
    pop_f_mat,
    births_f,
    sr_f_mat
  )

  # First year is empty, so we exclude
  list(
    mig_m = net_mig_m,
    mig_f = net_mig_f
  )
}

#' @rdname mig_resid_stock
#' @export
mig_resid_cohort <- function(pop_m_mat,
                             pop_f_mat,
                             sr_m_mat,
                             sr_f_mat,
                             asfr_mat,
                             srb_vec,
                             ages = NULL,
                             ages_fertility = NULL,
                             years_pop = NULL,
                             years_sr = NULL,
                             years_asfr = NULL,
                             years_srb = NULL) {
  # this arg list can feed into the checker
  args_list <- as.list(match.call())
  
  mig_resid_dim_checker(args_list)
  
  pop_m_mat <- args_list$pop_m_mat
  pop_f_mat <- args_list$pop_f_mat
  sr_m_mat  <- args_list$sr_m_mat
  sr_f_mat  <- args_list$sr_f_mat
  asfr_mat  <- args_list$asfr_mat
  srb_vec   <- args_list$srb_vec

  # Estimate stock method
  mig_res <-
    mig_resid_stock(
      pop_m_mat = pop_m_mat,
      pop_f_mat = pop_f_mat,
      sr_m_mat = sr_m_mat,
      sr_f_mat = sr_f_mat,
      asfr_mat = asfr_mat,
      srb_vec = srb_vec,
      ages = ages,
      ages_fertility = ages_fertility
    )

  net_mig_m <- mig_res$mig_m
  net_mig_f <- mig_res$mig_f

  # Estimate bounds for males
  mig_m_bounds <- migresid_bounds(net_mig_m, sr_m_mat)
  mig_upper_m <- mig_m_bounds$upper
  mig_lower_m <- mig_m_bounds$lower

  # Estimate bounds for females
  mig_f_bounds <- migresid_bounds(net_mig_f, sr_f_mat)
  mig_upper_f <- mig_f_bounds$upper
  mig_lower_f <- mig_f_bounds$lower

  # Adjust last age group in the bounds
  mig_bounds <- migresid_bounds_last_ageg(
    net_mig_m,
    net_mig_f,
    mig_upper_m,
    mig_lower_m,
    mig_upper_f,
    mig_lower_f
  )

  mig_upper_m <- mig_bounds$mig_upper_m
  mig_lower_m <- mig_bounds$mig_lower_m
  mig_upper_f <- mig_bounds$mig_upper_f
  mig_lower_f <- mig_bounds$mig_lower_f

  # Combine both upper/lower bound into a single rectangle
  mig_rectangle_m <- mig_upper_m + mig_lower_m
  mig_rectangle_f <- mig_upper_f + mig_lower_f

 list(
   mig_m = mig_rectangle_m[, -1],
   mig_f = mig_rectangle_f[, -1]
 )
 # TR: we prefer this, but somewhere earlier in processing
 # we get an extra column 1 full of NAs. So let's just not let that happen
  # list(
  #   mig_m = mig_rectangle_m,
  #   mig_f = mig_rectangle_f
  # )
}

#' @rdname mig_resid_stock
#' @export
mig_resid_time <- function(pop_m_mat,
                           pop_f_mat,
                           sr_m_mat,
                           sr_f_mat,
                           asfr_mat,
                           srb_vec,
                           ages = NULL,
                           ages_fertility = NULL,
                           years_pop = NULL,
                           years_sr = NULL,
                           years_asfr = NULL,
                           years_srb = NULL) {
  # this arg list can feed into the checker
  args_list <- as.list(match.call())
  
  mig_resid_dim_checker(args_list)
  
  pop_m_mat <- args_list$pop_m_mat
  pop_f_mat <- args_list$pop_f_mat
  sr_m_mat  <- args_list$sr_m_mat
  sr_f_mat  <- args_list$sr_f_mat
  asfr_mat  <- args_list$asfr_mat
  srb_vec   <- args_list$srb_vec

  # TR: add chunk (maybe a new function?) that
  # checks dimensions; names dimensions if necessary (and warns if so)
  # and trims dimensions if necessary (warning user if needed).
  # warning does mean warning() it just means cat("\nwatch out!\n")
  # Not important, but theese could be silenced using a new 'verbose' arg
  
  
  # Estimate stock method
  mig_res <-
    mig_resid_stock(
      pop_m_mat = pop_m_mat,
      pop_f_mat = pop_f_mat,
      sr_m_mat = sr_m_mat,
      sr_f_mat = sr_f_mat,
      asfr_mat = asfr_mat,
      srb_vec = srb_vec,
      ages = ages,
      ages_fertility = ages_fertility
    )

  # Separate male/female net migration
  net_mig_m <- mig_res$mig_m
  net_mig_f <- mig_res$mig_f

  # Adjust age group 0-4
  net_mig_m[1, ] <- 2 * net_mig_m[1, ]
  net_mig_f[1, ] <- 2 * net_mig_f[1, ]

  # Adjust age groups 5-10  to 100+ (of whatever maximum age groups)
  for (i in 2:nrow(net_mig_m)) {
    double_pop_m <- (2 * net_mig_m[i, ])
    double_pop_f <- (2 * net_mig_f[i, ])

    # Multiply net mig of i - 1 by survival rate of i
    # to get number of survived
    mig_sr_m <- net_mig_m[i - 1, ] * sr_m_mat[i, ]
    mig_sr_f <- net_mig_f[i - 1, ] * sr_f_mat[i, ]

    net_mig_m[i, ] <- double_pop_m - mig_sr_m
    net_mig_f[i, ] <- double_pop_f - mig_sr_f
  }

  list(
    mig_m = net_mig_m,
    mig_f = net_mig_f
  )
}


# Net migration is pop minus the people that survived from the previous
# age/cohort
migresid_net_surv <- function(pop_mat, sr_mat) {
  n                <- nrow(pop_mat)
  p                <- ncol(pop_mat)
  survived         <- pop_mat[-n, -p] * sr_mat[-1, ]
  res              <- pop_mat[-1, -1] - survived
  res[nrow(res), ] <- NA
  res              <- rbind(matrix(NA, nrow = 1, ncol = ncol(res)), res)
  #res              <- cbind(matrix(NA, nrow = nrow(res), ncol = 1), res)
  res              <- migresid_net_surv_last_ageg(res, pop_mat, sr_mat)
  rownames(res)    <- rownames(pop_mat)
  colnames(res)    <- colnames(pop_mat)[-p]
  res
}

# Net migration for last age group is pop for that age group in
# year j, minus the people from the previous age group the survived
migresid_net_surv_last_ageg <- function(net_mig, pop_mat, sr_mat) {
  # TR: this uses position indexing.
  n <- nrow(pop_mat)
  p <- ncol(pop_mat)
  previous_year <- 1:(p - 1)
  survived <-
    (pop_mat[n, previous_year] + pop_mat[n - 1, previous_year]) *
    sr_mat[n, previous_year]

  net_mig[nrow(net_mig), ] <- pop_mat[n, 2:p] - survived
  net_mig
}

migresid_births <- function(pop_f_mat,
                            asfr_mat,
                            #fertility_index,
                            age_interval) {
  p         <- ncol(pop_f_mat)
  asfr_ages <- rownames(asfr_mat)
  # Sum female pop from previous year and this year
  # f_pop <- pop_f_mat[asfr_ages, -1] + pop_f_mat[asfr_ages, -p]
  yrs     <- colnames(pop_f_mat) %>% as.numeric()
  yrs_out <- yrs[-p] + diff(yrs) / 2
  f_expos <-  interp(
                pop_f_mat[asfr_ages, ], 
                datesIn = yrs, 
                datesOut = yrs_out, 
                method = "linear")
  asfr_years  <- yrs[-p] %>% as.character()
  # Births that occurred for all age groups for all years
  # based on the age-specific fertility rate (asfr) from
  # previous years to the population
  these_births <- age_interval * (f_expos * asfr_mat[ , asfr_years]) # / 1000
  these_births <- colSums(these_births)
  names(these_births) <- asfr_years
  # all_births <- c(NA, colSums(these_births))
  # col_names  <- attr(pop_f_mat, "dimnames")[[2]]
  # all_births <- stats::setNames(all_births, col_names)
  # all_births
  these_births
}

migresid_net_surv_first_ageg <- function(net_mig, pop_mat, births, sr_mat) {
  # 20 yrs of births
  # 21 yrs of population
  # 20 yrs of sr
  p    <- ncol(net_mig)
  pyrs <- colnames(pop_mat)[-1] 
  
  # TR: a little hack
  D    <- pyrs %>% as.numeric() %>% diff() %>% '['(1)
  byrs <- pyrs %>% as.numeric() %>% '-'(D ) %>% as.character()
  # TR: note net_mig col labels seem to be one too high
  # we want byrs indexing on the left
  net_mig[1, ] <- pop_mat[1, pyrs] - births[byrs] * sr_mat[1, byrs]
  net_mig
}


# Returns age/year matrices with upper/lower bounds
# for net migration based on the net migration and
# survival rates. These, I believe are the upper/lower
# bounds of a lexis surfave (which is why we do ^0.5).
migresid_bounds <- function(net_mig, sr_mat) {
  n <- nrow(net_mig)
  p <- ncol(net_mig)

  # Upper bound is net mig / 2 times the survival ratio ^ 0.5
  mig_upper      <- net_mig / (2 * sr_mat^0.5)
  mig_upper      <- cbind(matrix(NA, ncol = 1, nrow = n), mig_upper)
  mig_lower      <- mig_upper
  mig_upper[1, ] <- NA
  mig_upper[n, ] <- NA
  mig_lower[n, ] <- NA
  mig_lower      <- mig_lower[-1, ]
  empty_matrix   <- matrix(NA, ncol = ncol(mig_lower), nrow = 1)
  mig_lower      <- rbind(mig_lower, empty_matrix)

  # Estimate upper bounds for the first age group. Why
  # no lower bound for the first age group? because we have
  # no previous age group.
  p_upper        <- ncol(mig_upper)
  mig_upper[1, 2:p_upper] <- net_mig[1, -p_upper] / (sr_mat[1, -p_upper]^0.5)

  list(upper = mig_upper, lower = mig_lower)
}

# Updates last age group for all upper/lower bounds
migresid_bounds_last_ageg <- function(net_mig_m,
                                      net_mig_f,
                                      mig_upper_m,
                                      mig_lower_m,
                                      mig_upper_f,
                                      mig_lower_f) {


  # last age group
  n <- nrow(mig_upper_m)
  p <- ncol(mig_upper_m)

  mig_lower_m[n - 1, ] <- mig_upper_m[n - 1, ]
  mig_lower_f[n - 1, ] <- mig_upper_f[n - 1, ]
  mig_upper_m[n, 2:p] <- net_mig_m[n, -p] * 0.5
  mig_upper_f[n, 2:p] <- net_mig_f[n, -p] * 0.5
  mig_lower_m[n, 2:p] <- net_mig_m[n, -p] * 0.5
  mig_lower_f[n, 2:p] <- net_mig_f[n, -p] * 0.5

  list(
    mig_lower_m = mig_lower_m,
    mig_upper_m = mig_upper_m,
    mig_lower_f = mig_lower_f,
    mig_upper_f = mig_upper_f
  )
}


mig_resid_dim_checker <- function(arg_list){

  # TR: objectives, either we get args from a properly captured arg_list,
  # or we simply pass in all args by name (maybe the easiest to be certain of)
  # ground rules:
  # age ranges should match for sr and pop. If they don't then we should trim to the
  # lowest common denominator, right?
  # year ranges depend on the input:
  # sr, asfr, srb need to have same years, but pop needs one extra year on the right side.
  
  # Each data argument should be given adequate dimnames for purposes of named selection
  # Each data argument should be trimmed as appropriate for conformable computations
  # If trimming happens, we warn if verbose.
  # This function basically just needs to return data inputs whose dimensions are
  # guaranteed to not cause problems in downstream mig_resid*() calcs.
  # the reason why we do this here is so that these many lines of code aren't repeated.
  
  pop_m_mat       <- arg_list$pop_m_mat
  pop_f_mat       <- arg_list$pop_f_mat
  sr_m_mat        <- arg_list$sr_m_mat
  sr_f_mat        <- arg_list$sr_f_mat
  asfr_mat        <- arg_list$asfr_mat
  srb_vec         <- arg_list$srb_vec
  ages            <- arg_list$ages
  ages_fertility  <- arg_list$ages_fertility
  
  # Make sure to add these year args to top level mig_resid* funcions.
  years_pop       <- arg_list$years_pop
  years_sr        <- arg_list$years_sr   
  years_asfr      <- arg_list$years_asfr
  years_srb       <- arg_list$years_srb
  
  # Make sure to add verbose arg to top level mig_resid*() functions.
  verbose         <- arg_list$verbose  
  # These are easier to insist on:
  stopifnot(all(dim(pop_m_mat) == dim(pop_f_mat)))
  stopifnot(all(dim(sr_m_mat) == dim(sr_f_mat)))
  
  # These args, could be NULL, so look to dimnames:
  if (is.null(ages)){
    ages             <- rownames(pop_m_mat) %>% as.numeric()
  }
  if (is.null(years_pop)){
    ages_fertility   <- rownames(asfr_mat) %>% as.numeric()
  }
  if (is.null(years_pop)){
    years_pop        <- colnames(pop_m_mat) %>% as.numeric()
  }
  if (is.null(years_asfr)){
    # TR: let's be careful that this doesn't end up hard coded at 15-45 or 15-49
    # when used throughout the functions. Hypothetically, it could have same ages
    # as pop or mort, but have 0s in non-fertile ages, make sense? This note
    # may be out of place, but came to mind here.
    years_asfr       <- colnames(asfr_mat) %>% as.numeric()
  }
  if (is.null(years_sr)){
    years_sr         <- colnames(sr_m_mat) %>% as.numeric()
  }
  if (is.null(years_srb)){
    years_srb        <- names(srb_vec) %>% as.numeric()
  }
 
  # Note, after the above, the years/ ages could still be NULL,
  # In this case we demand that dimensions already conform with expectations
  
  # For ages, we can guess from dims. For years, we can't guess from dims.
  # Therefore at least one of the year vectors needs to be non-NULL, AND
  # the dims of matrices to which NULL years correspond must already be correct.
  
  np    <- ncol(pop_f_mat)
  nsr   <- ncol(sr_m_mat)
  nfert <- ncol(asfr_mat)
  nsrb  <- length(srb_vec)
  
  dims_already_correct <- all(diff(c(np-1,nsr,nfert,nsrb) == 0))
  
  ind_nulls <- c(years_pop = is.null(years_pop), 
                 years_asfr = is.null(years_asfr), 
                 years_srb = is.null(years_srb), 
                 years_sr = is.null(years_sr))
  
  # it's easiest to just force users to give year ranges via args
  # or dimnames. If neither is available, just make them do it.
  if (any(ind_nulls)){
    stop("Year references must be given, either via function args or dimnames. Following references missing:\n",paste(names(ind_nulls)[ind_nulls],collapse=", "))
  }
  

  # 1) assign names
  colnames(pop_m_mat)    <- years_pop
  colnames(pop_f_mat)    <- years_pop
  colnames(asfr_mat)     <- years_fertility
  colnames(sr_m_mat)     <- years_sr
  colnames(sr_f_mat)     <- years_sr
  names(srb_vec)         <- years_srb
  
  # maybe there should be more thorough checks on age? 
  # we might be assigning NULL here...
  rownames(pop_m_mat)    <- ages
  rownames(pop_f_mat)    <- ages
  rownames(sr_m_mat)     <- ages
  rownames(sr_f_mat)     <- ages
  rownames(asfr_mat)     <- ages_fertility
  
  # 2) determine ranges
  # if dims aren't already correct
  yr1    <- max(c(min(years_pop),
                  min(years_sr),
                  min(years_asfr),
                  min(years_srb)))
  yrlast <- min(c(max(years_pop[-np]),
                  max(years_sr),
                  max(years_asfr),
                  max(years_srb)))
  
  interval      <- diff(years_asfr)[1] %>% as.integer()
  
  # just remember we need 1 more for pops!
  years_final   <- seq(yr1, yrlast, by = interval)
  years_final_p <- c(years_final, max(years_final) + interval)
  # trim
  pop_m_mat     <- pop_m_mat[, years_final_p ]
  
  if (ncol(pop_m_mat)!=np){
    # TR: should have one of these per trim step?
    if (verbose){
      cat("\npop_m_mat years have trimmed\n")
    }
  }
  
  pop_f_mat     <- pop_f_mat[, years_final_p ]
  sr_m_mat      <- sr_m_mat[, years_final ]
  sr_f_mat      <- sr_f_mat[, years_final ]
  asfr_mat      <- asfr_mat[, years_final ]
  srb_vec       <- srb_vec[ years_final ]
  
  out <- list(pop_m_mat = pop_m_mat,
              pop_f_mat = pop_f_mat,
              sr_m_mat = sr_m_mat,
              sr_f_mat = sr_f_mat,
              asfr_mat = asfr_mat,
              srb_vec = srb_vec)
  
  
}


# match.call.defaults <- function(...) {
#   call <- evalq(match.call(expand.dots = FALSE), parent.frame(1))
#   formals <- evalq(formals(), parent.frame(1))
#   
#   for(i in setdiff(names(formals), names(call)))
#     call[i] <- list( formals[[i]] )
#   
#   
#   match.call(sys.function(sys.parent()), call)
# }

# foo <- function(a,b=NULL,...){match.call.defaults()}
# foo(a=1,x=5)
# foo(a = 1, b = NULL, ... = pairlist(x = 5))
