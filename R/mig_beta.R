#' Estimate intercensal migration by comparing census population, by age and
#' sex, to the results of a RUP projection.
#'
#' @description
#' This methods projects population from the first starting point to next census
#' without migration and computes the "Net Census Error" (NCE) which is
#' Census - Estimate by age from projection. It then distributes the NCE over
#' the cohort parallelogram assuming uniform distribution assuming it is all
#' migration. It finalizes by summing the estimate by age groups across the entire
#' intercensal period to have a total migration during the entire period.
#'
#' @param c1 numeric vector. The first (left) census in single age groups
#' @param c2 numeric vector. The second (right) census in single age groups
#' @param date1 reference date of c1`. Either a Date class object or an unambiguous character string in the format "YYYY-MM-DD".
#' @param date2 reference date of c2`. Either a Date class object or an unambiguous character string in the format "YYYY-MM-DD".
#' @param age1 integer vector. single ages of `c1`
#' @param age2 integer vector. single ages of `c2`
#' @param dates_out vector of desired output dates coercible to numeric using `dec.date()`
#' @param lxMat numeric matrix containing lifetable survivorship, `l(x)`. Each row is an age group and each column a time point. At least two intercensal time points needed.
#' @param age_lx integer vector. Age classes in `lxMat`
#' @param dates_lx date, character, or numeric vector of the column time points for `lxMat`. If these are calendar-year estimates, then you can choose mid-year time points
#' @param births integer vector. Raw birth counts for the corresponding (sub)-population, one value per each year of the intercensal period including both census years. The first and last years should include all births in the given year; don't discount them in advance.
#' @param years_births numeric vector of calendar years of births.
#' @param location country name or LocID
#' @param sex character string, either `"male"`, `"female"`, or `"both"`
#' @param midyear logical. `FALSE` means all Jan 1 dates between `date1` and `date2` are returned. `TRUE` means all July 1 intercensal dates are returned.
#' @param verbose logical. Shall we send informative messages to the console?
#' @param ... optional arguments passed to \code{lt_single_qx}
#' @export
#'
#' @return a numeric vector of the total migration in the intercensal period
#' for each age. Ages are set as names of each migration estimate.
#'
#' @importFrom data.table := dcast
#'
#' @examples
#'
#' \dontrun{
#'   mig_beta(
#'   location = "Russian Federation",
#'   sex = "male",
#'   c1 = pop1m_rus2002,
#'   c2 = pop1m_rus2010,
#'   date1 = "2002-10-16",
#'   date2 = "2010-10-25",
#'   age1 = 0:100,
#'   births = c(719511L, 760934L, 772973L, 749554L, 760831L, 828772L, 880543L, 905380L, 919639L)
#' )
#' }
mig_beta <- function(
                     c1,
                     c2,
                     date1,
                     date2,
                     age1 = 1:length(c1) - 1,
                     age2 = 1:length(c2) - 1,
                     dates_out = NULL,
                     lxMat = NULL,
                     age_lx = NULL,
                     dates_lx = NULL,
                     births = NULL,
                     years_births = NULL,
                     location = NULL,
                     sex = "both",
                     midyear = FALSE,
                     verbose = TRUE,
                     ...
                     ) {

  # convert the dates into decimal numbers
  date1 <- dec.date(date1)
  date2 <- dec.date(date2)

  res_list <- rup(
    c1 = c1,
    c2 = c2,
    date1 = date1,
    date2 = date2,
    age1 = age1,
    age2 = age2,
    dates_out = dates_out,
    lxMat = lxMat,
    age_lx = age_lx,
    dates_lx = dates_lx,
    births = births,
    years_births = years_births,
    location = location,
    sex = sex,
    midyear = midyear,
    verbose = verbose,
    ... = ...
  )

  pop_jan1 <- res_list$pop_jan1
  dates_out <- res_list$dates_out

  age <- NULL
  year <- NULL
  cum_resid <- NULL
  decum_resid <- NULL
  discount <- NULL
  resid <- NULL

  # add "cumulative" residual to the RUP (pop_jan1_pre)
  pop_jan1[, `:=`(cum_resid = resid * discount)]
  pop_jan1 <- pop_jan1[!is.na(cohort)]
  # Group by cohort and decumulate the residual with the first
  # value being the first year of the cohort
  pop_jan1 <- pop_jan1[, decum_resid := c(cum_resid[1], diff(cum_resid)), key = cohort]

  # Transform the long data frame to wide with ages on rows, years on columns
  # and the decum_resid on the values.
  mat_resid <-
    data.table::dcast(
                  pop_jan1[, list(year, age, decum_resid)],
                  age ~ year,
                  value.var = "decum_resid"
                )
  # Sum over all ages to get a total decum_resid over all years for each age.
  mig <- stats::setNames(rowSums(mat_resid, na.rm = TRUE), mat_resid$age)
  mig
}





mig_beta_cwr <- function(mig, 
                         c1_females, 
                         c2_females, 
                         date1, 
                         date2, 
                         maternal_window = 30, 
                         maternal_min = 15){
  age <- names2age(mig)
  
  # conservative guess at how many child ages to cover:
  n_cohs <- as.integer(ceiling(date2) - floor(date1))
  
  mig_out <- mig
  for (i in 1:n_cohs){
    # index maternal ages
    a_min      <- i + maternal_min
    a_max      <- min(i + maternal_min + maternal_window, 49)
    mat_ind    <- a_min:a_max
    cwr_i      <- (c1_females[i] / sum(c1_females[mat_ind]) + c2_females[i] / sum(c2_females[mat_ind])) / 2
    # proportional to maternal neg mig.
    mig_out[i] <- cwr_i * sum(mig[mat_ind])
  }
  
  mig_out
}

# rough stb at constant child adjustment
mig_beta_constant_child <- function(mig, c1, c2, ageMax = 14){
  age <- names2age(mig)
  
  denom <- (c1 + c2) / 2
  
  ind <- age <= ageMax
  mig_rate_const <- sum(mig[ind]) / sum(denom[ind])
  
  mig[ind] <- denom[ind] * mig_rate_const
  
  mig
}


# One more function for optional old age smoothing




# TR: prep for constant child. Need denom for rates though,
# but ideally without re-calculating an intermediate object
# that was already needed. Maybe be we can get exposures
# from RUP. Hmm.
# mig_beta_constant_child <- function(mig, 
#                                     c1,
#                                     c2,
#                                     date1, 
#                                     date2, 
#                                     maternal_window = 30, 
#                                     maternal_min = 15){
#   age <- names2age(mig)
#   
#   # conservative guess at how many child ages to cover:
#   n_cohs  <- as.integer(ceiling(date2) - floor(date1))
#   
#   mig_out <- mig
#   
#   
#   
#   mig_out
# }
# 

