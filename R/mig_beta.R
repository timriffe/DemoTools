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
#' Alternatively, a child adjustment and an old age adjustment can be applied.
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
#'
#' @param child_adjust The method with which to adjust the youngest age groups.
#' If \code{"none"}, no adjustment is applied (default). If
#' child-woman ratio (\code{"cwr"}) is chosen, the first cohorts reflecting the
#' difference between \code{date2 - date1} are adjusted (plus age 0). If
#' child constant ratio (\code{"constant"}) is chosen, the first 15 age groups
#' are adjusted.
#'
#' @param childage_max The maximum age from which to apply \code{child_adjust}.
#' By default, set to \code{NULL}, which gets translated into all the cohorts
#' between \code{date2} and \code{date1}. If \code{date2} is 2010 and
#' \code{date1} is 2002, the first 8 cohorts are adjusted. Otherwise, the user
#' can supply an integer.
#'
#' @param cwr_factor A numeric between 0 and 1 to which adjust the CWR method
#' for the young ages from \code{child_adjust}. \strong{This is only used
#' when \code{child_adjust} is \code{"cwr"}}.
#'
#' @param oldage_adjust The type of adjustment to apply to ages at and above
#' \code{oldage_min}. \code{'beers'} applies a beers graduation method
#' while \code{'mav'} applies a moving average with cascading on the tails.
#' For more information see \code{?mav} and \code{?graduation_beers}.
#'
#' @param oldage_min The minimum age from which to apply \code{oldage_adjust}.
#' By default, set to 65, so any adjustment from \code{oldage_adjust} will be
#' applied for 65+.
#'
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
#'
#' mig_beta(
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
                     child_adjust = c("none", "cwr", "constant"),
                     childage_max = NULL,
                     cwr_factor = 0.3,
                     oldage_adjust = c("none", "beers", "mav"),
                     oldage_min = 65,
                     ...) {
  child_adjust <- match.arg(child_adjust)
  oldage_adjust <- match.arg(oldage_adjust)

  # convert the dates into decimal numbers
  date1 <- dec.date(date1)
  date2 <- dec.date(date2)

  # If null, assume, the cohorts between censuses date2 and dates2
  if (is.null(childage_max)) {
    childage_max <- as.integer(ceiling(date2) - floor(date1))
  }

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

  # Child adjustment
  mig <-
    switch(
      child_adjust,
      "none" = mig,
      "cwr" = mig_beta_cwr(mig, c1, c2, date1, date2, n_cohs = childage_max, cwr_factor = cwr_factor),
      "constant" = mig_beta_constant_child(mig, c1, c2, ageMax = childage_max)
    )

  # Old age adjustment
  mig_oldage <-
    switch(
      oldage_adjust,
      "none" = mig,
      "beers" = graduate_beers(mig, as.integer(names(mig)), AgeInt = 1),
      "mav" = mav(mig, names(mig), tails = TRUE)
    )

  # Only apply the old age adjustment on ages above oldage_min
  ages_oldages <- as.integer(names(mig_oldage))
  mig[ages_oldages >= oldage_min] <- mig_oldage[ages_oldages >= oldage_min]

  mig
}


mig_beta_cwr <- function(mig,
                         c1_females,
                         c2_females,
                         date1,
                         date2,
                         maternal_window = 30,
                         maternal_min = 15,
                         n_cohs = NULL,
                         cwr_factor = 0.3) {

  age <- names2age(mig)

  # conservative guess at how many child ages to cover:
  if (is.null(n_cohs)) n_cohs <- as.integer(ceiling(date2) - floor(date1))

  mig_out <- mig
  for (i in 1:n_cohs) {
    # index maternal ages
    a_min <- i + maternal_min
    a_max <- min(i + maternal_min + maternal_window, 49)
    mat_ind <- a_min:a_max
    cwr_i <- (c1_females[i] / sum(c1_females[mat_ind]) + c2_females[i] / sum(c2_females[mat_ind])) / 2
    # proportional to maternal neg mig.
    mig_out[i] <- cwr_factor * cwr_i * sum(mig[mat_ind])
  }

  mig_out
}

# rough stb at constant child adjustment
mig_beta_constant_child <- function(mig, c1, c2, ageMax = 14) {
  age <- names2age(mig)

  denom <- (c1 + c2) / 2

  ind <- age <= ageMax
  mig_rate_const <- sum(mig[ind]) / sum(denom[ind])

  mig[ind] <- denom[ind] * mig_rate_const

  mig
}
