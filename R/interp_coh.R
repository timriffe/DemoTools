#' shift census populations to match single year cohorts
#' @description Matches the (single) ages of a census to single cohorts. For use in intercensal interpolations. Ages are potentially blended to match single cohort line assuming that the population in each age is uniformly distributed over the age group.
#' @param pop numeric vector. Population counts in age groups, presumably from a census with an exact reference date.
#' @param age integer vector. Lower bound of single age groups
#' @param date Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}.
#'@export
#' @examples
#' pop <- seq(10000,100,length.out = 101)
#' age <- 0:100
#' d1 <- "2020-01-01"
#' d2 <- "2020-07-01"
#' d3 <- "2020-12-21"
#'
#' census_cohort_adjust(pop, age, d1)
#' census_cohort_adjust(pop, age, d2)
#' census_cohort_adjust(pop, age, d3)
#' census_cohort_adjust(pop, age, 2020.5)
census_cohort_adjust <- function(pop, age, date){

  stopifnot(is_single(age))

  date       <- dec.date(date)
  yr         <- floor(date)

  f1         <- date - yr

  upper_part_of_cohort <- pop * f1
  lower_part_of_cohort <- pop * (1 - f1)

  shift      <- ceiling(f1)
  pop_out    <- shift.vector(lower_part_of_cohort,shift) + upper_part_of_cohort

  cohorts    <- yr - age - 1 + shift

  list(pop = pop_out, cohort = cohorts, date = date, f1 = f1)
}

#' component-free intercensal cohort interpolation
#' @description Cohorts between two censuses are interpolated flexibly using linear, exponential, or power rules. The lower and upper intercensal triangles are filled using within-age interpolation. This function is experimental and still in development.
#' @seealso interp
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
#' @param country text string. The country of the census
#' @param sex character string, either `"male"`, `"female"`, or `"both"`
#' @param midyear logical. `FALSE` means all Jan 1 dates between `date1` and `date2` are returned. `TRUE` means all July 1 intercensal dates are returned.
#' @param verbose logical. Shall we send informative messages to the console?
#' @param ... optional arguments passed to
#' @export
#' @importFrom countrycode countrycode
#' @importFrom data.table := as.data.table melt data.table dcast between
#' @examples
#'
#' \dontrun{
#' interp_coh(
#' country = "Russian Federation",
#' sex = "male",
#' c1 = pop1m_rus2002,
#' c2 = pop1m_rus2010,
#' date1 = "2002-10-16",
#' date2 = "2010-10-25",
#' age1 = 0:100,
#' births = c(719511L, 760934L, 772973L, 749554L, 760831L, 828772L, 880543L, 905380L, 919639L)
#' )
#' }
interp_coh <- function(
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
                       country = NULL,
                       sex = "both",
                       midyear = FALSE,
                       verbose = TRUE,
                       ...
                       ) {

  stopifnot(length(age1) == length(c1))
  stopifnot(length(age2) == length(c2))

  stopifnot(is_single(age1))
  stopifnot(is_single(age2))

  if (length(age1) != length(age2) & verbose){
    cat("\nFYI: age ranges are different for c1 and c2\nWe'll still get intercensal estimates,\nbut returned data will be chopped off after age", max(age1), "\n")
  }

  # If lxMat or births are missing -- message requiring country and sex
  if (is.null(lxMat) & is.null(country)) {
    stop("lxMat not specified, please specify country and sex\n")
  }
  if (is.null(births) & is.null(country)) {
    stop("births not specified, please specify country and sex\n")
  }

  # convert the dates into decimal numbers
  date1 <- dec.date(date1)
  date2 <- dec.date(date2)

  if (is.na(date1) | is.na(date2)){
    stop("\nCensus dates didn't parse\n")
  }

  if (!is.null(lxMat) && ncol(lxMat) == 1) {
    stop("lxMat should have at least two or more dates as columns. lxMat contains only one column") #nolintr
  }

  # TR: resolve dates_out
  # if some dates were given, let's coerce to numeric and ensure valid
  if (!is.null(dates_out)){
    dates_out          <- sapply(dates_out, dec.date)
    if (any(is.na(dates_out))){
      cat("\nSome dates_out didn't parse, FYI, you should have a look\n")
      dates_out          <- dates_out[!is.na(dates_out)]
    }
    if (length(dates_out) == 0){
      stop("\nno valid dates to interpolate to\n")
    }


    # if we still have valid dates, then check we're not extrapolating
    dates_out_keep     <- data.table::between(dates_out,
                                              date1,
                                              date2,
                                              incbounds = FALSE)
    dates_out_for_real <- dates_out[dates_out_keep]

    # warn about any dates lost due to extrap request:
    if (length(dates_out_for_real) != length(dates_out) & verbose){
      cat("\nFollowing dates requested, but not returned\nbecause they'd require extrapolation:\n",paste(dates_out[!dates_out_keep],collapse = ", "),"\n")
    }
    if (length(dates_out) == 0){
      stop("\nuh oh! This method is strictly for cohort component interpolation\nYour requested dates_out didn't have anything between date1 and date2\n")
    }
  }

  if (any(c1 < 0)) stop("No negative values allowed in `c1`")
  if (any(c2 < 0)) stop("No negative values allowed in `c2`")
  if (any(lxMat < 0)) stop("No negative values allowed in `lxMat`")

  # If dates_out not given, then we resolve using the midyear argument.
  # If FALSE (default) we return intermediate Jan 1, not including c1 and c2
  # If TRUE we return intermediate July 1 (.5) dates, not including c1 and c2
  if (is.null(dates_out)){
    if (! midyear){
      # jan 1 dates
      left_date  <- floor(date1) + 1
      right_date <- ceiling(date2) - 1
      dates_out  <- left_date:right_date
    }
    if (midyear){
      left_date  <- floor(date1) + .5
      right_date <- ceiling(date2) - .5
      dates_out  <- left_date:right_date
      dates_out_lgl  <- data.table::between(dates_out,
                                            date1,
                                            date2,
                                            incbounds = FALSE)
      dates_out <- dates_out[dates_out_lgl]
    }
  }


  DD <- date2 - date1
  if (DD >= 15 & verbose){
    cat("\nFYI, there are",DD,"years between c1 and c2\nBe wary.\n")
  }

  # let's store the proportions separately
  f1    <- date1 %>% magrittr::subtract(date1 %>% floor)
  f2    <- date2 %>% magrittr::subtract(date2 %>% floor)

  # get the lexis surface of survival probabilities
  if (is.null(lxMat)){
    if (verbose) cat(paste0("\nlxMat not provided. Downloading lxMat for ", country, ", gender: ", "`", sex, "`, for years between ", round(date1, 1), " and ", round(date2, 1), "\n"))

    pxt <- suppressMessages(
      interp_coh_download_mortality(country, sex, date1, date2, OAnew = max(age1) + 1)
    )
  } else {

    if (is.null(dates_lx)){
      # if lx dates not given we assume dates evenly distributed from date1 to date2?
      dates_lx <- seq(date1,date2,length.out = ncol(lxMat))
      if (verbose) {
        cat("lxMat specified, but not dates_lx\nAssuming:",paste(dates_lx,collapse=", "),"\n")
      }
    }

    available_dates <- data.table::between(dates_lx, date1, date2)
    if (!all(available_dates)) stop("All `dates_lx` must be within the range of `date1` and `date2`")

    # if the shortest distance from dates_lx to date1 or date2 is greater than 7
    # warn
    dates_df <- expand.grid(dates_lx = dates_lx, dates = c(date1, date2))
    dates_df$diff <- with(dates_df, abs(dates_lx - dates))
    if (min(dates_df$diff) > 7 && verbose) {
      d_lx <- dates_df$dates_lx[which.min(dates_df$dif)]
      date_compare <- dates_df$dates[which.min(dates_df$dif)]
      cat(
        "The shortest distance from `dates_lx` (",
        d_lx,
        ") to `date1/date2`(",
        date_compare,
        ") is greater than 7 years. Be wary."
      )
    }

    ic_period <- date2 - date1
    lx_mm <- range(dates_lx)
    overlap <- min(c(lx_mm[2], date2)) - c(max(lx_mm[1], date1))
    extrap_low <- lx_mm[1] - min(lx_mm[1],date1)
    extrap_high <- max(lx_mm[2],date2) - lx_mm[2]
    t1 <- overlap / ic_period < .25
    t2 <- extrap_low > 6
    t3 <- extrap_high > 6
    if (any(c(t1, t2, t3))) cat("\nRange between `date1` and `date2` must overlap with `lx_dates` for at least 25% of the range or 6 years.\n")

    if (is.null(age_lx)){
      if (nrow(lxMat)  < 26){

        N      <- nrow(lxMat)
        age_lx <- c(0,1,seq(5,5*(N-2),by=5))
      } else {
        age_lx <- 1:nrow(lxMat) - 1
      }
      if (verbose) {
        cat("lxMat specified, but Age_lx missing\nAssuming:",paste(age_lx,collapse=", "),"\n")
      }
    }

    # ensure lx fills timepoints.
    # would like to pass ... here for the lifetable part
    pxt <- interp_coh_lxMat_pxt(
      lxMat = lxMat,
      dates_lx = dates_lx,
      age_lx = age_lx,
      date1 = date1,
      date2 = date2,
      OAnew = max(age1) + 1,
      ...)
  }

  yrs_births   <- seq(floor(date1), floor(date2), 1)

  # TR: if right-side is jan 1 then we can cut it off of pxt.
  if (f2 == 0){
    pxt        <- pxt[, -ncol(pxt)]
    yrs_births <- yrs_births[-length(yrs_births)]
    f2         <- 1
  }

  # fetch WPP births if not provided by user
  if (is.null(births)) {

    # load WPP births
    requireNamespace("DemoToolsData", quietly = TRUE)
    WPP2019_births <- DemoToolsData::WPP2019_births

    # format filtering criteria -- country and years
    cntr_iso3 <- countrycode::countrycode(
                                country,
                                origin = "country.name",
                                destination  = "iso3c")

    # filter out country and years
    ind       <- WPP2019_births$ISO3 == cntr_iso3 & WPP2019_births$Year %in% yrs_births
    b_filt    <- WPP2019_births[ind, ]
    bt        <- b_filt$TBirths
    SRB       <- b_filt$SRB

    # extract births depending on sex
    if (sex == "both")  births  <- bt
    if (sex == "male")   births  <- bt * SRB / ( 1 + SRB)
    if (sex == "female") births  <- bt / (SRB + 1)

    if (verbose){
      cat("Births fetched from WPP for:", paste(country, sex), "population, years", paste(yrs_births, collapse = ", "), "\n")
    }

  }

  # check length of births, also filter using provided dates if necessary
  if (!is.null(years_births)){
    stopifnot(length(births) == length(years_births))

    years_births <- floor(years_births)
    yrs_keep     <- data.table::between(years_births,
                                        min(yrs_births),
                                        max(yrs_births),
                                        incbounds = TRUE)

    births       <- births[yrs_keep]
  }

  # now that births should be available we can do this check.
  stopifnot(length(births) == length(yrs_births))

  # a note for future: interp_coh_download_mortality should use {countrycode} to
  # better match the country names. As of now, just Russia won't work
  # [ISSUE #166]

  # Since we're using data.table, we need to create these empty
  # variables to avoid having R CMD checks with no visible binding
  # for global variable.

  age      <- NULL
  year     <- NULL
  px       <- NULL
  lower    <- NULL
  upper    <- NULL
  triangle <- NULL
  adj      <- NULL
  value    <- NULL
  pop      <- NULL
  coh_p    <- NULL
  coh_lx   <- NULL
  pop_c2_obs <- NULL
  x        <- NULL
  y        <- NULL
  discount <- NULL
  .N       <- NULL
  .        <- NULL

  px_triangles <-
    pxt %>%
    data.table::as.data.table(keep.rownames = "age") %>%
    data.table::melt(
                  id.vars = "age",
                  variable.name = "year",
                  value.name = "px",
                  variable.factor = FALSE
                )

  # No need for assignment: data.table assigns without creating a copy
  px_triangles[, `:=`(age = as.numeric(age),
                      year = as.numeric(year),
                      lower = magrittr::raise_to_power(px, 0.5),
                      upper = magrittr::raise_to_power(px, 1 - 0.5))]

  px_triangles <-
    px_triangles[, list(age, year, lower, upper)] %>%
    data.table::melt(
                  id.vars = c("age", "year"),
                  measure.vars = c("lower", "upper"),
                  variable.name = "triangle",
                  value.name = "value",
                  variable.factor = FALSE
                )

  px_triangles[, `:=`(adj = ifelse(triangle == "upper", 1, 0))]
  px_triangles[, `:=`(cohort = magrittr::subtract(year, age) %>% magrittr::subtract(adj) %>% floor())]

  # cohort changes over the whole period
  # px_cum1 <- px_triangles[, list(n_triangles = .N, coh_p = prod(value)),
  #                          keyby = list(cohort)]

  # adjust the census population vectors
  c1c <- census_cohort_adjust(c1, age1, date1)
  c2c <- census_cohort_adjust(c2, age2, date2)

  # correction for the first year age 0 -- only take first for the remaining of
  # the year
  births[1] <- births[1] * (1 - f1)

  # correction for the last year age 0
  n_yrs         <- length(births)
  births[n_yrs] <- births[n_yrs] * f2

  cohort_dt <-
    data.table::data.table(
                  cohort = yrs_births,
                  pop = births
                )

  input <-
    data.table::data.table(cohort = c1c$cohort, pop = c1c$pop) %>%
    .[order(cohort)] %>%
    rbind(cohort_dt) %>%
    .[, list(pop = sum(pop)), keyby = list(cohort)]

  # population c2 observed
  pop_c2 <- data.frame(
    cohort = c2c$cohort,
    pop_c2_obs = c2c$pop
  )


  pop_jan1_pre <-
    px_triangles %>%
    .[, list(n_triangles = .N, coh_p = prod(value)), keyby = list(year, cohort)] %>%
    .[order(cohort, year)]

  pop_jan1_pre[, `:=`(coh_lx = cumprod(coh_p)), keyby = list(cohort)]

  pop_jan1_pre <- pop_jan1_pre[input, on = "cohort"]

  pop_jan1_pre[, `:=`(
    pop_jan1_pre = pop * coh_lx,
    age = floor(year) - cohort,
    year = floor(year) + 1
  )]

  pop_jan1_pre[, `:=`(year = ifelse(year == max(year), year + f2 - 1, year))]

  # calculate the discrepancy (migration) -- to be disrtibuted uniformly in
  # cohorts
  resid <-
    pop_jan1_pre %>%
    .[year == max(year)] %>%
    .[pop_c2, on = "cohort"]

  resid[, `:=`(resid = pop_c2_obs - pop_jan1_pre)]
  # Only used in the process for diagnostics
  # resid[, `:=`(rel_resid = resid / pop_c2_obs)]
  resid <- resid[, list(cohort, resid)]

  # This should just be one value per cohort.

  # determine uniform error discounts:

  resid_discounts <-
    stats::approx(
             x = c(date1, date2),
             y = c(0, 1),
             xout = yrs_births
           ) %>%
    data.table::as.data.table() %>%
    .[, list(year = x, discount = y)]

  # output
  pop_jan1 <-
    pop_jan1_pre %>%
    merge(resid, by = "cohort", all = TRUE) %>%
    merge(resid_discounts, by = "year", all = TRUE)

  pop_jan1[, `:=`(
    resid = ifelse(is.na(resid), 0, resid),
    discount = ifelse(year == max(year), 1, discount)
  )]

  pop_jan1[, `:=`(pop_jan1 = pop_jan1_pre + resid * discount)]
  pop_jan1 <- pop_jan1[!is.na(cohort)]

  PopAP <-
    pop_jan1 %>%
    .[, list(age, year, pop_jan1)] %>%
    data.table::dcast(age ~ year, value.var = "pop_jan1") %>%
    .[order(age)]

  matinterp <- PopAP[age <= max(age1), -1] %>% as.matrix()
  rownames(matinterp) <- age1

  # Handle negatives (small pops, or large negative residuals relative to pop size)
  ind <- matinterp < 0
  if (any(ind) & verbose){
    cat("\n",sum(ind),"negatives detected in output.\nThese have been imputed with 0s.\n")
    matinterp[ind] <- 0
  }


  yrsIn     <- as.numeric(colnames(matinterp))
  if (all(yrsIn > date1)){
    matinterp <- cbind(c1, matinterp)
    yrsIn     <- c(date1, yrsIn)
  }
  if (all(yrsIn < date2)){
    matinterp <- cbind(matinterp, c2[1:length(c2)])
    yrsIn     <- c(yrsIn, date2)
  }
  colnames(matinterp) <- yrsIn
  # now we either return Jan1 dates or July 1 dates.

  out <- interp(
    matinterp,
    datesIn = yrsIn,
    datesOut = dates_out,
    rule = 1
  )

  out
}

# old code kept for now ---------------------------------------------------

# c1 <- seq(10000,10,length.out = 10)
# c2 <- seq(15000,10,length.out = 10)
#
# d1 <- "2020-07-01"
# d2 <- "2025-10-14"
#
# c1_coh <-census_cohort_adjust(c1, 0:9, d1)
# c2_coh <-census_cohort_adjust(c2, 0:9, d2)
#
# cohs_match <-
#
# matrix(c1_coh$pop)
#
# interp()
#
## commenting out interp_coh_bare won't be used
# interp_coh_bare <- function(c1, c2, date1, date2, age1, age2, ...){
#
#   date1 <- dec.date(date1)
#   date2 <- dec.date(date2)
#
#   # !!! do we plan to allow age1 != age2 ?
#
#   c1c <-census_cohort_adjust(c1, age1, date1)
#   c2c <-census_cohort_adjust(c2, age2, date2)
#
#   # Connect cohorts observed (completely) in both censuses
#   obs_coh <- intersect(c1c$cohort, c2c$cohort)
#
#   # remove first cohort if not observed in full
#   if(c1c$date - c1c$cohort[1] != 1){
#     obs_coh <- obs_coh[-1]
#   }
#
#   # Tim: select, make some intermediate data objects as necessary
#
#   # fully observed cohorts in a pop matrix
#   obs_coh_mat <- cbind(
#     c1c$pop[match(obs_coh, c1c$cohort)],
#     c2c$pop[match(obs_coh, c2c$cohort)]
#   )
#   # set names
#   dimnames(obs_coh_mat) <- list(obs_coh, c(c1c$date, c2c$date))
#
#   # Tim: then use interp()
#
#   # interpolate
#   dates_in <- dimnames(obs_coh_mat)[[2]] %>% as.numeric()
#   dates_out <- seq(floor(dates_in[1]), ceiling(dates_in[-1]), 1)
#   # what should be the default behavior here?
#   # I start off with the one year step period, inclusive
#
#   interpolated_coh_mat <- interp(
#     popmat = obs_coh_mat,
#     datesIn = dates_in,
#     datesOut = dates_out,
#     method = "linear",
#     rule = 2
#   )
#
#
#
#   ########
#
#   # Do something to fill in the lower triangle
#   # The most basic thing (suggested by Patrick)
#   # Just do a between-age interpolation and select out
#   # that triangle.
#
#   # !!! between-age interpolation
#   period_mat <-  cbind(c1, c2)
#   # set names
#   dimnames(period_mat) <- list(age1, c(date1, date2))
#
#
#   interpolated_period_mat <- interp(
#     popmat = period_mat,
#     datesIn = dimnames(period_mat)[[2]] %>% as.numeric(),
#     datesOut =  seq(floor(dates_in[1]), ceiling(dates_in[-1]), 1) ,
#     method = "linear",
#     rule = 2
#   )
#   ########
#
#   # Now fill in the upper triangle, doing something
#   # simple and robust
#
#   ########
#
#   # now the task is to take interpolated_period_mat
#   # and overwrite the matching values from interpolated_coh_mat
#
#   # achieved using a for loop that iterates across columns
#   # so I use the period interpolated matrix as canvas
#   # and overwrite the matching values from the cohort matrix
#
#   # dup interpolated_period_mat
#   out <- interpolated_period_mat
#
#   for (i in dimnames(interpolated_coh_mat)[[2]]) {
#     # take the i-th column from cohort interpolated matrix
#     replacement <- interpolated_coh_mat[,i]
#     # calculate the corresponding ages fo the interpolated values
#     ages <- as.numeric(i) - as.numeric(names(replacement))
#     # overwrite the cohort values in the period matrix
#     out[ages,i] <- replacement
#   }
#
#
#   # The remaining task is to frame the output
#   return(out)
#
# }
#
#
# # above rudimentary code works; below goes the new development
#
# # c1 = seq(10000,10,length.out = 10); c2 = seq(15000,10,length.out = 10); date1 = "2020-07-01"; date2 = "2025-10-14"; age1 = 0:9; age2 = 0:9
#
#
#
# canvas <-  interpolated_period_mat %>%
#   as_tibble(rownames = "age") %>%
#   pivot_longer(names_to = "year", values_to = "value", cols = -age) %>%
#   mutate(
#     age = age %>% as.numeric,
#     year = year %>% as.numeric,
#     cohort = year - age
#   )
#
# patch <- interpolated_coh_mat %>%
#   as_tibble(rownames = "cohort") %>%
#   pivot_longer(
#     names_to = "year", values_to = "value", cols = -1,
#     values_drop_na = TRUE
#   )%>%
#   mutate(
#     cohort = cohort %>% as.numeric,
#     year = year %>% as.numeric,
#     age = year - (cohort+1)
#   )
#
# final <- canvas %>%
#   rows_upsert(
#     patch %>% filter(!year %in% c(2020, 2026)),
#     by = c("age", "year")
#   ) %>%
#   select(-cohort) %>%
#   pivot_wider(names_from = year)
#
#
#
#
# view_ap <- function(long_apc_df) {
#   long_apc_df %>%
#     select(-cohort) %>%
#     pivot_wider(names_from = year)
# }
#
# patch %>% view_ap
#
# canvas %>% view_ap
#
# final %>% view_ap
#
#
#
#
#
#
# # the survival probabilities approach -------------------------------------
# load_this <- FALSE
# if (load_this) {
#   # blocking this off lets us to
#   devtools::load_all()
#   library(magrittr)
#   library(tidyverse)
# pxt <- suppressMessages(interp_coh_download_mortality("Russian Federation","male","2002-10-16","2010-10-25"))
# # a note for future: interp_coh_download_mortality should use {countrycode} to better match the country names. As of now, just Russia won't work
#
# # convert the AP output to CP
# px_triangles <- pxt %>%
#   as_tibble(rownames = "age") %>%
#   pivot_longer(
#     names_to = "year", values_to = "px", cols = -1,
#     values_drop_na = TRUE
#   ) %>%
#   mutate(
#     age = age %>% as.numeric,
#     year = year %>% as.numeric,
#     # cohort = floor(year) - age
#   ) %>%
#   # in triangles
#   mutate(
#     # year_frac = year - floor(year) # for now just .5 ~ sqrt
#     lower = px %>% raise_to_power(.5),
#     upper = px %>% raise_to_power(1 - .5) # .5 to be changed to year_frac
#   ) %>%
#   select(-px) %>%
#   pivot_longer(
#     names_to = "triangle", values_to = "value", cols = lower:upper
#   ) %>%
#   mutate(
#     adj = case_when(triangle=="upper" ~ 1, TRUE ~ 0),
#     cohort = year %>% subtract(age) %>% subtract(adj) %>% floor
#   )
#
#
#
# # cohort changes over the whole period
# px_cum <- px_triangles %>%
#   group_by(cohort) %>%
#   summarise(
#     n_triangles = n(),
#     coh_p = value %>% prod
#   ) %>%
#   ungroup()
#
# # foo %>% interp_coh_tidy_pc("1971-01-14","1978-02-01") %>% view
#
# # # generate two census populations -- single years of age
# # set.seed(911)
# # c1 <- spline(c(6,7,9,8,7,6,4,2,1)*1e3,n = 101)$y * runif(101, 1, 1.1)
# # set.seed(444)
# # c2 <- spline(c(6,7,9,8,7,6,4,2,1)*1e3,n = 101)$y * runif(101, 1.05, 1.15)
# # # births as random +-10% of the c1 and c2 age 0 average
# # births <- runif(6, .9*mean(c1[1], c2[1]), 1.1*mean(c1[1], c2[1])) %>% round
#
# # EXAMPLE DATA: Russian male population from the last two censuses
# # 2002 -- http://www.demoscope.ru/weekly/ssp/rus2002_01.php
# # 2020 -- http://www.demoscope.ru/weekly/ssp/rus_age1_10.php
# rus2002m <- c(682698L, 641551L, 644671L, 644652L, 662998L, 659306L, 678341L, 717053L, 740366L, 753300L, 875113L, 963123L, 1081671L, 1145059L, 1247787L, 1314341L, 1291147L, 1266227L, 1306873L, 1325599L, 1234028L, 1162951L, 1170248L, 1115312L, 1100598L, 1088833L, 1092321L, 1070733L, 1045802L, 1016461L, 1061391L, 994896L, 1007712L, 933628L, 916902L, 929632L, 957895L, 981477L, 1039571L, 1116279L, 1195521L, 1210704L, 1278766L, 1216728L, 1182385L, 1167289L, 1123058L, 1117150L, 1087663L, 998307L, 1035886L, 951627L, 960428L, 963751L, 730354L, 798841L, 604983L, 382611L, 298788L, 280702L, 493677L, 625270L, 694930L, 741777L, 695339L, 693911L, 559111L, 467811L, 358252L, 364999L, 427681L, 405822L, 435844L, 385155L, 379150L, 317841L, 258185L, 193023L, 154406L, 112987L, 89944L, 73858L, 63570L, 54955L, 47194L, 30300L, 28748L, 29419L, 26635L, 20166L, 16673L, 10857L, 8189L, 4839L, 3333L, 2287L, 1458L, 984L, 644L, 488L, 967L)
# rus2010m <- c(842354L, 859562L, 849138L, 788376L, 744105L, 750282L, 748514L, 746626L, 709493L, 675127L, 683827L, 656887L, 678395L, 669374L, 696685L, 743449L, 774172L, 800765L, 923952L, 1035555L, 1167860L, 1187193L, 1252421L, 1300116L, 1262584L, 1247974L, 1230926L, 1249086L, 1156502L, 1125283L, 1182017L, 1088248L, 1073221L, 1038733L, 1051852L, 1046293L, 1008882L, 983045L, 985075L, 949072L, 980924L, 881915L, 866214L, 859808L, 885432L, 926771L, 951739L, 1015812L, 1051749L, 1093184L, 1155128L, 1076307L, 1043777L, 1005283L, 967830L, 964217L, 919814L, 837341L, 841362L, 789019L, 787516L, 775999L, 585545L, 624976L, 471186L, 295668L, 222526L, 205594L, 336318L, 431670L, 471562L, 485883L, 446533L, 438107L, 337694L, 273086L, 198303L, 190828L, 210878L, 195219L, 200564L, 162820L, 151191L, 120794L, 93394L, 66247L, 48072L, 32932L, 23840L, 18087L, 13839L, 10228L, 7790L, 4327L, 3544L, 3137L, 2380L, 1666L, 1137L, 687L, 1379L)
# # MALE BIRTHS IN RUSSIA 2002--2010 (https://www.fedstat.ru/indicator/31606)
# births <- c(
#   719511L, 760934L, 772973L, 749554L, 760831L,
#   828772L, 880543L, 905380L, 919639L
# )
#
#
# c1 = rus2002m; c2 = rus2010m
#
# date1 = "2002-10-16"; date2 = "2010-10-25"; age1 = 0:100; age2 = 0:100
#
# date1 <- dec.date(date1)
# date2 <- dec.date(date2)
#
# # let's store the proportions separately
# f1 <- date1 %>% subtract(date1 %>% floor)
# f2 <- date2 %>% subtract(date2 %>% floor)
#
# # IK: do we plan to allow age1 != age2 ?
# # TR: for now we force them to be equal. Later a wrapper can take care of cleaning up these details.
# # we have OPAG() to extend open ages; graduate() to spit to single-
# # any other adjustments should be done in advance (smoothing, __ )
#
# c1c <-census_cohort_adjust(c1, age1, date1)
# c2c <-census_cohort_adjust(c2, age2, date2)
#
# # correction for the first year age 0 -- only take first for the remaining of the year
# births[1] <- births[1] * (1 - f1) # TR: good
#
# # TR: correction for the last year age 0
# n_yrs <- length(births)
# births[n_yrs] <- births[n_yrs] * f2
#
# # input
# input <- tibble(
#   cohort = c1c$cohort,
#   pop = c1c$pop
# ) %>%
#   arrange(cohort) %>%
#   bind_rows(
#     tibble(
#       cohort = 1:length(births) + floor(date1) - 1,
#       pop = births
#     )
#   ) %>%
#   # treat the duplicated cohort of the first census year, 2002
#   group_by(cohort) %>%
#   summarise(
#     pop = pop %>% sum,
#     .groups = "drop"
#   )
#
# # population c2 observed
# pop_c2 <- tibble(
#   cohort = c2c$cohort,
#   pop_c2_obs = c2c$pop
# )
#
# # # cohort survival to the second census
# # input %>%
# #   left_join(px_cum, by = "cohort") %>%
# #   mutate(pop_c2_prj = pop * coh_p) %>%
# #   left_join(pop_c2, by = "cohort") %>%
# #   mutate(
# #     discrepancy = pop_c2_obs - pop_c2_prj,
# #     disc_rel = discrepancy / pop_c2_obs * 100
# #   )
#
#
# # estimates of jan 1 population,
# # prior to redistribution of the residual
# # includes partial year estimate on the right-hand side,
# # excludes c1.
#
# pop_jan1_pre <-
#   px_triangles %>%
#   group_by(year, cohort) %>%
#   summarise(
#     n_triangles = n(),
#     coh_p = value %>% prod,
#     .groups = "drop"
#   ) %>%
#   arrange(cohort, year) %>%
#   group_by(cohort) %>%
#   mutate(coh_lx = cumprod(coh_p)) %>%
#   ungroup() %>%
#   left_join(input, by = "cohort") %>%
#   mutate(
#     pop_jan1_pre = pop * coh_lx,
#     age = floor(year) - cohort,
#     year = floor(year) + 1,
#     year = ifelse(year == max(year), year + f2 - 1, year)
#   )
#
# resid <-
#   pop_jan1_pre %>%
#   dplyr::filter(year == max(year)) %>%
#   left_join(pop_c2, by = "cohort") %>%
#   mutate(
#     resid = pop_c2_obs - pop_jan1_pre,
#     rel_resid = resid / pop_c2_obs
#   ) %>%
#   select(cohort, resid)
#
# # determine uniform error discounts:
#
# resid_discounts <-
#   approx(
#     x=c(date1, date2),
#     y=c(0,1),
#     xout=seq(ceiling(date1),floor(date2))
#   ) %>%
#   as.data.frame() %>%
#   select(year = x, discount= y)
#
# pop_jan1 <-
#   pop_jan1_pre %>%
#   left_join(resid, by = "cohort") %>%
#   left_join(resid_discounts, by = "year") %>%
#   mutate(
#     resid = ifelse(is.na(resid),0,resid),
#     discount = ifelse(year == max(year),1,discount),
#     pop_jan1 = pop_jan1_pre + resid * discount
#   )
#
# pop_jan1 %>%
#   # reshape2::acast(age~year, value.var = "pop_jan1") %>%
#   select(age, year, pop_jan1) %>%
#   pivot_wider(names_from = year, values_from = "pop_jan1") %>%
#   view()
#
# }
#
#

# This script does nothing yet, still in development, deciding how to
# graduate abridged lifetables

# This is a temporary script to hold a utility function for
# interp_coh()

# goal will be to fill a mortality surface between two censuses.
# args should be date1, date2, country.

# A few temporary functions internal to interp_coh(). These can be replaced as better
# or more efficient options become available.

lt_a2s_chunk <- function(chunk, OAnew, ...){
  ndx <- chunk$dx
  nLx <- chunk$Lx
  Age <- chunk$x

  lt_abridged2single(ndx = ndx,
                     nLx = nLx,
                     Age = Age,
                     OAnew = OAnew,
                     ...)
}

interp_coh_download_mortality <- function(country, sex, date1, date2, OAnew = 100){

  . <- NULL

  date1      <- dec.date(date1)
  date2      <- dec.date(date2)

  year1      <- floor(date1) + 1
  year2      <- floor(date2)

  year_seq   <- year1:year2

  dates_out  <- c(dec.date(date1), year_seq)


  PX <- suppressMessages(lapply(dates_out,fertestr::FetchLifeTableWpp2019,
                                locations = country,
                                sex = sex)) %>%
    lapply(function(X){
      X[,c("year","x","dx","Lx")]
    }) %>%
    lapply(lt_a2s_chunk, OAnew = OAnew) %>%
    lapply(function(X){
      1 - X$nqx
    }) %>%
    do.call("cbind",.)


  dimnames(PX)   <- list(0:OAnew, dates_out)

  PX[PX > 1]     <- 1
  # discount first and last periods.

  f1             <- diff(dates_out)[1]
  f2             <- date2 - floor(date2)

  # assume linear px change within age class
  PX[, 1]        <- PX[, 1] ^f1
  PX[,ncol(PX)]  <- PX[, ncol(PX)] ^f2

  PX
}

# lxMat <-suppressMessages(lapply(dates_out,fertestr::FetchLifeTableWpp2019,
#                              locations = country,
#                              sex = sex) %>%
#                         lapply("[[","lx") %>%
#                         dplyr::bind_cols() %>%
#                         as.matrix())

interp_coh_lxMat_pxt <- function(lxMat,
                                 dates_lx,
                                 age_lx,
                                 date1,
                                 date2,
                                 OAnew, ...){
  # TR: this is a temp functin, a stop-gap. Some redundant code with
  # interp_coh_download_mortality(), which it itself temporary.
  # the age graduation will move to lt_abridged2single() as soon as it's
  # fixed.
  date1      <- dec.date(date1)
  date2      <- dec.date(date2)

  year1      <- floor(date1) + 1
  year2      <- floor(date2)

  year_seq   <- year1:year2

  dates_out  <- c(dec.date(date1), year_seq)

  # get ndx andnLx from lt_abridged()

  a1  <- 0:OAnew
  qx1 <- matrix(ncol = ncol(lxMat),
                nrow = length(a1),
                dimnames = list(a1,
                                dates_lx))
  for (i in 1:ncol(lxMat)){

    if (is_abridged(age_lx)){
      LTA     <- lt_abridged(Age = age_lx,
                             lx = lxMat[, i],
                             OAnew = OAnew,
                             radix = 1e6,
                             ...)
      LT1     <- lt_abridged2single(ndx = LTA$ndx,
                                    nLx = LTA$nLx,
                                    Age = LTA$Age,
                                    OAnew = OAnew,
                                    ...)
      qx1[, i] <- LT1$nqx
    } else {
      qx             <- lt_id_l_q(lxMat[, i])

      LT1 <- lt_single_qx(nqx = qx,
                          Age=1:length(qx)-1,
                          OAnew = OAnew,
                          ...)


      qx1[, i] <- LT1$nqx
    }

  }

  # We do linear interpolation of the logit-transformed qx.
  logit_qx  <- log(qx1 / (1 - qx1))

  logit_qx_interp     <-
    interp(
      popmat = logit_qx,
      datesIn = dates_lx,
      datesOut = dates_out,
      rule = 2)
  # transform back
  QX            <- exp(logit_qx_interp) / (1 + exp(logit_qx_interp))

  QX[nrow(QX), ]  <- 1


  f1            <- diff(dates_out)[1]
  f2            <- date2 - floor(date2)

  # assume linear px change within age class
  PX            <- 1 - QX
  PX[,1]        <- PX[, 1] ^f1
  PX[,ncol(PX)] <- PX[, ncol(PX)] ^f2


  PX
}
