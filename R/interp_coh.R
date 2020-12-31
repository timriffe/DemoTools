#' shift census populations to match single year cohorts
#' @description Matches the (single) ages of a census to single cohorts. For use in intercensal interpolations. Ages are potentially blended to match single cohort line assuming that the population in each age is uniformly distributed over the age group.
#' @param age integer. Lower bound of single age groups
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

# c1 = seq(10000,10,length.out = 10); c2 = seq(15000,10,length.out = 10); date1 = "2020-07-01"; date2 = "2025-10-14"; age1 = 0:9; age2 = 0:9

#' component-free intercensalcohort interpolation
#' @description Cohorts between two censuses are interpolated flexibly using linear, exponential, or power rules. The lower and upper intercensal triangles are filled using within-age interpolation. This function is experimental and still in development.
#' @seealso interp
#' @param c1 numeric vector. The first (left) census in single age groups
#' @param c1 numeric vector. The second (right) census in single age groups
#' @param date1 reference date of c1`. Either a Date class object or an unambiguous character string in the format "YYYY-MM-DD".
#' @param date2 reference date of c2`. Either a Date class object or an unambiguous character string in the format "YYYY-MM-DD".
#' @param age1 integer vector. single ages of `c1`
#' @param age2 integer vector. single ages of `c2` 
#' @export
interp_coh_bare <- function(c1, c2, date1, date2, age1, age2, ...){
  
  date1 <- dec.date(date1)
  date2 <- dec.date(date2)
  
  # !!! do we plan to allow age1 != age2 ?
  
  c1c <-census_cohort_adjust(c1, age1, date1)
  c2c <-census_cohort_adjust(c2, age2, date2)
  
  # Connect cohorts observed (completely) in both censuses
  obs_coh <- intersect(c1c$cohort, c2c$cohort)
  
  # remove first cohort if not observed in full
  if(c1c$date - c1c$cohort[1] != 1){
    obs_coh <- obs_coh[-1]
  }
  
  # Tim: select, make some intermediate data objects as necessary
  
  # fully observed cohorts in a pop matrix
  obs_coh_mat <- cbind(
    c1c$pop[match(obs_coh, c1c$cohort)],
    c2c$pop[match(obs_coh, c2c$cohort)]
  ) 
  # set names 
  dimnames(obs_coh_mat) <- list(obs_coh, c(c1c$date, c2c$date))
  
  # Tim: then use interp()
  
  # interpolate
  dates_in <- dimnames(obs_coh_mat)[[2]] %>% as.numeric()
  dates_out <- seq(floor(dates_in[1]), ceiling(dates_in[-1]), 1) 
  # what should be the default behavior here? 
  # I start off with the one year step period, inclusive
  
  interpolated_coh_mat <- interp(
    popmat = obs_coh_mat, 
    datesIn = dates_in, 
    datesOut = dates_out, 
    method = "linear", 
    rule = 2
  )
  
  
  
  ########
  
  # Do something to fill in the lower triangle
  # The most basic thing (suggested by Patrick)
  # Just do a between-age interpolation and select out
  # that triangle.
  
  # !!! between-age interpolation
  period_mat <-  cbind(c1, c2) 
  # set names 
  dimnames(period_mat) <- list(age1, c(date1, date2))
  
  
  interpolated_period_mat <- interp(
    popmat = period_mat, 
    datesIn = dimnames(period_mat)[[2]] %>% as.numeric(), 
    datesOut =  seq(floor(dates_in[1]), ceiling(dates_in[-1]), 1) , 
    method = "linear", 
    rule = 2
  )
  ########
  
  # Now fill in the upper triangle, doing something
  # simple and robust
  
  ########
  
  # now the task is to take interpolated_period_mat  
  # and overwrite the matching values from interpolated_coh_mat
  
  # achieved using a for loop that iterates across columns
  # so I use the period interpolated matrix as canvas
  # and overwrite the matching values from the cohort matrix
  
  # dup interpolated_period_mat
  out <- interpolated_period_mat
  
  for (i in dimnames(interpolated_coh_mat)[[2]]) {
    # take the i-th column from cohort interpolated matrix
    replacement <- interpolated_coh_mat[,i]
    # calculate the corresponding ages fo the interpolated values
    ages <- as.numeric(i) - as.numeric(names(replacement))
    # overwrite the cohort values in the period matrix
    out[ages,i] <- replacement
  }
  
  
  # The remaining task is to frame the output
  return(out)
  
}
# 
# 
# canvas <-  interpolated_period_mat %>%
#   as_tibble(rownames = "age") %>%
#   pivot_longer(names_to = "year", values_to = "value", cols = -age) %>%
#   mutate(
#     age = age %>% as.numeric,
#     year = year %>% as.numeric,
#     cohort = year - age
#   ) %>%
#   view
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
#   ) %>%
#   view
# 
# final <- canvas %>%
#   # rows_upsert(
#   #   patch %>% filter(!year %in% c(2020, 2026)),
#   #   by = c("age", "year")
#   # ) %>%
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
# 
# foo <- interp_coh_download_mortality("France","male","1971-02-14","1978-07-01")
# 
# # try out
# #
# # boo <- interpolated_period_mat
# #
# # boo[is.numeric(boo)] <- 0
# #
# # foo <- interpolated_coh_mat
# #
# #
# #
# # for (i in dimnames(foo)[[2]]) {
# #
# #   replacement <- foo[,i]
# #   ages <- as.numeric(i) - as.numeric(names(replacement))
# #
# #   boo[ages,i] <- replacement
# # }
# 
# microbenchmark::microbenchmark(
#   base = {
#     for (i in dimnames(interpolated_coh_mat)[[2]]) {
#       # take the i-th column from cohort interpolated matrix
#       replacement <- interpolated_coh_mat[,i]
#       # calculate the corresponding ages fo the interpolated values
#       ages <- as.numeric(i) - as.numeric(names(replacement))
#       # overwrite the cohort values in the period matrix
#       interpolated_period_mat[ages,i] <- replacement
#     }
#   },
#   tidy = {
# 
#     final <- canvas %>%
#       rows_upsert(patch, by = c("age", "year")) %>%
#       select(-cohort) %>%
#       pivot_wider(names_from = year)
#   },
#   join = {
#     jouned <- canvas %>%
#       select(-cohort) %>%
#       left_join(patch %>% select(-cohort),
#                 by = c("age", "year")) %>%
#       transmute(age, year,
#                 value = case_when(!is.na(value.y)~value.y, TRUE~value.x)) %>%
#       pivot_wider(names_from = year)
#   }
# )
