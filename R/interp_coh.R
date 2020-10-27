#' shift census populations to match single year cohorts
#' @description Matches the (single) ages of a census to single cohorts. For use in intercensal interpolations. Ages are potentially blended to match single cohort line assuming that the population in each age is uniformly distributed over the age group.
#' @param Age integer. Lower bound of single age groups
#' @param date Either a \code{Date} class object or an unambiguous character string in the format \code{"YYYY-MM-DD"}.

#' @examples
#' Pop <- seq(10000,100,length.out = 101)
#' Age <- 0:100
#' d1 <- "2020-01-01"
#' d2 <- "2020-07-01"
#' d3 <- "2020-12-21"
#' 
#' census_cohort_adjust(Pop,Age,d1)
#' census_cohort_adjust(Pop,Age,d2)
#' census_cohort_adjust(Pop,Age,d3)
#' census_cohort_adjust(Pop,Age,2020.5)
census_cohort_adjust <- function(Pop, Age, date){
  
  stopifnot(is_single(Age))
  
  date       <- dec.date(date)
  yr         <- floor(date)
  
  f1         <- date - yr
  
  upper_part_of_cohort <- Pop * f1
  lower_part_of_cohort <- Pop * (1 - f1)
  
  shift      <- ceiling(f1)
  pop_out    <- shift.vector(lower_part_of_cohort,shift) + upper_part_of_cohort

  cohorts    <- yr - Age - 1 + shift
  
  list(Pop = pop_out, Cohort = cohorts, Date = date, f1 = f1)
}

# C1 <- seq(10000,10,length.out = 10)
# C2 <- seq(15000,10,length.out = 10)
# 
# d1 <- "2020-07-01"
# d2 <- "2025-10-14"
# 
# C1_coh <-census_cohort_adjust(C1, 0:9, d1)
# C2_coh <-census_cohort_adjust(C2, 0:9, d2)
# 
# cohs_match <- 
# 
# matrix(C1_coh$Pop)
# 
# interp()

interp_coh_bare <- function(C1, C2, date1, date2, Age1, Age2, ...){
  
  C1_coh <-census_cohort_adjust(C1, Age1, date1)
  C2_coh <-census_cohort_adjust(C2, Age2, date2)
  
  # Connect cohorts observed (completely) in both censuses
  
  # select, make some intermediate data objects as necessary
  # then use interp()
  
  
  ########
  
  # Do something to fill in the lower triangle
  # The most basic thing (suggested by Patrick)
  # Just do a between-age interpolation and select out
  # that triangle.
  
  
  ########
  
  # Now fill in the upper triangle, doing something
  # simple and robust
  
  ########
  
  
  
}





