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

census_cohort_adjust <- function(Pop, Age, date){
  
  stopifnot(is_single(Age))
  
  date       <- dec.date(date)
  yr         <- floor(date)
  
  f1         <- date - yr
  
  upper_part_of_cohort <- Pop * f1
  lower_part_of_cohort <- Pop * (1 - f1)
  
  pop_out <- shift.vector(lower_part_of_cohort,1) + upper_part_of_cohort
  
  shift <- ceiling(f1)
  cohorts <- yr - Age - 1 + shift
  
  list(Pop = pop_out, Cohort = cohorts, Date = date, f1 = f1)
}
















