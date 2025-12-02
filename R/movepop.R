#' Project Population Between Two Dates Using Age-Specific Rates
#'
#' The `movepop()` function projects a population from an initial date to a desired date
#' using age-specific fertility (ASFR) and mortality rates. It supports single-year or 
#' five-year age group formats and optionally accounts for net migration. 
#' The function produces both projected population by age and sex and summary statistics.
#'
#' @param initial_date Numeric or Date. The starting year/date of the population.
#' @param desired_date Numeric or Date. The year/date to which the population is projected.
#' @param male_pop Numeric vector of male population counts by age.
#' @param female_pop Numeric vector of female population counts by age.
#' @param male_mx Numeric vector of male mortality rates by age.
#' @param female_mx Numeric vector of female mortality rates by age.
#' @param Age Numeric vector of ages Age that are present in `female_mx`, `male_mx`, `female_pop`, `male_pop`.

#' @param asfr Numeric vector of age-specific fertility rates corresponding to female ages.
#' @param annual_net_migrants Numeric scalar for annual net migration to include in growth (default 0).
#' @param age_out Character specifying returning age structure: `"single"`, `"five_year"`, or `"abridged"` (default `"auto"`).
#' @param first_asfr_age Numeric Optional first age at which the fertility starts e.g. 15
#' 
#' @details
#' - The function ensures that population and mortality vectors have equal lengths.
#' - Age groups are automatically generated if not provided, and ASFR values are aligned starting at the age group `"15"` (or equivalent).
#' - Projected populations are computed using exponential growth based on crude birth rate, crude death rate, and net migration rate.
#' - The function returns both age-specific projected populations and summary statistics.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{initial_summaries}}{A tibble summarising total births, deaths, crude rates, growth rate, total population, and adjustment factor.}
#'   \item{\code{projected_summaries}}{A tibble summarising projected total population by sex and overall sex ratio.}
#'   \item{\code{data}}{A tibble with age-specific populations, ASFR, births, deaths, projected populations by sex, both-sexes totals, and sex ratios.}
#' }
#'
#' @importFrom dplyr mutate summarise case_when
#' @importFrom stringr str_detect
#' @importFrom magrittr %>% 
#' @importFrom tibble tibble
#' @importFrom tidyr replace_na
#' @importFrom rlang .data
#' @importFrom readr parse_number
#' @examples
#' example
#'male_pop <- c(48875, 164390, 173551, 130297, 101143, 73615, 60594, 55175,
#'              49530, 46562, 39028, 27837, 22110, 18066, 15340, 13318,
#'              12002, 6424)
#'
#'female_pop <- c(47105, 159546, 168760, 119437, 92080, 70515, 58801, 53381,
#'                46757, 41164, 33811, 24121, 19315, 16319, 14058, 12302,
#'                11047, 5922)
#'
#'male_mx <- c(0.12427, 0.01639, 0.00274, 0.00167, 0.00251, 0.00380, 0.00382,
#'             0.00442, 0.00506, 0.00663, 0.00872, 0.01240, 0.01783, 0.02700,
#'             0.04126, 0.06785, 0.11287, 0.21015)
#'
#'female_mx <- c(0.11050, 0.01577, 0.00254, 0.00159, 0.00232, 0.00304, 0.00344,
#'               0.00370, 0.00418, 0.00492, 0.00592, 0.00831, 0.01182, 0.01942,
#'              0.03221, 0.05669, 0.09771, 0.19385)
#'
#'asfr <- c(0.199, 0.478, 0.418, 0.321, 0.163, 0.071, 0.028)
#'names(female_mx) <- seq(0, length.out = 18, by = 5)
#'
#'res <- movepop(
#'  initial_date = 1973.58, 
#'  desired_date = 1973,
#'  male_pop = male_pop, 
#'  female_pop = female_pop,
#'  male_mx  = male_mx,  
#'  female_mx  = female_mx,
#'  asfr = asfr,
#'  annual_net_migrants = -50000,
#'  age_out = "5-year",
#'  first_asfr_age = 15
#')
#'
#' @export

movepop <- function(initial_date        = NULL, 
                    desired_date        = NULL,
                    male_pop            = NULL,
                    female_pop          = NULL,
                    male_mx             = NULL,
                    female_mx           = NULL,
                    asfr                = NULL,
                    Age                 = NULL,
                    first_asfr_age      = NULL,
                    annual_net_migrants = 0,
                    age_out             = "single") {
  
  
  # first sanity check. negative values in data
  if(
  any(c(male_pop,
  female_pop,
  male_mx,
  female_mx,
  asfr,
  Age) < 0)) { 
    
    stop("One of your incoming vectors (population, asfr or mx) has negative values.")
    
    }
  
  
  
  
  # Here we detect the format of the outcoming data
  age_out <- match.arg(age_out, c("single", "5-year", "abridged"))
  
  # if age is not provided we try to find it in the names of at least one of 4 incoming vectors
  # we also demand that all these vectors should be the same length
  if(is.null(Age)) {
    
    # length of incoming vectors 
    lengths_vec <- c(
      length(male_pop),
      length(female_pop),
      length(male_mx),
      length(female_mx)
    )
    
    if(length(unique(lengths_vec)) != 1) {
      stop("All population and mortality vectors must have the same length (ages)")
    }
    
    # Now we can assume that if names are provided in more than one of the vectors
    # they are the same.
    Age <- unique(
      c(names2age(male_pop),
        names2age(female_pop),
        names2age(male_mx),
        names2age(female_mx)))
    
    # and we remove NA
    Age <- Age[!is.na(Age)]
    
    # now if age is provided
  } else {
    
    # length of incoming vecrtors 
    lengths_vec <- c(
      length(male_pop),
      length(female_pop),
      length(male_mx),
      length(female_mx),
      length(Age)
    )
    
    if(length(unique(lengths_vec)) != 1) {
      stop("All population and mortality vectors must have the same length (ages)")
    }
    
    if(is.character(Age)) { 
      
      Age <- parse_number(Age)
      
    }
    
    if(any(is.na(Age))) { 
      
      stop("You provided age in character format. Automatic conversion to numeric failed in one or more cases. Please reassign the Age vector with unambiguous values.")
      
      
      }
    
  }
  
  # Extreme case. No age provided non of the vectors have age as names too
  # we stop and demand Age
  if(length(Age) < 1) {
    
    stop("You should either provide an Age argument to the function or pass at least one of the mx / population vectors with Age as names.")
    
    }
  
  
  # incoming age type detection
  if(is_single(Age)) { 
    
    inc_age <- "single"
    
  } else if(is_abridged(Age)) { 
    
    inc_age <- "abridged"
    
  } else if(all(diff(Age[-length(Age)]) == 5)) { 
    
    inc_age <- "5-year"
    
  } else { 
    
    stop("Incoming Age is not in supporterd format (single, 5-year, or abridged). Please change the format.")
    
  }
  
  # collect all data into tibble
  data <- tibble(
    male_pop   = male_pop,
    female_pop = female_pop,
    male_mx    = male_mx,
    female_mx  = female_mx,
    age        = Age
  )
  
  # now we want to join asfr data if first_asfr_age is not provided we assume Age is in names
  # and again we have 2 options having ages as names in vector or now
  # if ages are there all good
  
  if(is.null(first_asfr_age)) { 
    
    if(sum(is.na(names2age(asfr))) == length(asfr)) { 
      
      stop("You did not provide the asfr starting age (first_asfr_age), and your asfr vector has no names. Please either provide first_asfr_age or assign correct ages as names to asfr vector.")
      
      }
    
    asfr_age <- names2age(asfr)
    asfr_age <- asfr_age[!is.na(asfr_age)]
    
  } else {
    
    if(inc_age == "single") { 
      
      asfr_age <- first_asfr_age:(first_asfr_age + length(asfr))
      
    } else { 
      
      asfr_age <- seq(from = first_asfr_age, length.out = length(asfr), by = 5)
      
    }  
  }
  
  asfr <- tibble(asfr = asfr,
                 age  = asfr_age)
  
  data <- data %>% 
    left_join(asfr, by = "age") %>% 
    mutate(asfr   = replace_na(.data$asfr, 0),
           births = .data$asfr * .data$female_pop,
           deaths = .data$male_mx * .data$male_pop + .data$female_mx * .data$female_pop)
  
  summaries <- data %>%
    summarise(total_births        = sum(.data$births),
              total_deaths        = sum(.data$deaths),
              total_pop_initial   = sum(.data$male_pop) + sum(.data$female_pop),
              crude_birth_rate    = .data$total_births / .data$total_pop_initial,
              crude_death_rate    = .data$total_deaths / .data$total_pop_initial,
              net_migration_rate  = annual_net_migrants / .data$total_pop_initial,
              growth_rate         = .data$crude_birth_rate - .data$crude_death_rate + .data$net_migration_rate,
              inc_age             = inc_age,
              tfr                 = ifelse(.data$inc_age %in% c("five_year", "abridged"), 5 * sum(.data$asfr), sum(.data$asfr)),
              desired_date        = desired_date,
              initial_date        = initial_date,
              time_diff           = .data$desired_date - .data$initial_date,
              total_pop_projected = .data$total_pop_initial * exp(.data$growth_rate * .data$time_diff),
              adjustment_factor   = .data$total_pop_projected / .data$total_pop_initial)
  
  
  projected_data <- data %>% 
    mutate(male_pop_projected   = summaries$adjustment_factor * .data$male_pop,
           female_pop_projected = summaries$adjustment_factor * .data$female_pop,
           both_sexes_projected = .data$male_pop_projected + .data$female_pop_projected,
           sex_ratio            = .data$male_pop_projected / .data$female_pop_projected)
  
  summaries_proj <- projected_data %>%
    summarise(age_group  = "All ages",
              both_sexes = sum(.data$both_sexes_projected),
              male       = sum(.data$male_pop_projected),
              female     = sum(.data$female_pop_projected),
              sex_ratio  = sum(.data$male_pop_projected) / sum(.data$female_pop_projected))
  
  projected_data1 <- projected_data %>%
    select("age", contains("proj"))
  
  # if user requested single there 3 options
  if(age_out == "single") {
    
    if(inc_age == "single") {
      
      projected_data1 <- projected_data1
      
    }
    
    # we have graduate with pclm
    if(inc_age %in% c("5-year", "abridged")) {
      
      projected_data1 <- projected_data1 %>% 
        select(-"age") %>%
        as.list() %>%
        map(~ { names(.x) <- projected_data1$age; .x}) %>%
        map(~ graduate_pclm(Value  = .x,
                            Age    = as.numeric(names(.x)),
                            OAG    = TRUE)) %>% 
        bind_rows(.id = "result") %>% 
        pivot_longer(-"result",
                     names_to   = "age",
                     values_to  = "val",
                     names_transform = list(age = as.numeric)) %>% 
        pivot_wider(names_from  = "result",
                    values_from = "val")
      
    }
  }
  
  # same for abridged 3 options
  if(age_out == "abridged") {
    
    if(inc_age == "single") {
      
      projected_data1 <- projected_data1 %>%
        select(-"age") %>%
        as.list() %>%
        map(~ { names(.x) <- projected_data1$age; .x}) %>% 
        map(~ single2abridged(.x)) %>% 
        bind_rows(.id = "result") %>% 
        pivot_longer(-"result",
                     names_to   = "age",
                     values_to  = "val",
                     names_transform = list(age = as.numeric)) %>% 
        pivot_wider(names_from  = "result",
                    values_from = "val")
      
    }
    
    if(inc_age == "5-year") {
      
      projected_data1 <- projected_data1 %>% 
        select(-"age") %>%
        as.list() %>%
        map(~ { names(.x) <- projected_data1$age; .x}) %>%
        map(~ graduate_pclm(Value  = .x,
                            Age    = as.numeric(names(.x)),
                            OAG    = TRUE)) %>% 
        map(~ single2abridged(.x)) %>% 
        bind_rows(.id = "result") %>% 
        pivot_longer(-"result",
                     names_to   = "age",
                     values_to  = "val",
                     names_transform = list(age = as.numeric)) %>% 
        pivot_wider(names_from  = "result",
                    values_from = "val") 
      
    }
    
    # the most complecated case.
    # first ungroup with pclm and then group together
    if(inc_age == "abridged") {
      
      projected_data1 <- projected_data1
      
    }
  }
  
  # same for 5 years
  if(age_out == "5-year") {
    
    if(inc_age == "single") {
      
      projected_data1 <- projected_data1 %>%
        select(-"age") %>%
        as.list() %>%
        map(~ { names(.x) <- projected_data1$age; .x}) %>% 
        map(~ groupAges(.x, N = 5)) %>% 
        bind_rows(.id = "result") %>% 
        pivot_longer(-"result",
                     names_to   = "age",
                     values_to  = "val",
                     names_transform = list(age = as.numeric)) %>% 
        pivot_wider(names_from  = "result",
                    values_from = "val")
      
    }
    
    if(inc_age == "5-year") {
      
      projected_data1 <- projected_data1
      
    }
    
    if(inc_age == "abridged") {
      
      projected_data1 <- projected_data1 %>% 
        mutate(age = replace(.data$age, .data$age %in% 0:1, 0)) %>%
        summarise(
          across(everything(), sum),
          .by = "age"
        )
      
    }
    
  }
  
  # Assemble results
  results <- list(
    initial_summaries   = summaries,
    projected_summaries = summaries_proj,
    data                = data,
    projected_data      = projected_data1)
  
  return(results)
}
