library(DemoTools)
## library(tibble)

## Ages <- seq(0, 80, by = 5)

## female_pop_test <- c(
##   58366,
##   55812, 
##   35268, 
##   33672, 
##   31352, 
##   33038, 
##   24029, 
##   16120, 
##   14679, 
##   8831, 
##   9289, 
##   4172, 
##   6174, 
##   2715, 
##   3344, 
##   1455, 
##   4143 
## )

## female_smoothed_arriaga <- c(
##   61688,
##   55882, 
##   35268, 
##   33672, 
##   31352, 
##   33038, 
##   24029, 
##   16120, 
##   14679, 
##   8831, 
##   9289, 
##   4172, 
##   6174, 
##   2715, 
##   3344, 
##   1455, 
##   4143 
## )

## male_pop_test <-
##   c(58422,
##     55639, 
##     37514, 
##     29398, 
##     27187, 
##     27770, 
##     20920, 
##     16973, 
##     14999, 
##     11330, 
##     10415, 
##     6164, 
##     7330, 
##     3882, 
##     3882, 
##     1840, 
##     4200
##     )

## male_smoothed_arriaga <-
##   c(61002,
##     54397, 
##     34738, 
##     32174, 
##     29399, 
##     25558, 
##     20622, 
##     17271, 
##     14459, 
##     11870, 
##     9143, 
##     7436, 
##     6288, 
##     4924, 
##     3552, 
##     2170, 
##     4200
##     )

## smooth_pop <-
##   round(
##     smooth_age_5(Value = female_pop_test, Age = Ages, method = "Arriaga", OAG = TRUE),
##     3
##   )

## dt <-
##   tibble(ag = Ages,
##          pop = female_pop_test,
##          smooth_pop = as.integer(round(smooth_pop, 0)),
##          smoothed_arriaga = as.integer(female_smoothed_arriaga),
##          equal = smooth_pop - smoothed_arriaga
##          )


# Only for males
# (1) Reference year of the counts
reference_year <- 1986.21

# (2) Reported population by 5-year age groups and sex in the base year (Include unknowns).
pop_male_counts <- c(`0` = 11684,
                     `1` = 46738,
                     `5` = 55639,
                     `10` = 37514,
                     `15` = 29398, 
                     `20` = 27187,
                     `25` = 27770,
                     `30` = 20920,
                     `35` = 16973,
                     `40` = 14999, 
                     `45` = 11330,
                     `50` = 10415,
                     `55` = 6164,
                     `60` = 7330,
                     `65` = 3882, 
                     `70` = 3882,
                     `75` = 1840,
                     `80` = 4200
                     )

pop_female_counts <- c(`0` = 11673,
                       `1` = 46693,
                       `5` = 55812,
                       `10` = 35268,
                       `15` = 33672, 
                       `20` = 31352,
                       `25` = 33038,
                       `30` = 24029,
                       `35` = 16120,
                       `40` = 14679, 
                       `45` = 8831,
                       `50` = 9289,
                       `55` = 4172,
                       `60` = 6174,
                       `65` = 2715, 
                       `70` = 3344,
                       `75` = 1455,
                       `80` = 4143)

# (4) Sex ratio at birth (m/f)
sex_ratio <- 1.0300

# (6) The male and female nLx functions for ages under 1 year, 1 to 4 years, and 5 to 9
# years, pertaining to an earlier and later date
nlx_earlydate <- 1977.31
nlx_male_earlydate <- c(87732,
                        304435,
                        361064)

nlx_female_earlydate <- c(89842,
                          314521,
                          372681,
                          353053,
                          340650,
                          326588,
                          311481,
                          295396,
                          278646,
                          261260,
                          241395,
                          217419)
nlx_latedate <- 1986.50
nlx_male_latedate <- c(88451, 310605, 370362)
nlx_female_latedate <- c(90478,
                         320755,
                         382531,
                         364776, 
                         353538, 
                         340687, 
                         326701, 
                         311573, 
                         295501, 
                         278494, 
                         258748,
                         234587)

# (7) A set of age-specific fertility rates pertaining to an earlier and later
# date
asfr_earlydate <- 1977.81
asfr_interp_earlydate <- c("15-19" = 0.2000,
                           "20-24" = 0.3000,
                           "25-29" = 0.3000, 
                           "30-34" = 0.2500, 
                           "35-39" = 0.2000, 
                           "40-44" = 0.1500, 
                           "45-49" = 0.0500)

asfr_latedate <- 1985.71
asfr_interp_latedate <- c("15-19" = 0.1500,
                          "20-24" = 0.2000,
                          "25-29" = 0.2750, 
                          "30-34" = 0.2250, 
                          "35-39" = 0.1750, 
                          "40-44" = 0.1250, 
                          "45-49" = 0.0500)

## reference_year
## pop_gender_counts <- pop_counts_male
ages <- names(pop_female_counts)
## sex_ratio
## nlx_gender_earlydate <- nlx_male_earlydate
## nlx_earlydate
## nlx_gender_latedate <- nlx_male_latedate
## nlx_latedate
## asfr_interp_earlydate
## asfr_earlydate
## asfr_interp_latedate
## asfr_latedate

adjust_female_pop <- function(x, asfr_pastyear, sex_ratio, male = TRUE) {
  births_age_female <- x * asfr_pastyear
  est_tot <- sum(births_age_female) * sex_ratio / (1 + sex_ratio)

  if (!male) est_tot <- sum(births_age_female) - est_tot

  est_tot
}

estimates_nlx <- function(nlx_m_t1,
                          nlx_t1,
                          nlx_m_t2,
                          nlx_t2,
                          t0,
                          t_before = 0.5) {
  res <- nlx_m_t1*(t0-t_before-nlx_t2)/(nlx_t1-nlx_t2)+nlx_m_t2*(nlx_t1-t0+t_before)/(nlx_t1-nlx_t2)

  res
}

smooth_youngtail_5_arriaga <- function(reference_year,
                                       pop_male_counts,
                                       pop_female_counts,
                                       ages,
                                       sex_ratio,
                                       nlx_male_earlydate,
                                       nlx_female_earlydate,
                                       nlx_earlydate,
                                       nlx_male_latedate,
                                       nlx_female_latedate,
                                       nlx_latedate,
                                       asfr_interp_earlydate,
                                       asfr_earlydate,
                                       asfr_interp_latedate,
                                       asfr_latedate,
                                       male = TRUE) {

  if (male) {
    nlx_gender_earlydate <- nlx_male_earlydate
    nlx_gender_latedate <- nlx_male_latedate    
  } else {
    nlx_gender_earlydate <- nlx_female_earlydate
    nlx_gender_latedate <- nlx_female_latedate    
  }

  past_year <- c("past_year" = 0.5, "twohalf_before" = 2.5, "sevenhalf_before" = 7.5)
  # Matches with female
  # Matches with male
  nlx <- lapply(past_year,
                estimates_nlx,
                nlx_m_t1 = nlx_gender_earlydate,
                nlx_t1 = nlx_earlydate,
                nlx_m_t2 = nlx_gender_latedate,
                nlx_t2 = nlx_latedate,
                t0 = reference_year
                )

  # Matches with female
  nlx_female <- lapply(past_year,
                       estimates_nlx,
                       nlx_m_t1 = nlx_female_earlydate,
                       nlx_t1 = nlx_earlydate,
                       nlx_m_t2 = nlx_female_latedate,
                       nlx_t2 = nlx_latedate,
                       t0 = reference_year
                       )

  # Matches with baseline asfr
  asfr <- lapply(past_year,
                 estimates_nlx,
                 nlx_m_t1 = asfr_interp_earlydate,
                 nlx_t1 = asfr_earlydate,
                 nlx_m_t2 = asfr_interp_latedate,
                 nlx_t2 = asfr_latedate,
                 t0 = reference_year
                 )

  ## Intermediatte calculations
  smoothed_middleages <- smooth_age_5(Value = pop_female_counts,
                                      Age = as.numeric(ages),
                                      method = "Arriaga",
                                      OAG = TRUE)[as.character(seq(15, 55, by = 5))]

  ## smoothed_middleages[1] <- 34721

  nlx_div <- nlx_female$twohalf_before[4:length(nlx_female$twohalf_before)]
  iter <- c(1, rep(seq(2, length(nlx_div) - 1), each = 2), length(nlx_div))
  nlx_seq <- lapply(seq(1, length(iter), by = 2), function(i) nlx_div[iter[i:(i+1)]])

  # Matches with baseline except for first age group
  est_pop_middledate <- mapply(function(x, y) x * y[1] / y[2],
                             smoothed_middleages[-1],
                             nlx_seq
                             )

  # Matches with baseline except for first age group
  est_pop_earlydate <- 0.2 * est_pop_middledate + 0.8 * smoothed_middleages[-length(smoothed_middleages)]

  nlx_div <- nlx_female$sevenhalf_before[4:(length(nlx_female$sevenhalf_before) - 1)]
  iter <- c(1, rep(seq(2, length(nlx_div) - 1), each = 2), length(nlx_div))
  nlx_seq <- lapply(seq(1, length(iter), by = 2), function(i) nlx_div[iter[i:(i+1)]])

  est_pop_latedate <- mapply(function(x, y) x * y[1] / y[2],
                             est_pop_middledate[-1],
                             nlx_seq
                             )

  restricted_middleages <- smoothed_middleages[as.character(seq(15, 45, by = 5))]

  # All est_pop_* match except for the first age group
  est_pop_past_year <- vapply(seq_along(restricted_middleages),
                              function(i) mean(c(smoothed_middleages[i], est_pop_earlydate[i])),
                              FUN.VALUE = numeric(1))

  est_pop_twohalfbefore <- vapply(seq_along(restricted_middleages),
                                  function(i) mean(c(smoothed_middleages[i], est_pop_middledate[i])),
                                  FUN.VALUE = numeric(1))

  est_pop_sevenhalfbefore <- vapply(seq_along(restricted_middleages),
                                    function(i) mean(c(est_pop_middledate[i], est_pop_latedate[i])),
                                    FUN.VALUE = numeric(1))

  ## End intermediate calculations

  est_pop_years <- list("past_year" = est_pop_past_year,
                        "twohalf_before" = est_pop_twohalfbefore,
                        "sevenhalf_before" = est_pop_sevenhalfbefore)

  est_tot <- mapply(
    adjust_female_pop,
    est_pop_years,
    list(asfr$past_year, asfr$twohalf_before, asfr$sevenhalf_before),
    sex_ratio,
    male,
    SIMPLIFY = FALSE
  )

  pop_gender_counts <- if (male) pop_male_counts else pop_female_counts

  # Age 0
  pop_gender_counts[1] <- est_tot$past_year * nlx$past_year[1] / 100000
  # Age 1-4
  pop_gender_counts[2] <- est_tot$twohalf_before * 5 * (sum(nlx$twohalf_before[1:2])) / 500000 - pop_gender_counts[1]
  # Age 5-9
  pop_gender_counts[3] <- est_tot$sevenhalf_before * 5 * (sum(nlx$sevenhalf_before[1:2])) / 500000 * nlx$twohalf_before[3] / sum(nlx$twohalf_before[1:2])

  pop_gender_counts
}

arriaga_male <-
  smooth_youngtail_5_arriaga(
    reference_year = reference_year,
    pop_male_counts = pop_male_counts,
    pop_female_counts = pop_female_counts,
    ages = names(pop_male_counts),
    sex_ratio = sex_ratio,
    nlx_female_earlydate = nlx_female_earlydate,
    nlx_male_earlydate = nlx_male_earlydate,
    nlx_earlydate = nlx_earlydate,
    nlx_female_latedate = nlx_female_latedate,
    nlx_male_latedate = nlx_male_latedate,
    nlx_latedate = nlx_latedate,
    asfr_interp_earlydate = asfr_interp_earlydate,
    asfr_earlydate = asfr_earlydate,
    asfr_interp_latedate = asfr_interp_latedate,
    asfr_latedate = asfr_latedate,
    male = TRUE
  )

## all(round(arriaga_male[1:3], 0) == c(13559, 47444, 54397))

difference <- arriaga_male[1:3] - c(13559, 47444, 54397)
difference

arriaga_female <-
  smooth_youngtail_5_arriaga(
    reference_year = reference_year,
    pop_male_counts = pop_male_counts,
    pop_female_counts = pop_female_counts,
    ages = names(pop_female_counts),
    sex_ratio = sex_ratio,
    nlx_female_earlydate = nlx_female_earlydate,
    nlx_male_earlydate = nlx_male_earlydate,
    nlx_earlydate = nlx_earlydate,
    nlx_female_latedate = nlx_female_latedate,
    nlx_male_latedate = nlx_male_latedate,
    nlx_latedate = nlx_latedate,
    asfr_interp_earlydate = asfr_interp_earlydate,
    asfr_earlydate = asfr_earlydate,
    asfr_interp_latedate = asfr_interp_latedate,
    asfr_latedate = asfr_latedate,
    male = FALSE
  )

## all(round(arriaga_female[1:3], 0) == c(13467, 47576, 54554))

difference <- arriaga_female[1:3] - c(13467, 47576, 54554)
difference
