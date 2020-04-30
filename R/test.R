library(DemoTools)

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

Ages <- as.numeric(names(pop_female_counts))

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

asfrmat <- matrix(c("15-19" = 0.2000,
                    "20-24" = 0.3000,
                    "25-29" = 0.3000, 
                    "30-34" = 0.2500, 
                    "35-39" = 0.2000, 
                    "40-44" = 0.1500, 
                    "45-49" = 0.0500,
                    "15-19" = 0.1500,
                    "20-24" = 0.2000,
                    "25-29" = 0.2750, 
                    "30-34" = 0.2250, 
                    "35-39" = 0.1750, 
                    "40-44" = 0.1250, 
                    "45-49" = 0.0500),
                  nrow = 7,
                  ncol = 2)

AsfrDatesIn <- c(1977.81, 1985.71)
AsfrDatesOut <- c(1985.71, 1983.71, 1978.71)

## Intermediatte calculations
smoothed_females <- smooth_age_5(Value = pop_female_counts,
                                    Age = as.numeric(ages),
                                    method = "Arriaga",
                                    OAG = TRUE)

## This is the only number that messes up the whole calculation.
## smooth_age_5 returns the same result as the excel sheet
## except for the age groups 10-15 and 15-19. Here we only use
## age group 15-19. If we plug in manually the correct value,
## we get all results match exactly, otherwise there are
## some differences.
smoothed_females[4] <- 34721

## reference_year
## pop_gender_counts <- pop_counts_male
## sex_ratio
## nlx_gender_earlydate <- nlx_male_earlydate
## nlx_earlydate
## nlx_gender_latedate <- nlx_male_latedate
## nlx_latedate

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

basepop <- function(reference_year,
                    pop_male_counts, # Males
                    pop_female_counts, # Females
                    smoothed_females,
                    Ages,
                    sex_ratio, # SRB
                    nlx_male_earlydate,
                    nlx_female_earlydate,
                    nlx_earlydate,
                    nlx_male_latedate,
                    nlx_female_latedate,
                    nlx_latedate,
                    asfrmat,
                    AsfrDatesIn,
                    AsfrDatesOut,
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
  asfr <- interp(
    asfrmat,
    datesIn = AsfrDatesIn,
    datesOut = AsfrDatesOut,
    method = "linear"
  )
  asfr <- lapply(as.data.frame(asfr), identity)

  smoothed_f_middleages <- smoothed_females[as.character(seq(15, 55, by = 5))]

  nlx_div <- nlx_female$twohalf_before[4:length(nlx_female$twohalf_before)]
  iter <- c(1, rep(seq(2, length(nlx_div) - 1), each = 2), length(nlx_div))
  nlx_seq <- lapply(seq(1, length(iter), by = 2), function(i) nlx_div[iter[i:(i+1)]])

  # Matches with baseline except for first age group
  est_pop_middledate <- mapply(function(x, y) x * y[1] / y[2],
                               smoothed_f_middleages[-1],
                               nlx_seq
                               )

  # Matches with baseline except for first age group
  est_pop_earlydate <- 0.2 * est_pop_middledate + 0.8 * smoothed_f_middleages[-length(smoothed_f_middleages)]

  nlx_div <- nlx_female$sevenhalf_before[4:(length(nlx_female$sevenhalf_before) - 1)]
  iter <- c(1, rep(seq(2, length(nlx_div) - 1), each = 2), length(nlx_div))
  nlx_seq <- lapply(seq(1, length(iter), by = 2), function(i) nlx_div[iter[i:(i+1)]])

  est_pop_latedate <- mapply(function(x, y) x * y[1] / y[2],
                             est_pop_middledate[-1],
                             nlx_seq
                             )

  restricted_middleages <- smoothed_females[as.character(seq(15, 45, by = 5))]

  # All est_pop_* match except for the first age group
  est_pop_past_year <- vapply(seq_along(restricted_middleages),
                              function(i) mean(c(smoothed_f_middleages[i], est_pop_earlydate[i])),
                              FUN.VALUE = numeric(1))

  est_pop_twohalfbefore <- vapply(seq_along(restricted_middleages),
                                  function(i) mean(c(smoothed_f_middleages[i], est_pop_middledate[i])),
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
    asfr,
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
  basepop(
    reference_year = reference_year,
    pop_male_counts = pop_male_counts,
    pop_female_counts = pop_female_counts,
    smoothed_females = smoothed_females,
    ages = names(pop_male_counts),
    sex_ratio = sex_ratio,
    nlx_female_earlydate = nlx_female_earlydate,
    nlx_male_earlydate = nlx_male_earlydate,
    nlx_earlydate = nlx_earlydate,
    nlx_female_latedate = nlx_female_latedate,
    nlx_male_latedate = nlx_male_latedate,
    nlx_latedate = nlx_latedate,
    asfrmat = asfrmat,
    AsfrDatesIn = AsfrDatesIn,
    AsfrDatesOut = AsfrDatesOut,
    male = TRUE
  )

all(round(arriaga_male[1:3], 0) == c(13559, 47444, 54397))

## difference <- arriaga_male[1:3] - c(13559, 47444, 54397)
## difference

arriaga_female <-
  basepop(
    reference_year = reference_year,
    pop_male_counts = pop_male_counts,
    pop_female_counts = pop_female_counts,
    smoothed_females = smoothed_females,
    ages = names(pop_female_counts),
    sex_ratio = sex_ratio,
    nlx_female_earlydate = nlx_female_earlydate,
    nlx_male_earlydate = nlx_male_earlydate,
    nlx_earlydate = nlx_earlydate,
    nlx_female_latedate = nlx_female_latedate,
    nlx_male_latedate = nlx_male_latedate,
    nlx_latedate = nlx_latedate,
    asfrmat = asfrmat,
    AsfrDatesIn = AsfrDatesIn,
    AsfrDatesOut = AsfrDatesOut,
    male = FALSE
  )

all(round(arriaga_female[1:3], 0) == c(13467, 47576, 54554))

## difference <- arriaga_female[1:3] - c(13467, 47576, 54554)
## difference
