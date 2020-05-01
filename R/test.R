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
nlxDatesIn <- c(1977.31, 1986.50)
nlxDatesOut <- c(1985.71, 1983.71, 1978.71)

nlxMale <- matrix(c(87732,
                    304435,
                    361064,
                    88451,
                    310605,
                    370362
                    ),
                  nrow = 3, ncol = 2)

nlxFemale <- matrix(c(89842,
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
                      217419,
                      90478,
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
                      234587),
                    nrow = 12,
                    ncol = 2)

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
                                 Age = as.numeric(Ages),
                                 method = "Arriaga",
                                 OAG = TRUE,
                                 young.tail = "Original")

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

basepop <- function(reference_year,
                    pop_male_counts, # Males
                    pop_female_counts, # Females
                    smoothed_females,
                    Ages,
                    sex_ratio, # SRB
                    nlxFemale,
                    nlxDatesIn,
                    nlxDatesOut,
                    nlxMale = NULL,
                    asfrmat,
                    AsfrDatesIn,
                    AsfrDatesOut) {

  # We calculate basepop for males **ONLY** if the nlxMale
  # argument is not null
  nlxGender <- if (!is.null(nlxMale)) nlxMale else nlxFemale
  male <- if (!is.null(nlxMale)) TRUE else FALSE

  # Interpolate the gender specific nlx to the requested
  # dates out
  nlx <- interp(
    nlxGender,
    datesIn = nlxDatesIn,
    datesOut = nlxDatesOut,
    method = "linear"
  )
  # Turn the columns from the matrix into a list
  nlx <- lapply(as.data.frame(nlx), identity)

  # Interpolate only the female nlx to the requested dates.
  # This is independent of the nlx above, since to apply
  # basepop to either gender you need to have the female nlx
  nlx_female <- interp(
    nlxFemale,
    datesIn = nlxDatesIn,
    datesOut = nlxDatesOut,
    method = "linear"
  )
  # Turn the columns from the matrix into a list
  nlx_female <- lapply(as.data.frame(nlx_female), identity)

  # Interpolate the asfr to the requested dates.
  # This is gender agnostic.
  asfr <- interp(
    asfrmat,
    datesIn = AsfrDatesIn,
    datesOut = AsfrDatesOut,
    method = "linear"
  )
  # Turn the columns from the matrix into a list
  asfr <- lapply(as.data.frame(asfr), identity)

  
  EarliestDate <- which.max(names(nlx_female))
  OlderNlxFemale <- nlx_female[-EarliestDate]
  OlderNlxFemale <- OlderNlxFemale[sort(names(OlderNlxFemale), decreasing = TRUE)]

  # We assume that smoothed_females is returned by smooth_age_5, since
  # we use the smoothed vector names to only get certain age groups
  smoothed_f_middleages <- smoothed_females[as.character(seq(15, 55, by = 5))]

  ## Currently, we assume that
  for (i in seq_along(OlderNlxFemale)) {
    .x <- OlderNlxFemale[[i]]
    vec_mult <- if (i == 1) smoothed_f_middleages else OlderNlxFemale[[i-1]]
    ind_length <- if (i == 1) 0 else 1

    nlx_div <- .x[4:(length(.x) - ind_length)]
    iter <- c(1, rep(seq(2, length(nlx_div) - 1), each = 2), length(nlx_div))
    nlx_seq <- lapply(seq(1, length(iter), by = 2), function(i) nlx_div[iter[i:(i+1)]])

    OlderNlxFemale[[i]] <- mapply(function(.y, .z) .y * .z[1] / .z[2],
                                  vec_mult[-1],
                                  nlx_seq
                                  )
  }

  FirstDate <- list(0.2 * OlderNlxFemale[[1]] + 0.8 * smoothed_f_middleages[-length(smoothed_f_middleages)])
  names(FirstDate) <- max(names(nlx_female))
  OlderNlxFemale <- c(
    FirstDate,
    OlderNlxFemale
  )
  
  restricted_middleages <- smoothed_females[as.character(seq(15, 45, by = 5))]

  ## Currently, we assume that for estimating the population
  ## For other years, after the second year, we carry forward
  ## the mean between the current year the previous one
  est_pop <- vector("list", length(OlderNlxFemale))
  names(est_pop) <- names(OlderNlxFemale)
  for (i in seq_along(OlderNlxFemale)) {
    .x <- OlderNlxFemale[[i]]
    mean_vec <- if (i > 2) OlderNlxFemale[[i-1]] else smoothed_f_middleages

    est_pop[[i]] <- vapply(seq_along(restricted_middleages),
                           function(i) mean(c(mean_vec[i], .x[i])),
                           FUN.VALUE = numeric(1))
  }
  ## End intermediate calculations

  est_tot <- mapply(
    adjust_female_pop,
    est_pop,
    asfr,
    sex_ratio,
    male,
    SIMPLIFY = FALSE
  )

  pop_gender_counts <- if (male) pop_male_counts else pop_female_counts

  ## Currently, this assumes that there can only be 3 dates. How
  ## would we multiply the above if we had 5 dates?
  ## We only have 3 age groups to adjust and 3 dates

  # Age 0
  pop_gender_counts[1] <- est_tot[[1]] * nlx[[1]][1] / 100000
  # Age 1-4
  pop_gender_counts[2] <- est_tot[[2]] * 5 * (sum(nlx[[2]][1:2])) / 500000 - pop_gender_counts[1]
  # Age 5-9
  pop_gender_counts[3] <- est_tot[[3]] * 5 * (sum(nlx[[3]][1:2])) / 500000 * nlx[[2]][3] / sum(nlx[[2]][1:2])

  pop_gender_counts
}

arriaga_male <-
  basepop(
    reference_year = reference_year,
    pop_male_counts = pop_male_counts,
    pop_female_counts = pop_female_counts,
    smoothed_females = smoothed_females,
    Ages = names(pop_male_counts),
    sex_ratio = sex_ratio,
    nlxFemale = nlxFemale,
    nlxMale = nlxMale,
    nlxDatesIn = nlxDatesIn,
    nlxDatesOut = nlxDatesOut,
    asfrmat = asfrmat,
    AsfrDatesIn = AsfrDatesIn,
    AsfrDatesOut = AsfrDatesOut
  )

all(round(arriaga_male[1:3], 0) == c(13559, 47444, 54397))

arriaga_female <-
  basepop(
    reference_year = reference_year,
    pop_male_counts = pop_male_counts,
    pop_female_counts = pop_female_counts,
    smoothed_females = smoothed_females,
    Ages = names(pop_female_counts),
    sex_ratio = sex_ratio,
    nlxFemale = nlxFemale,
    nlxDatesIn = nlxDatesIn,
    nlxDatesOut = nlxDatesOut,
    asfrmat = asfrmat,
    AsfrDatesIn = AsfrDatesIn,
    AsfrDatesOut = AsfrDatesOut
  )

all(round(arriaga_female[1:3], 0) == c(13467, 47576, 54554))

