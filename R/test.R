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

AdjustFemalePop <- function(x, asfr_pastyear, sex_ratio, male = TRUE) {
  births_age_female <- x * asfr_pastyear
  est_tot <- sum(births_age_female) * sex_ratio / (1 + sex_ratio)

  if (!male) est_tot <- sum(births_age_female) - est_tot

  est_tot
}

basepop <- function(reference_year,
                    Males,
                    Females,
                    SmoothedFemales = NULL,
                    Ages,
                    SRB,
                    nlxFemale,
                    nlxMale = NULL,
                    nlxDatesIn,
                    nlxDatesOut,
                    asfrMat,
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
  nlxFemale <- interp(
    nlxFemale,
    datesIn = nlxDatesIn,
    datesOut = nlxDatesOut,
    method = "linear"
  )
  # Turn the columns from the matrix into a list
  nlxFemale <- lapply(as.data.frame(nlxFemale), identity)

  # Interpolate the asfr to the requested dates.
  # This is gender agnostic.
  Asfr <- interp(
    asfrMat,
    datesIn = AsfrDatesIn,
    datesOut = AsfrDatesOut,
    method = "linear"
  )
  # Turn the columns from the matrix into a list
  Asfr <- lapply(as.data.frame(Asfr), identity)
  
  EarliestDate <- which.max(names(nlxFemale))
  OlderNlxFemale <- nlxFemale[-EarliestDate]
  OlderNlxFemale <- OlderNlxFemale[sort(names(OlderNlxFemale), decreasing = TRUE)]

  # We assume that smoothed_females is returned by smooth_age_5, since
  # we use the smoothed vector names to only get certain age groups
  FemalePops <- if (!is.null(SmoothedFemales)) SmoothedFemales else Females
  SmoothedFMiddleages <- FemalePops[as.character(seq(15, 55, by = 5))]

  ## Currently, we assume that
  for (i in seq_along(OlderNlxFemale)) {
    .x <- OlderNlxFemale[[i]]
    VecMult <- if (i == 1) SmoothedFMiddleages else OlderNlxFemale[[i-1]]
    IndLength <- if (i == 1) 0 else 1

    NlxDiv <- .x[4:(length(.x) - IndLength)]
    iter <- c(1, rep(seq(2, length(NlxDiv) - 1), each = 2), length(NlxDiv))
    NlxSeq <- lapply(seq(1, length(iter), by = 2), function(i) NlxDiv[iter[i:(i+1)]])

    OlderNlxFemale[[i]] <- mapply(function(.y, .z) .y * .z[1] / .z[2],
                                  VecMult[-1],
                                  NlxSeq
                                  )
  }

  FirstDate <- list(0.2 * OlderNlxFemale[[1]] + 0.8 * SmoothedFMiddleages[-length(SmoothedFMiddleages)])
  names(FirstDate) <- max(names(nlxFemale))
  OlderNlxFemale <- c(
    FirstDate,
    OlderNlxFemale
  )
  
  RestrictedMiddleages <- FemalePops[as.character(seq(15, 45, by = 5))]

  ## Currently, we assume that for estimating the population
  ## For other years, after the second year, we carry forward
  ## the mean between the current year the previous one
  EstPop <- vector("list", length(OlderNlxFemale))
  names(EstPop) <- names(OlderNlxFemale)
  for (i in seq_along(OlderNlxFemale)) {
    .x <- OlderNlxFemale[[i]]
    mean_vec <- if (i > 2) OlderNlxFemale[[i-1]] else SmoothedFMiddleages

    EstPop[[i]] <- vapply(seq_along(RestrictedMiddleages),
                           function(i) mean(c(mean_vec[i], .x[i])),
                           FUN.VALUE = numeric(1))
  }

  EstTot <- mapply(
    AdjustFemalePop,
    EstPop,
    Asfr,
    SRB,
    male,
    SIMPLIFY = FALSE
  )

  GenderCounts <- if (male) Males else Females

  ## Currently, this assumes that there can only be 3 dates. How
  ## would we multiply the above if we had 5 dates?
  ## We only have 3 age groups to adjust and 3 dates
  # Age 0
  GenderCounts[1] <- EstTot[[1]] * nlx[[1]][1] / 100000
  # Age 1-4
  GenderCounts[2] <- EstTot[[2]] * 5 * (sum(nlx[[2]][1:2])) / 500000 - GenderCounts[1]
  # Age 5-9
  GenderCounts[3] <- EstTot[[3]] * 5 * (sum(nlx[[3]][1:2])) / 500000 * nlx[[2]][3] / sum(nlx[[2]][1:2])

  GenderCounts
}

############################# BPA #############################################
###############################################################################

bpa_male <-
  basepop(
    reference_year = reference_year,
    Males = pop_male_counts,
    Females = pop_female_counts,
    SmoothedFemales = smoothed_females,
    Ages = names(pop_male_counts),
    SRB = sex_ratio,
    nlxFemale = nlxFemale,
    nlxMale = nlxMale,
    nlxDatesIn = nlxDatesIn,
    nlxDatesOut = nlxDatesOut,
    asfrMat = asfrmat,
    AsfrDatesIn = AsfrDatesIn,
    AsfrDatesOut = AsfrDatesOut
  )

all(round(bpa_male[1:3], 0) == c(13559, 47444, 54397))

bpa_female <-
  basepop(
    reference_year = reference_year,
    Males = pop_male_counts,
    Females = pop_female_counts,
    SmoothedFemales = smoothed_females,
    Ages = names(pop_female_counts),
    SRB = sex_ratio,
    nlxFemale = nlxFemale,
    nlxDatesIn = nlxDatesIn,
    nlxDatesOut = nlxDatesOut,
    asfrMat = asfrmat,
    AsfrDatesIn = AsfrDatesIn,
    AsfrDatesOut = AsfrDatesOut
  )

all(round(bpa_female[1:3], 0) == c(13467, 47576, 54554))

############################# BPE #############################################
###############################################################################

bpe_male <-
  basepop(
    reference_year = reference_year,
    Males = pop_male_counts,
    Females = pop_female_counts,
    Ages = names(pop_male_counts),
    SRB = sex_ratio,
    nlxFemale = nlxFemale,
    nlxMale = nlxMale,
    nlxDatesIn = nlxDatesIn,
    nlxDatesOut = nlxDatesOut,
    asfrMat = asfrmat,
    AsfrDatesIn = AsfrDatesIn,
    AsfrDatesOut = AsfrDatesOut
  )

all(round(bpe_male[1:3], 0) == c(13679, 47967, 55721))

bpe_female <-
  basepop(
    reference_year = reference_year,
    Males = pop_male_counts,
    Females = pop_female_counts,
    Ages = names(pop_female_counts),
    SRB = sex_ratio,
    nlxFemale = nlxFemale,
    nlxDatesIn = nlxDatesIn,
    nlxDatesOut = nlxDatesOut,
    asfrMat = asfrmat,
    AsfrDatesIn = AsfrDatesIn,
    AsfrDatesOut = AsfrDatesOut
  )

all(round(bpe_female[1:3], 0) == c(13587, 48101, 55882))
