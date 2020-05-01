# (1) RefYear
refYear <- 1986.21

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
nlxDatesIn <- c(1977.31, 1986.50)

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

smoothed_females <- smooth_age_5(Value = pop_female_counts,
                                 Age = as.numeric(names(pop_female_counts)),
                                 method = "Arriaga",
                                 OAG = TRUE,
                                 young.tail = "Original")

## This is the only number that messes up the whole calculation.
## smooth_age_5 returns the same result as the PASS excel sheet
## except for the age groups 10-15 and 15-19. Here we only use
## age group 15-19. If we plug in manually the correct value,
## we get all results match exactly, otherwise there are
## some differences.
smoothed_females[4] <- 34721

test_that("basepop - bpa matches the expected result from PASS", {

  bpa_male <-
    basepop(
      refYear = refYear,
      Males = pop_male_counts,
      Females = pop_female_counts,
      SmoothedFemales = smoothed_females,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxMale = nlxMale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "m",
      method = "linear"      
    )

  expect_true(all(round(bpa_male[1:3], 0) == c(13559, 47444, 54397)))

  bpa_female <-
    basepop(
      refYear = refYear,
      Females = pop_female_counts,
      SmoothedFemales = smoothed_females,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"      
    )

  expect_true(all(round(bpa_female[1:3], 0) == c(13467, 47576, 54554)))
})

test_that("basepop - bpe matches the expected result from PASS", {

  bpe_male <-
    basepop(
      refYear = refYear,
      Males = pop_male_counts,
      Females = pop_female_counts,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxMale = nlxMale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "m",
      method = "linear"      
    )

  expect_true(all(round(bpe_male[1:3], 0) == c(13679, 47967, 55721)))

  bpe_female <-
    basepop(
      refYear = refYear,
      Females = pop_female_counts,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"      
    )

  expect_true(all(round(bpe_female[1:3], 0) == c(13587, 48101, 55882)))
})

test_that("basepop Sex argument controls the result regardless of the gender arguments ", {

  bpa_male <-
    basepop(
      refYear = refYear,
      Males = pop_male_counts,
      Females = pop_female_counts,
      SmoothedFemales = smoothed_females,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxMale = nlxMale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"      
    )

  bpa_female <-
    basepop(
      refYear = refYear,
      Females = pop_female_counts,
      SmoothedFemales = smoothed_females,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"      
    )

  expect_true(all(bpa_male == bpa_female))

  # Works as well for smoothing
  bpe_male <-
    basepop(
      refYear = refYear,
      Males = pop_male_counts,
      Females = pop_female_counts,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxMale = nlxMale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"      
    )

  bpe_female <-
    basepop(
      refYear = refYear,
      Females = pop_female_counts,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"
    )

  expect_true(all(bpe_male == bpe_female))
})

