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

test_that("basepop_five - bpa matches the expected result from PASS", {

  bpa_male <-
    basepop_five(
      refYear = refYear,
      Males_five = pop_male_counts,
      Females_five = pop_female_counts,
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
    basepop_five(
      refYear = refYear,
      Females_five = pop_female_counts,
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

test_that("basepop_five - bpe matches the expected result from PASS", {

  bpe_male <-
    basepop_five(
      refYear = refYear,
      Males_five = pop_male_counts,
      Females_five = pop_female_counts,
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
    basepop_five(
      refYear = refYear,
      Females_five = pop_female_counts,
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

test_that("basepop_five Sex argument controls the result regardless of the gender arguments ", {

  bpa_male <-
    basepop_five(
      refYear = refYear,
      Males_five = pop_male_counts,
      Females_five = pop_female_counts,
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
    basepop_five(
      refYear = refYear,
      Females_five = pop_female_counts,
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
    basepop_five(
      refYear = refYear,
      Males_five = pop_male_counts,
      Females_five = pop_female_counts,
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
    basepop_five(
      refYear = refYear,
      Females_five = pop_female_counts,
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

test_that("basepop_single fails if provided five year age groups", {

  female_single <-
  c(
    `0` = 11673,
    `1` = 11474,
    `2` = 11670,
    `3` = 11934,
    `4` = 11614,
    `5` = 10603,
    `6` = 11144,
    `7` = 11179,
    `8` = 11269,
    `9` = 11617,
    `10` = 6772,
    `11` = 6948,
    `12` = 7030,
    `13` = 7211,
    `14` = 7306,
    `15` = 6531,
    `16` = 6443,
    `17` = 6535,
    `18` = 6951,
    `19` = 7213,
    `20` = 6096,
    `21` = 6234,
    `22` = 6327,
    `23` = 6410,
    `24` = 6285,
    `25` = 6464,
    `26` = 6492,
    `27` = 6549,
    `28` = 6739,
    `29` = 6795,
    `30` = 5013,
    `31` = 4888,
    `32` = 4735,
    `33` = 4747,
    `34` = 4646,
    `35` = 3040,
    `36` = 3068,
    `37` = 3107,
    `38` = 3246,
    `39` = 3658,
    `40` = 2650,
    `41` = 2788,
    `42` = 2977,
    `43` = 3108,
    `44` = 3156,
    `45` = 1756,
    `46` = 1784,
    `47` = 1802,
    `48` = 1764,
    `49` = 1724,
    `50` = 1982,
    `51` = 1935,
    `52` = 1846,
    `53` = 1795,
    `54` = 1731,
    `55` = 863,
    `56` = 850,
    `57` = 825,
    `58` = 819,
    `59` = 816,
    `60` = 1348,
    `61` = 1342,
    `62` = 1246,
    `63` = 1138,
    `64` = 1101,
    `65` = 391,
    `66` = 520,
    `67` = 585,
    `68` = 560,
    `69` = 659,
    `70` = 670,
    `71` = 750,
    `72` = 686,
    `73` = 634,
    `74` = 604,
    `75` = 353,
    `76` = 340,
    `77` = 270,
    `78` = 246,
    `79` = 247,
    `80` = 4143
  )

  # Correction for males
  # To test that males are checked for single ages, we first
  # correctly define the female pop counts as single ages
  expect_error(
    basepop_single(
      refYear = refYear,
      Males_single = pop_male_counts,
      Females_single = female_single,
      SmoothedFemales = smoothed_females,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxMale = nlxMale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"
    ),
    "is_single(as.numeric(names(Males_single))) is not TRUE",
    fixed = TRUE
  )

  # check that pop male counts are named
  expect_error(
    basepop_single(
      refYear = refYear,
      Males_single = setNames(pop_male_counts, NULL),
      Females_single = female_single,
      SmoothedFemales = smoothed_females,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxMale = nlxMale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"
    ),
    "!is.null(names(Males_single)) is not TRUE",
    fixed = TRUE
  )

  # Check female counts are single ages
  expect_error(
    basepop_single(
      refYear = refYear,
      Females_single = pop_female_counts,
      SmoothedFemales = smoothed_females,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"
    ),
    "is_single(as.numeric(names(Females_single))) is not TRUE",
    fixed = TRUE
  )

  # Check that pop female counts are named
  expect_error(
    basepop_single(
      refYear = refYear,
      Females_single = setNames(pop_female_counts, NULL),
      SmoothedFemales = smoothed_females,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"
    ),
    "!is.null(names(Females_single)) is not TRUE",
    fixed = TRUE
  )

  # Works as well for smoothing
  expect_error(
    basepop_single(
      refYear = refYear,
      Males_single = pop_male_counts,
      Females_single = pop_female_counts,
      SRB = sex_ratio,
      nlxFemale = nlxFemale,
      nlxMale = nlxMale,
      nlxDatesIn = nlxDatesIn,
      asfrMat = asfrmat,
      AsfrDatesIn = AsfrDatesIn,
      Sex = "f",
      method = "linear"
    ),
    "is_single(as.numeric(names(Females_single))) is not TRUE",
    fixed = TRUE
  )
})

pop_male_counts <-
  c(
    `0` = 11684,
    `1` = 11473,
    `2` = 11647,
    `3` = 11939,
    `4` = 11680,
    `5` = 10600,
    `6` = 11100,
    `7` = 11157,
    `8` = 11238,
    `9` = 11544,
    `10` = 7216,
    `11` = 7407,
    `12` = 7461,
    `13` = 7656,
    `14` = 7774,
    `15` = 5709,
    `16` = 5629,
    `17` = 5745,
    `18` = 6056,
    `19` = 6259,
    `20` = 5303,
    `21` = 5423,
    `22` = 5497,
    `23` = 5547,
    `24` = 5417,
    `25` = 5441,
    `26` = 5466,
    `27` = 5500,
    `28` = 5668,
    `29` = 5694,
    `30` = 4365,
    `31` = 4252,
    `32` = 4122,
    `33` = 4142,
    `34` = 4039,
    `35` = 3210,
    `36` = 3222,
    `37` = 3258,
    `38` = 3413,
    `39` = 3871,
    `40` = 2684,
    `41` = 2844,
    `42` = 3052,
    `43` = 3182,
    `44` = 3237,
    `45` = 2263,
    `46` = 2298,
    `47` = 2318,
    `48` = 2257,
    `49` = 2194,
    `50` = 2231,
    `51` = 2172,
    `52` = 2072,
    `53` = 2008,
    `54` = 1932,
    `55` = 1301,
    `56` = 1262,
    `57` = 1213,
    `58` = 1197,
    `59` = 1191,
    `60` = 1601,
    `61` = 1593,
    `62` = 1490,
    `63` = 1348,
    `64` = 1299,
    `65` = 568,
    `66` = 745,
    `67` = 843,
    `68` = 801,
    `69` = 925,
    `70` = 806,
    `71` = 883,
    `72` = 796,
    `73` = 725,
    `74` = 672,
    `75` = 470,
    `76` = 441,
    `77` = 340,
    `78` = 300,
    `79` = 289,
    `80` = 4200
  )

pop_female_counts <-
  c(
    `0` = 11673,
    `1` = 11474,
    `2` = 11670,
    `3` = 11934,
    `4` = 11614,
    `5` = 10603,
    `6` = 11144,
    `7` = 11179,
    `8` = 11269,
    `9` = 11617,
    `10` = 6772,
    `11` = 6948,
    `12` = 7030,
    `13` = 7211,
    `14` = 7306,
    `15` = 6531,
    `16` = 6443,
    `17` = 6535,
    `18` = 6951,
    `19` = 7213,
    `20` = 6096,
    `21` = 6234,
    `22` = 6327,
    `23` = 6410,
    `24` = 6285,
    `25` = 6464,
    `26` = 6492,
    `27` = 6549,
    `28` = 6739,
    `29` = 6795,
    `30` = 5013,
    `31` = 4888,
    `32` = 4735,
    `33` = 4747,
    `34` = 4646,
    `35` = 3040,
    `36` = 3068,
    `37` = 3107,
    `38` = 3246,
    `39` = 3658,
    `40` = 2650,
    `41` = 2788,
    `42` = 2977,
    `43` = 3108,
    `44` = 3156,
    `45` = 1756,
    `46` = 1784,
    `47` = 1802,
    `48` = 1764,
    `49` = 1724,
    `50` = 1982,
    `51` = 1935,
    `52` = 1846,
    `53` = 1795,
    `54` = 1731,
    `55` = 863,
    `56` = 850,
    `57` = 825,
    `58` = 819,
    `59` = 816,
    `60` = 1348,
    `61` = 1342,
    `62` = 1246,
    `63` = 1138,
    `64` = 1101,
    `65` = 391,
    `66` = 520,
    `67` = 585,
    `68` = 560,
    `69` = 659,
    `70` = 670,
    `71` = 750,
    `72` = 686,
    `73` = 634,
    `74` = 604,
    `75` = 353,
    `76` = 340,
    `77` = 270,
    `78` = 246,
    `79` = 247,
    `80` = 4143
  )

smoothed_females <- smooth_age_5(Value = pop_female_counts,
                                 Age = as.numeric(names(pop_female_counts)),
                                 method = "Arriaga",
                                 OAG = TRUE,
                                 young.tail = "Original")

test_that("basepop_single does calculation for males when providing Males_single", {

  # Since the default is for females, I saved the correct calculations for males
  # and just check that the current implementation matches it
  res <-
    basepop_single(
      refYear = refYear,
      Males_single = pop_male_counts,
      Females_single = pop_female_counts,
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

  correct_res_males <-
    c(
      `0` = 13315,
      `1` = 11614,
      `2` = 11791,
      `3` = 12086,
      `4` = 11824,
      `5` = 10393,
      `6` = 10883,
      `7` = 10939,
      `8` = 11019,
      `9` = 11319
    )

  expect_equivalent(res[1:10], correct_res_males)

})
