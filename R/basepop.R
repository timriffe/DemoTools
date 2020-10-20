#' BPA and BPE methods for adjusting age groups under 10
#' @description Adjust population counts for the age groups 0, 1-4 and 5-9
#' @details
#'
#' \code{basepop_five} and \code{basepop_single} can estimate both the BPA and
#' BPE methods. If the user specifies \code{SmoothedFemales}, both \code{basepop_*}
#' functions will return the BPA method. If \code{SmoothedFemales} is left empty,
#' both \code{basepop_*} functions will adjust using the BPE method.
#'
#' For \code{basepop_five}, adjusting the female population counts is the
#' default. For this, all arguments except \code{Males_five} and
#' \code{nLxMale} are needed. However, for adjusting the male population
#' counts, the user needs to specify all arguments as well as \code{Males_five},
#' \code{nLxMale} and \code{nLxFemale} and set \code{female = FALSE}.
#'
#' For \code{basepop_single}, the same procedure applies. The only difference
#' is that the vector \code{Males_five} is named \code{Males_single} and accepts
#' a vector of single ages rather than five year abridged age groups. Similarly,
#' the vector for females is \code{Females_single} rather than
#' \code{Females_five} and accepts single age groups.
#'
#' Currently, \code{basepop_five} works only with five year abridged age groups
#' while \code{basepop_single} works only with single year age groups.
#'
#' @section BPA:
#'
#' Description:
#'
#' The method estimates a smoothed population ages 10 and over and adjusts
#' the population under age 10 using the smoothed population and estimates
#' of fertility and mortality.
#'
#' Based on the smoothed female population counts, it rejuvenates the female
#' "reported" population 20 to 59 years of age for the two 5 year periods prior
#' to the census date to represent the female population in reproductive ages
#' 5 and 10 years earlier.  Based on the rejuvenated population and fertility
#' and mortality levels, the method then estimates the male and female births
#' during the two 5 year periods prior to the census date. Next, it projects
#' the two 5-year birth cohorts to the census date. The projected figures
#' represent the adjusted population ages 0 to 4 years and 5 to 9 years
#' at the census date.
#'
#' Advantages:
#'
#' (1) The method adjusts under-10 population to be consistent with fertility
#' and mortality levels and adjusted adult female population.
#'
#' Limitations:
#'
#' (1) BPA assumes a linear change in fertility and mortality during the decade
#' prior to the reference year.
#'
#' (2) The procedure ignores migration, which can lead to misleading results.
#' There are two issues. First, age groups 0-4 and 5-9 are subject to migration,
#' which would affect the comparability of estimated and reported populations
#' in the base year.  Second, the estimated size of age groups 0-4 and 5-9 are
#' calculated from numbers of women of reproductive age in the base year
#' rejuvenated to points in the past.  With migration, rejuvenated number of
#' women may exceed or be smaller than the number present, and giving
#' birth to children, in the decade prior to the base year.
#'
#' (3) BPA’s smoothing calculations may mask unusual, but real, variations
#' in population age group size.  Smoothing irregularities in age structure
#' not attributable to age misreporting will distort estimated births and
#' survived children in the base year.
#'
#' Assumptions:
#'
#' (1) No significant international migration took place within the
#' reference periods for the population, mortality, and fertility input.
#'
#' (2) The data input as the "reported" population is not affected by
#' underenumeration of persons in certain ages, nor by age misreporting.
#'
#' @section BPE:
#'
#' Description:
#'
#' The method adjusts the population under age 10 using the reported population
#' ages 10 and above and estimates of fertility and mortality.
#'
#' The method rejuvenates the reported female population 20 to 59 years of age
#' for the two 5 year periods prior to the census date to represent the female
#' population in reproductive ages 5 and 10 years earlier. Based on the
#' rejuvenated population and fertility and mortality levels, the method then
#' estimates the male and female births during the two 5 year periods prior to
#' the census date.  Next, it projects the two 5-year birth cohorts to the
#' census date.  The projected figures represent the adjusted population ages
#' 0 to 4 years and 5 to 9 years at the census date.
#'
#' Advantages:
#'
#' (1) The method adjusts the under-10 population to be consistent with
#' fertility and mortality levels and adult female population.
#'
#' Limitations:
#'
#' (1) BPE assumes a linear change in fertility and mortality during the decade
#' prior to the reference year.
#'
#' (2) The procedure ignores migration, which can lead to misleading results.
#' There are two issues.  First, age groups 0-4 and 5-9 are subject to
#' migration, which would affect the comparability of estimated and reported
#' populations in the base year.  Second, the estimated size of age groups
#' 0-4 and 5-9 are calculated from numbers of women of reproductive age in
#' the base year rejuvenated to points in the past.  With migration, rejuvenated
#' number of women may exceed or be smaller than the number present, and
#' giving birth to children, in the decade prior to the base year.
#'
#' (3) The method does not adjust for possible underenumeration and age
#' misreporting errors in the over-10 “reported” population. If the
#' reported population is subject to age-misreporting or age-sex-specific
#' underenumeration, the over-10 population should be smoothed or otherwise
#' corrected prior to use.
#'
#' Assumptions:
#'
#' (1) No significant international migration took place within the reference
#' periods for the population, mortality, and fertility input.
#'
#' (2) The data input as the “reported” population is not affected by
#' underenumeration of persons in certain ages, nor by age misreporting.
#'
#' @return When using \code{basepop_five}, a numeric vector of adjusted counts
#' for the age groups 0, 1-4 and 5-9 for either \code{Females} or \code{Males},
#' depending on whether the user specified \code{TRUE} or \code{FALSE} in
#' \code{female}. The remaining age counts are left as is. When
#' \code{basepop_single} is used, the return value issues a numeric vector with
#' **single year age groups** where the counts between 0 and 10 are adjusted.
#'
#' @param refYear An integer. The decimal reference year for which the reported population pertain
#' @param Females_five A named numeric vector. The names of the vector should reflect the age groups. See examples. Reported population by abridged 5-year lower bound age groups for females pertaining to \code{refYear}. This means that age 0 must be kept as a single age, ages 1-4 must grouped together as abridged age 1, and thereafter 5-year age groups.
#' @param SRB A numeric. Sex ratio at birth (males / females)
#' @param nLxFemale A numeric matrix. The female nLx function of two life tables with ages in the rows and time in columns. The earlier date should be at least 7.5 years before the reference date of the "reported" population. The later date should be no earlier than one-half year before the reference date of the "reported" population
#' @param nLxDatesIn Vector of dates. The dates which pertain to the columns in \code{nLxFemale} and \code{nLxMale}. For details on the format of the dates, see the details of \code{\link{interp}}
#' @param asfrMat A numeric matrix. An age-period matrix of age specific fertility rates with age in rows, time in columns. 
#' @param AsfrDatesIn A vector of dates. The dates which pertain to the columns in \code{asfrMat}. For details on the format of the dates, see the details of \code{\link{interp}}
#' @param ... Arguments passed to \code{\link{interp}}. In particular, users might be interested in changing the interpolation method. By default, it's linearly interpolated
#' @param female A logical \code{TRUE} or \code{FALSE}. See the details section for the mandatory arguments for adjusting the male counts.
#' @param SmoothedFemales A numeric vector. This should be the result of passing the \code{Females} arguments through \code{smooth_age_5}. It is assumed that this argument is a numeric vector with the lower bound age group definitions as names. By default, this is what \code{smooth_age_5} returns. See the example section for an example.
#' @param Males_five A named numeric vector. The names of the vector should reflect the age groups. See examples. Reported population by abridged 5-year lower bound age groups for males pertaining to \code{refYear}. This means that age 0 must be kept as a single age, ages 1-4 must grouped together as abridged age 1, and thereafter 5-year age groups. This argument is only used when \code{female} is set to \code{FALSE}.
#' @param nLxMale A numeric matrix. The male nLx function of two life tables with ages in the rows and time in columns. The male nLx functions should be for age groups 0, 1 to 4, and 5 to 9. The dates which are represented in the columns are assumed to be the same as \code{nLxDatesIn}. This argument is only used when \code{female} is set to \code{FALSE}.
#'
#' @export
#' @examples
#'
#'  \dontrun{
#'  # (1) RefYear
#'  refYear <- 1986.21
#'
#'  # (2) Reported population by 5-year age groups and sex in the base year
#'  # (Include unknowns).
#'
#'  pop_male_counts <- c(`0` = 11684,
#'                       `1` = 46738,
#'                       `5` = 55639,
#'                       `10` = 37514,
#'                       `15` = 29398,
#'                       `20` = 27187,
#'                       `25` = 27770,
#'                       `30` = 20920,
#'                       `35` = 16973,
#'                       `40` = 14999,
#'                       `45` = 11330,
#'                       `50` = 10415,
#'                       `55` = 6164,
#'                       `60` = 7330,
#'                       `65` = 3882,
#'                       `70` = 3882,
#'                       `75` = 1840,
#'                       `80` = 4200
#'                       )
#'
#'  pop_female_counts <- c(`0` = 11673,
#'                         `1` = 46693,
#'                         `5` = 55812,
#'                         `10` = 35268,
#'                         `15` = 33672,
#'                         `20` = 31352,
#'                         `25` = 33038,
#'                         `30` = 24029,
#'                         `35` = 16120,
#'                         `40` = 14679,
#'                         `45` = 8831,
#'                         `50` = 9289,
#'                         `55` = 4172,
#'                         `60` = 6174,
#'                         `65` = 2715,
#'                         `70` = 3344,
#'                         `75` = 1455,
#'                         `80` = 4143)
#'
#'  # (4) Sex ratio at birth (m/f)
#'  sex_ratio <- 1.0300
#'
#'  # (6) The male and female nLx functions for ages under 1 year, 1 to 4 years, and 5 to 9
#'  # years, pertaining to an earlier and later date
#'  nLxDatesIn <- c(1977.31, 1986.50)
#'
#'  nLxMale <- matrix(c(87732,
#'                      304435,
#'                      361064,
#'                      88451,
#'                      310605,
#'                      370362
#'                      ),
#'                    nrow = 3, ncol = 2)
#'
#'  nLxFemale <- matrix(c(89842,
#'                        314521,
#'                        372681,
#'                        353053,
#'                        340650,
#'                        326588,
#'                        311481,
#'                        295396,
#'                        278646,
#'                        261260,
#'                        241395,
#'                        217419,
#'                        90478,
#'                        320755,
#'                        382531,
#'                        364776, 
#'                        353538, 
#'                        340687, 
#'                        326701, 
#'                        311573, 
#'                        295501, 
#'                        278494, 
#'                        258748,
#'                        234587),
#'                      nrow = 12,
#'                      ncol = 2)
#'
#'  # (7) A set of age-specific fertility rates pertaining to an earlier and later
#'  # date
#'
#'  asfrmat <- matrix(c("15-19" = 0.2000,
#'                      "20-24" = 0.3000,
#'                      "25-29" = 0.3000, 
#'                      "30-34" = 0.2500, 
#'                      "35-39" = 0.2000, 
#'                      "40-44" = 0.1500, 
#'                      "45-49" = 0.0500,
#'                      "15-19" = 0.1500,
#'                      "20-24" = 0.2000,
#'                      "25-29" = 0.2750, 
#'                      "30-34" = 0.2250, 
#'                      "35-39" = 0.1750, 
#'                      "40-44" = 0.1250, 
#'                      "45-49" = 0.0500),
#'                    nrow = 7,
#'                    ncol = 2)
#'
#'  AsfrDatesIn <- c(1977.81, 1985.71)
#'
#'  smoothed_females <- smooth_age_5(Value = pop_female_counts,
#'                                   Age = as.numeric(names(pop_female_counts)),
#'                                   method = "Arriaga",
#'                                   OAG = TRUE,
#'                                   young.tail = "Original")
#'
#'  ## This is the only number that messes up the whole calculation.
#'  ## smooth_age_5 returns the same result as the PASS excel sheet
#'  ## except for the age groups 10-15 and 15-19. Here we only use
#'  ## age group 15-19. If we plug in manually the correct value,
#'  ## we get all results match exactly, otherwise there are
#'  ## some differences.
#'  smoothed_females[4] <- 34721
#'
#'  # For adjusting using BPA for males, we need to specify
#'  # female = FALSE with Males and nLxMale.
#'  bpa_male <-
#'    basepop_five(
#'      refYear = refYear,
#'      Males_five = pop_male_counts,
#'      Females_five = pop_female_counts,
#'      SmoothedFemales = smoothed_females,
#'      SRB = sex_ratio,
#'      nLxFemale = nLxFemale,
#'      nLxMale = nLxMale,
#'      nLxDatesIn = nLxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn,
#'      female = FALSE
#'    )
#'
#'  # See adjustments?
#'  pop_male_counts[1:3]
#'  bpa_male[1:3]
#'
#'  # Adjusting the BPA for females requires less arguments
#'
#'  bpa_female <-
#'    basepop_five(
#'      refYear = refYear,
#'      Females_five = pop_female_counts,
#'      SmoothedFemales = smoothed_females,
#'      SRB = sex_ratio,
#'      nLxFemale = nLxFemale,
#'      nLxDatesIn = nLxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn
#'    )
#'
#'  pop_female_counts[1:3]
#'  bpa_female[1:3]
#'
#'  # For adjustment using BPE, we use exactly the same definitions as above
#'  # but remove SmoothedFemales.
#'  bpe_male <-
#'    basepop_five(
#'      refYear = refYear,
#'      Males_five = pop_male_counts,
#'      Females_five = pop_female_counts,
#'      SRB = sex_ratio,
#'      nLxFemale = nLxFemale,
#'      nLxMale = nLxMale,
#'      nLxDatesIn = nLxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn,
#'      female = FALSE
#'      )
#'
#'  # See adjustments?
#'  pop_male_counts[1:3]
#'  bpa_male[1:3]
#'  bpe_male[1:3]
#'
#'  # Adjusting the BPA for females requires less arguments
#'
#'  bpe_female <-
#'    basepop_five(
#'      refYear = refYear,
#'      Females_five = pop_female_counts,
#'      SRB = sex_ratio,
#'      nLxFemale = nLxFemale,
#'      nLxDatesIn = nLxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn
#'    )
#'
#'  pop_female_counts[1:3]
#'  bpa_female[1:3]
#'  bpe_female[1:3]
#'
#' # basepop_single for single ages
#' # Single ages for males and females
#'
#' pop_male_counts <-
#'   c(
#'     `0` = 11684,
#'     `1` = 11473,
#'     `2` = 11647,
#'     `3` = 11939,
#'     `4` = 11680,
#'     `5` = 10600,
#'     `6` = 11100,
#'     `7` = 11157,
#'     `8` = 11238,
#'     `9` = 11544,
#'     `10` = 7216,
#'     `11` = 7407,
#'     `12` = 7461,
#'     `13` = 7656,
#'     `14` = 7774,
#'     `15` = 5709,
#'     `16` = 5629,
#'     `17` = 5745,
#'     `18` = 6056,
#'     `19` = 6259,
#'     `20` = 5303,
#'     `21` = 5423,
#'     `22` = 5497,
#'     `23` = 5547,
#'     `24` = 5417,
#'     `25` = 5441,
#'     `26` = 5466,
#'     `27` = 5500,
#'     `28` = 5668,
#'     `29` = 5694,
#'     `30` = 4365,
#'     `31` = 4252,
#'     `32` = 4122,
#'     `33` = 4142,
#'     `34` = 4039,
#'     `35` = 3210,
#'     `36` = 3222,
#'     `37` = 3258,
#'     `38` = 3413,
#'     `39` = 3871,
#'     `40` = 2684,
#'     `41` = 2844,
#'     `42` = 3052,
#'     `43` = 3182,
#'     `44` = 3237,
#'     `45` = 2263,
#'     `46` = 2298,
#'     `47` = 2318,
#'     `48` = 2257,
#'     `49` = 2194,
#'     `50` = 2231,
#'     `51` = 2172,
#'     `52` = 2072,
#'     `53` = 2008,
#'     `54` = 1932,
#'     `55` = 1301,
#'     `56` = 1262,
#'     `57` = 1213,
#'     `58` = 1197,
#'     `59` = 1191,
#'     `60` = 1601,
#'     `61` = 1593,
#'     `62` = 1490,
#'     `63` = 1348,
#'     `64` = 1299,
#'     `65` = 568,
#'     `66` = 745,
#'     `67` = 843,
#'     `68` = 801,
#'     `69` = 925,
#'     `70` = 806,
#'     `71` = 883,
#'     `72` = 796,
#'     `73` = 725,
#'     `74` = 672,
#'     `75` = 470,
#'     `76` = 441,
#'     `77` = 340,
#'     `78` = 300,
#'     `79` = 289,
#'     `80` = 4200
#'   )
#' pop_female_counts <-
#'   c(
#'     `0` = 11673,
#'     `1` = 11474,
#'     `2` = 11670,
#'     `3` = 11934,
#'     `4` = 11614,
#'     `5` = 10603,
#'     `6` = 11144,
#'     `7` = 11179,
#'     `8` = 11269,
#'     `9` = 11617,
#'     `10` = 6772,
#'     `11` = 6948,
#'     `12` = 7030,
#'     `13` = 7211,
#'     `14` = 7306,
#'     `15` = 6531,
#'     `16` = 6443,
#'     `17` = 6535,
#'     `18` = 6951,
#'     `19` = 7213,
#'     `20` = 6096,
#'     `21` = 6234,
#'     `22` = 6327,
#'     `23` = 6410,
#'     `24` = 6285,
#'     `25` = 6464,
#'     `26` = 6492,
#'     `27` = 6549,
#'     `28` = 6739,
#'     `29` = 6795,
#'     `30` = 5013,
#'     `31` = 4888,
#'     `32` = 4735,
#'     `33` = 4747,
#'     `34` = 4646,
#'     `35` = 3040,
#'     `36` = 3068,
#'     `37` = 3107,
#'     `38` = 3246,
#'     `39` = 3658,
#'     `40` = 2650,
#'     `41` = 2788,
#'     `42` = 2977,
#'     `43` = 3108,
#'     `44` = 3156,
#'     `45` = 1756,
#'     `46` = 1784,
#'     `47` = 1802,
#'     `48` = 1764,
#'     `49` = 1724,
#'     `50` = 1982,
#'     `51` = 1935,
#'     `52` = 1846,
#'     `53` = 1795,
#'     `54` = 1731,
#'     `55` = 863,
#'     `56` = 850,
#'     `57` = 825,
#'     `58` = 819,
#'     `59` = 816,
#'     `60` = 1348,
#'     `61` = 1342,
#'     `62` = 1246,
#'     `63` = 1138,
#'     `64` = 1101,
#'     `65` = 391,
#'     `66` = 520,
#'     `67` = 585,
#'     `68` = 560,
#'     `69` = 659,
#'     `70` = 670,
#'     `71` = 750,
#'     `72` = 686,
#'     `73` = 634,
#'     `74` = 604,
#'     `75` = 353,
#'     `76` = 340,
#'     `77` = 270,
#'     `78` = 246,
#'     `79` = 247,
#'     `80` = 4143
#'   )
#'
#'  smoothed_females <- smooth_age_5(Value = pop_female_counts,
#'                                   Age = as.numeric(names(pop_female_counts)),
#'                                   method = "Arriaga",
#'                                   OAG = TRUE,
#'                                   young.tail = "Original")
#'
#'
#'  # For adjusting using BPA for males, we need to specify
#'  # female = FALSE with Males and nLxMale.
#'  bpa_male <-
#'    basepop_single(
#'      refYear = refYear,
#'      Males_single = pop_male_counts,
#'      Females_single = pop_female_counts,
#'      SmoothedFemales = smoothed_females,
#'      SRB = sex_ratio,
#'      nLxFemale = nLxFemale,
#'      nLxMale = nLxMale,
#'      nLxDatesIn = nLxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn,
#'      female = FALSE
#'      )
#'
#'  # See adjustments?
#'  pop_male_counts[1:10]
#'  bpa_male[1:10]
#'
#'  # Adjusting the BPA for females requires less arguments
#'  bpa_female <-
#'    basepop_single(
#'      refYear = refYear,
#'      Females_single = pop_female_counts,
#'      SmoothedFemales = smoothed_females,
#'      SRB = sex_ratio,
#'      nLxFemale = nLxFemale,
#'      nLxDatesIn = nLxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn
#'    )
#'
#'  pop_female_counts[1:10]
#'  bpa_female[1:10]
#'
#'  # For adjustment using BPE, we use exactly the same definitions as above
#'  # but remove SmoothedFemales.
#'  bpe_male <-
#'    basepop_single(
#'      refYear = refYear,
#'      Males_single = pop_male_counts,
#'      Females_single = pop_female_counts,
#'      SRB = sex_ratio,
#'      nLxFemale = nLxFemale,
#'      nLxMale = nLxMale,
#'      nLxDatesIn = nLxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn,
#'      female = FALSE
#'      )
#'
#'  # See adjustments?
#'  pop_male_counts[1:10]
#'  bpa_male[1:10]
#'  bpe_male[1:10]
#'
#'  # Adjusting the BPA for females requires less arguments
#'  bpe_female <-
#'    basepop_single(
#'      refYear = refYear,
#'      Females_single = pop_female_counts,
#'      SRB = sex_ratio,
#'      nLxFemale = nLxFemale,
#'      nLxDatesIn = nLxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn
#'    )
#'
#'  pop_female_counts[1:10]
#'  bpa_female[1:10]
#'  bpe_female[1:10]
#'
#'  }
#'
#' @references
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#'
basepop_five <- function(refYear,
                         Females_five,
                         SRB,
                         nLxFemale,
                         nLxDatesIn,
                         asfrMat,
                         AsfrDatesIn,
                         ...,
                         female = TRUE,
                         SmoothedFemales = NULL,
                         Males_five = NULL,
                         nLxMale = NULL,

  if (missing(nLxDatesIn)) {

  ## Check all arguments
  AllArgs <- as.list(environment())
  ArgsCheck(AllArgs)

  # We calculate basepop for males **ONLY** if the nLxMale
  # argument is not null
  nLxGender <- if (isFALSE(female)) nLxMale else nLxFemale
  male <- if (isFALSE(female)) TRUE else FALSE
  DatesOut <- sapply(c(0.5, 2.5, 7.5), function(x) refYear - x)

  # Interpolate the gender specific nLx to the requested
  # dates out
  nLx <- interp(
    nLxGender,
    datesIn = nLxDatesIn,
    datesOut = DatesOut,
    ...
  )
  # Turn the columns from the matrix into a list
  nLx <- lapply(as.data.frame(nLx), identity)

  # Interpolate only the female nLx to the requested dates.
  # This is independent of the nLx above, since to apply
  # basepop to either gender you need to have the female nLx
  nLxFemale <- interp(
    nLxFemale,
    datesIn = nLxDatesIn,
    datesOut = DatesOut,
    ...
  )

  # Turn the columns from the matrix into a list
  nLxFemale <- lapply(as.data.frame(nLxFemale), identity)

  # Interpolate the asfr to the requested dates.
  # This is gender agnostic.
  Asfr <- interp(
    asfrMat,
    datesIn = AsfrDatesIn,
    datesOut = DatesOut,
    ...
  )
  # Turn the columns from the matrix into a list
  Asfr <- lapply(as.data.frame(Asfr), identity)

  EarliestDate <- which.max(names(nLxFemale))
  OlderNLxFemale <- nLxFemale[-EarliestDate]
  # Reorder the years such that the earlier year is first
  OlderNLxFemale <- OlderNLxFemale[sort(names(OlderNLxFemale), decreasing = TRUE)]

  # If SmoothedFemales is specified, we assume that it is returned
  # by smooth_age_5, so it should be named with the age groups.
  # This is checked in ArgList
  FemalePops <- if (!is.null(SmoothedFemales)) SmoothedFemales else Females_five

  # We use the smoothed vector names to only get certain age groups
  SmoothedFMiddleages <- FemalePops[as.character(seq(15, 55, by = 5))]

  # Currently, we assume that for calculating the estimated
  # female population for the DatesOut, the first OlderNLxFemale
  # date is multiplied by the SmoothedFMiddleages. However,
  # all dates after that one are mutltiplied by the previous
  # OlderNLxFemale in the order.
  for (i in seq_along(OlderNLxFemale)) {
    .x <- OlderNLxFemale[[i]]
    VecMult <- if (i == 1) SmoothedFMiddleages else OlderNLxFemale[[i - 1]]
    IndLength <- if (i == 1) 0 else 1

    NLxDiv <- .x[4:(length(.x) - IndLength)]
    iter <- c(1, rep(seq(2, length(NLxDiv) - 1), each = 2), length(NLxDiv))
    NLxSeq <- lapply(seq(1, length(iter), by = 2), function(i) NLxDiv[iter[i:(i+1)]])

    # Here VecMult is either SmoothedFMiddleAges
    OlderNLxFemale[[i]] <- mapply(function(.y, .z) .y * .z[1] / .z[2],
                                  VecMult[-1],
                                  NLxSeq
                                  )
  }

  # We always assume that the FirstDate will be calculate with the first date
  # from the OlderNLxFemale, regardless of the number of dates.
  FirstDate <- list(0.2 * OlderNLxFemale[[1]] + 0.8 * SmoothedFMiddleages[-length(SmoothedFMiddleages)])
  names(FirstDate) <- max(names(nLxFemale))
  OlderNLxFemale <- c(
    FirstDate,
    OlderNLxFemale
  )

  # Currently, we assume that for estimating the population
  # for DatesOut, starting from the third year, we carry forward
  # the mean between the current year the previous one instead
  # of the SmoothedFMiddleages as fixed first year
  EstPop <- vector("list", length(OlderNLxFemale))
  names(EstPop) <- names(OlderNLxFemale)
  for (i in seq_along(OlderNLxFemale)) {
    .x <- OlderNLxFemale[[i]]
    MeanVec <- if (i > 2) OlderNLxFemale[[i-1]] else SmoothedFMiddleages

    # Loop along the 15-45 age groups
    EstPop[[i]] <- vapply(seq_along(seq(15, 45, by = 5)),
                           function(i) mean(c(MeanVec[i], .x[i])),
                           FUN.VALUE = numeric(1))
  }

  # This contains the estimate male/female population (depending on whether
  # the user specified nLxMale) for the DatesOut.
  EstTot <- mapply(
    AdjustFemalePop,
    EstPop,
    Asfr,
    SRB,
    male,
    SIMPLIFY = FALSE
  )

  GenderCounts <- if (male) Males_five else Females_five

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


#' @rdname basepop_five
#' @aliases basepop_five
#' @param Females_single A named numeric vector. Reported population by 1-year age groups for \code{refYear} for females. The names of the vector should reflect the age groups. See examples. The method assumes that the last age group is open (for example, the population ends at '80+' and '100+')
#' @param Males_single A named numeric vector. Reported population by 1-year age groups format \code{refYear} for males. The names of the vector should reflect the age groups. See examples. The method assumes that the last age group is open (for example, the population ends at '80+' and '100+')
#'
#' @export
#'
basepop_single <- function(refYear,
                           Females_single,
                           SRB,
                           nLxFemale,
                           nLxDatesIn,
                           asfrMat,
                           AsfrDatesIn,
                           ...,
                           female = TRUE,
                           SmoothedFemales = NULL,
                           Males_single = NULL,
                           nLxMale = NULL) {

  stopifnot(
    !is.null(names(Females_single)),
    is_single(as.numeric(names(Females_single)))
  )

  Females_abridged <- single2abridged(Females_single)
  males_present <- !is.null(Males_single)

  if (males_present) {
    stopifnot(
      !is.null(names(Males_single)),
      is_single(as.numeric(names(Males_single)))
    )

    Males_abridged <- single2abridged(Males_single)
    gender_single <- Males_single
  }  else {
    Males_abridged <- Males_single
    gender_single <- Females_single
  }

  res <-
    basepop_five(
      refYear = refYear,
      Females_five = Females_abridged,
      SRB = SRB,
      nLxFemale = nLxFemale,
      nLxDatesIn = nLxDatesIn,
      asfrMat = asfrMat,
      AsfrDatesIn = AsfrDatesIn,
      ... = ...,
      female = female,
      SmoothedFemales = SmoothedFemales,
      Males_five = Males_abridged,
      nLxMale = nLxMale
    )

  # Since diff always returns a vector of length `length(x) - 1`,
  # the 1 in the end is to reflct the the open ages for 80+ or 100+
  AgeBins1 <- c(diff(as.integer(names(gender_single))), 1)
  AgeBins2 <- c(diff(as.integer(names(res))), 1)

  rescaled_res <-
    rescaleAgeGroups(
      Value1 = gender_single,
      AgeInt1 = AgeBins1,
      Value2 = res,
      AgeInt2 = AgeBins2,
      splitfun = graduate_uniform
    )

  round(rescaled_res, 0)
}

single2abridged <- function(Age) {
  age_names <- as.integer(seq_along(Age) - 1)
  abridged_index <- calcAgeAbr(age_names)
  tapply(Age, abridged_index, sum)
}

AdjustFemalePop <- function(x, asfr_pastyear, sex_ratio, male = TRUE) {
  births_age_female <- x * asfr_pastyear
  est_tot <- sum(births_age_female) * sex_ratio / (1 + sex_ratio)

  if (!male) est_tot <- sum(births_age_female) - est_tot

  est_tot
}

ArgsCheck <- function(ArgList) {
  with(ArgList, {
    stopifnot(
      is.numeric(Females_five),
      length(SRB) == 1,
      is.numeric(SRB),
      is.matrix(nLxFemale),
      is.matrix(asfrMat),
      is.logical(female),
      ncol(nLxFemale) == length(nLxDatesIn)
    )

    if (isFALSE(female)) {
      stopifnot(
        is.numeric(Males_five),
        is.matrix(nLxMale),
        ncol(nLxMale) == length(nLxDatesIn),
        !is.null(Males_five),
        is.logical(female),
        !is.null(nLxMale))
    }

    is.named <- function(x) !is.null(names(x))
    if (!is.null(SmoothedFemales)) stopifnot(is.named(SmoothedFemales))

  })
}

downloadnLx <- function(nLx, country, refYear, gender, nLxDatesIn) {
  if (!missing(nLx)) return(nLx)
  if (is.null(nLx)) return(nLx)

  tmp <-
    lapply(nLxDatesIn, function(x) {
      print(paste0("Downloading nLx data for ", country, ", year ", x, ", gender ", gender))
      res <- fertestr::FetchLifeTableWpp2019(country, x, gender)["Lx"]
      names(res) <- NULL
      as.matrix(res)
    })

  nLx <- do.call(cbind, tmp)
  nLx
}

downloadAsfr <- function(Asfrmat, country, refYear, AsfrDatesIn) {

  if (!missing(Asfrmat)) {
    return(Asfrmat)
  }

  tmp <-
    lapply(AsfrDatesIn, function(x) {
      print(paste0("Downloading Asfr data for ", country, ", year ", x))
      res <- fertestr::FetchFertilityWpp2019(country, x)["asfr"]
      names(res) <- NULL
      as.matrix(res)
    })

  Asfrmat <- do.call(cbind, tmp)
  Asfrmat
}

## bpa_male <-
##   basepop_five(
##     country = "Spain",
##     refYear = refYear,
##     Females_five = pop_female_counts,
##     SmoothedFemales = smoothed_females
##   )
