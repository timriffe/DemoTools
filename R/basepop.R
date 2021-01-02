# TODO: 
# -[x] alwys do males and females together
# -[x] include SRB in DemoToolsData package
# -[x] allow SRB as scalar, vector, or lookup using country, pulling 
#      data from DemoToolsData, in which case has to be a time series.



#' BPA and BPE methods for adjusting age groups under 10
#' @description Adjust population counts for the age groups 0 to 10
#'
#' @details
#'
#' \code{basepop_five} and \code{basepop_single} can estimate both the BPA and
#' BPE methods. If the user specifies \code{SmoothedFemales}, both
#' \code{basepop_*} functions will return the BPA method.
#' If \code{SmoothedFemales} is left empty, both \code{basepop_*} functions will
#' adjust using the BPE method.
#'
#' For \code{basepop_five}, adjusting the female population counts is the
#' default. For this, only the \code{country}, \code{refDate} and
#' \code{Females_five} are needed. All other arguments are downloaded
#' or set to sensible defaults. For adjusting the male population
#' counts, the user needs to specify the \code{Males_five} population
#' counts and set \code{female = FALSE}.
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
#' The BPE method is used by default. To adjust the counts using
#' the BPA method, the user needs to provide the \code{SmoothedFemales}
#' argument. This is the female population counts passed through
#' a smoothing function such as \code{smooth_age_5}. See the examples
#' section for some examples.
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
#' for the age groups 0, 1-4 and 5-9 for either \code{females} or \code{males},
#' depending on whether the user specified \code{TRUE} or \code{FALSE} in
#' \code{female}. The remaining age counts are left as is. When
#' \code{basepop_single} is used, the return value is a numeric vector with
#' **single year age groups** where the counts between 0 and 10 are adjusted.
#'
#' @param country The country name or location code from which to download the n
#' Lx and asfr data. See \code{fertestr::locs_avail()} for all country
#' names/codes.
#'
#' @param refDate The reference year for which the reported population pertain
#' (these are the population counts in \code{Females_five} and
#' \code{Males_five}).
#'
#' @param Females_five A named numeric vector with the population counts for
#' five-year abridged age groups for females in \code{refDate}. The names of the
#' vector should reflect the age groups. See the example section for some
#' examples.
#'
#' @param nLxFemale A numeric matrix. The female nLx function of two life tables
#' with ages in the rows and time in columns. The earlier date should be at least
#' 7.5 years before the reference date of the "reported" population. The later
#' date should be no earlier than one-half year before the reference date of
#' the "reported" population. If not provided, it's automatically downloaded if
#' \code{country}, \code{refDate} and the equivalent population counts
#' \code{*_five} are provided.
#'
#' @param nLxDatesIn A vector of numeric years (for example, 1986). The dates
#' which pertain to the columns in \code{nLxFemale} and \code{nLxMale}. If not
#' provided, the function automatically determines two dates which are 8 years
#' before \code{refDate} and 0.5 years after \code{refDate}.
#'
#' @param AsfrMat A numeric matrix. An age-period matrix of age specific
#' fertility rates with age in rows, time in columns. If not provided, the
#' function automatically downloads the ASFR matrix based on the dates in
#' \code{AsfrDatesIn}.
#'
#' @param AsfrDatesIn A vector of numeric years (for example, 1986). These are
#' the dates which pertain to the columns in \code{AsfrMat}. If not provided,
#' the function automatically determines two dates which are 8 years before
#' \code{refDate} and 0.5 before \code{refDate}.
#'
#' @param ... Arguments passed to \code{\link{interp}}. In particular, users
#' might be interested in changing the interpolation method for the \code{nLx*}
#' matrices and the \code{Asfr} matrix. By default, it's linearly interpolated.
#'
#' @param female A logical \code{TRUE} or \code{FALSE} for whether to perform
#' the adjustment for females or males. If \code{female = FALSE},
#' \code{Males_five} needs to be provided.
#'
#' @param SmoothedFemales A numeric vector. This should be the result of
#' passing the \code{Females_five} argument through \code{smooth_age_5}. It
#' is assumed that this argument is a numeric vector with the lower bound age
#' group definitions as names. By default, this is what \code{smooth_age_5}
#' returns. See the example section for some examples.
#'
#' @param Males_five A named numeric vector with the population counts for
#' five-year abridged age groups for males in \code{refDate}. The names of
#' the vector should reflect the age groups. See the example section for
#' some examples.
#'
#' @param nLxMale A numeric matrix. The male nLx function of two life tables
#' with ages in the rows and time in columns. The male nLx function should be
#' only for age groups 0, 1 to 4, and 5 to 9. The dates which are represented
#' in the columns are assumed to be the same as \code{nLxDatesIn}. This
#' argument is only used when \code{female} is set to \code{FALSE} and
#' \code{Males_five} is provided. If \code{Males_five} is provided and
#' \code{female} set to \code{FALSE}, the \code{nLx} for males is
#' automatically downloaded for the dates in \code{nLxDatesIn}.
#'
#' @param SRB A numeric. Sex ratio at birth (males / females). Default is set
#' to 1.05
#'
#' @param radix starting point to use in the adjustment of the three first age
#' groups. Default is NULL. If not provided, the function extracts the \code{lx}
#' for age group \code{0} of the first year defined in \code{AsfrDatesIn}.
#'
#' @param verbose when downloading new data, should the function print details
#' about the download at each step? Defaults to \code{TRUE}. We recommend the
#' user to set this to \code{TRUE} at all times because the function needs to
#' make decisions (such as picking the dates for the Asfr and nLx) that the user
#' should be aware of.
#'
#' @export
#' @examples
#'
#'  \dontrun{
#'
#' ################ BPE for females (five year age groups) #####################
#'
#' # Grab population counts for females
#' refDate <- 1986
#' country <- "Brazil"
#' pop_female_single <- fertestr::FetchPopWpp2019(country, 
#'                                                refDate, 
#'                                                ages = 0:100, 
#'                                                sex = "female")
#' pop_female_counts <- single2abridged(setNames(pop_female_single$pop, 
#'                                               pop_female_single$ages))
#' pop_male_single   <- fertestr::FetchPopWpp2019(country, 
#'                                                refDate, 
#'                                                ages = 0:100, 
#'                                                sex = "male")
#' pop_male_counts   <- single2abridged(setNames(pop_male_single$pop, 
#'                                               pop_male_single$ages))
#' Age <- names2age(pop_male_counts)
#' # Automatically downloads the nLx and ASFR data
#' bpe <- basepop_five(
#'   country = country,
#'   refDate = refDate,
#'   Females_five = pop_female_counts,
#'   Males_five = pop_male_counts,
#'   Age = Age
#' )
#'
#' # The counts for the first three age groups have been adjusted:
#' bpe$Females_adjusted[1:3]
#' pop_female_counts[1:3]
#'
#' bpe$Males_adjusted[1:3]
#' pop_male_counts[1:3]
#'
#'
#' ################ BPE for single ages) ############################
#'
#' # pop_female_single <- setNames(pop_female_single$pop, pop_female_single$ages)
#' # 
#' # # Automatically downloads the nLx and ASFR data
#' # bpe_female <- basepop_single(
#' #   country = country,
#' #   refDate = refDate,
#' #   Females_single = pop_female_single
#' # )
#' # 
#' # # The counts for the first 10 age groups have been adjusted:
#' # bpe_female[1:10]
#' # pop_female_single[1:10]
#' ################ BPA (five year age groups) #####################
#' # for BPA, smooth counts in advance
#' smoothed_females <- smooth_age_5(Value = pop_female_counts,
#'                                  Age = Age,
#'                                  method = "Arriaga",
#'                                  OAG = TRUE,
#'                                  young.tail = "Original")
#' smoothed_females <- c(pop_female_counts[1:2], smoothed_females[-1])
#' smoothed_males <- smooth_age_5(Value = pop_male_counts,
#'                                Age = Age,
#'                                method = "Arriaga",
#'                                OAG = TRUE,
#'                                young.tail = "Original")
#' smoothed_males <- c(smoothed_males[1:2], smoothed_males[-1])
#'
#' # Automatically downloads the nLx and ASFR data
#' bpa <- basepop_five(
#'   country = country,
#'   refDate = refDate,
#'   Females_five = smoothed_females,
#'   Males_five = smoothed_males,
#'   Age = Age
#' )
#'
#' # The counts for the first three age groups have been adjusted:
#' bpa$Females_adjusted[1:3]
#' smoothed_females[1:3]
#' pop_female_counts[1:3]
#'
#' bpa$Males_adjusted[1:3]
#' smoothed_males[1:3]
#' pop_male_counts[1:3]
#'
#' ################ PAS example ###############################
#'
#'  # (1) refDate
#'  refDate <- 1986.21
#'
#'  # (2) Reported population by 5-year age groups and sex in the base year
#'  # (Include unknowns).
#'
#'  pop_male_counts <- c(11684, 46738, 55639, 37514, 29398, 27187, 27770, 20920, 16973, 
#'                       14999, 11330, 10415, 6164, 7330, 3882, 3882, 1840, 4200)
#'  
#'  pop_female_counts <- c(11673, 46693, 55812, 35268, 33672, 31352, 33038, 24029, 16120,
#'                         14679, 8831, 9289, 4172, 6174, 2715, 3344, 1455, 4143)
#'  Age <- c(0,1, seq(5, 80, by = 5))
#'
#'  # (4) Sex ratio at birth (m/f)
#'  sex_ratio <- 1.0300
#'
#'  # (6) The male and female nLx functions for ages under 1 year, 1 to 4 years, and 5 to 9
#'  # years, pertaining to an earlier and later date
#'  nLxDatesIn <- c(1977.31, 1986.50)
#'
#'  nLxMale <- matrix(c(87732, 304435, 361064, 88451, 310605, 370362),
#'                    nrow = 3, ncol = 2)
#'
#'  nLxFemale <- matrix(c(89842, 314521, 372681, 353053, 340650, 326588, 
#'                        311481, 295396, 278646, 261260, 241395,217419,
#'                        90478, 320755, 382531, 364776, 353538, 340687, 
#'                        326701, 311573, 295501, 278494, 258748,234587),
#'                      nrow = 12,
#'                      ncol = 2)
#'
#'  # (7) A set of age-specific fertility rates pertaining to an earlier and later
#'  # date
#'
#'  asfrmat <- structure(
#'     c(0.2, 0.3, 0.3, 0.25, 0.2, 0.15, 0.05, 0.15, 0.2, 
#'       0.275, 0.225, 0.175, 0.125, 0.05), .Dim = c(7L, 2L), 
#'     .Dimnames = list(
#'       c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"), 
#'       c("1977.81", "1985.71")))
#'
#'  # for BPA, smooth counts in advance
#'  smoothed_females <- smooth_age_5(Value = pop_female_counts,
#'                                   Age = Age,
#'                                   method = "Arriaga",
#'                                   OAG = TRUE,
#'                                   young.tail = "Original")
#'  smoothed_females <- c(pop_female_counts[1:2], smoothed_females[-1])
#'  smoothed_males <- smooth_age_5(Value = pop_male_counts,
#'                                   Age = Age,
#'                                   method = "Arriaga",
#'                                   OAG = TRUE,
#'                                   young.tail = "Original")
#'  smoothed_males <- c(pop_male_counts[1:2], smoothed_males[-1])
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
#'  bpa <-
#'    basepop_five(
#'      refDate = refDate,
#'      Males_five = smoothed_males,
#'      Females_five = smoothed_females,
#'      Age = Age,
#'      SRB = sex_ratio,
#'      nLxFemale = nLxFemale,
#'      nLxMale = nLxMale,
#'      nLxDatesIn = nLxDatesIn,
#'      AsfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn,
#'      radix = 1e5
#'    )
#'
#'  # See adjustments?
#'  pop_male_counts[1:3]
#'  bpa$Male_adjusted[1:3]
#'
#'  pop_female_counts[1:3]
#'  bpa$Female_adjusted[1:3]
#'
#'  # For adjustment using BPE, we use exactly the same definitions as above
#'  # but use the original inputs
#'
#'  bpe <-
#'    basepop_five(
#'      refDate = refDate,
#'      Females_five = pop_female_counts,
#'      Males_five = pop_male_counts,
#'      SRB = sex_ratio,
#'      nLxFemale = nLxFemale,
#'      nLxDatesIn = nLxDatesIn,
#'      AsfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn
#'    )
#'
#'  pop_female_counts[1:3]
#'  bpe$Females_adjusted[1:3]
#'  
#' # basepop_single for single ages
#' # Single ages for males and females
#'
#' # pop_male_counts <-
#' #   c(11684, 11473, 11647, 11939, 11680, 10600, 11100, 11157, 11238, 
#' #     11544, 7216, 7407, 7461, 7656, 7774, 5709, 5629, 5745, 6056, 
#' #     6259, 5303, 5423, 5497, 5547, 5417, 5441, 5466, 5500, 5668, 5694, 
#' #     4365, 4252, 4122, 4142, 4039, 3210, 3222, 3258, 3413, 3871, 2684, 
#' #     2844, 3052, 3182, 3237, 2263, 2298, 2318, 2257, 2194, 2231, 2172, 
#' #     2072, 2008, 1932, 1301, 1262, 1213, 1197, 1191, 1601, 1593, 1490, 
#' #     1348, 1299, 568, 745, 843, 801, 925, 806, 883, 796, 725, 672, 
#' #     470, 441, 340, 300, 289, 4200)
#' # 
#' # pop_female_counts <-
#' #   c(11673, 11474, 11670, 11934, 11614, 10603, 11144, 11179, 11269, 
#' #     11617, 6772, 6948, 7030, 7211, 7306, 6531, 6443, 6535, 6951, 
#' #     7213, 6096, 6234, 6327, 6410, 6285, 6464, 6492, 6549, 6739, 6795, 
#' #     5013, 4888, 4735, 4747, 4646, 3040, 3068, 3107, 3246, 3658, 2650, 
#' #     2788, 2977, 3108, 3156, 1756, 1784, 1802, 1764, 1724, 1982, 1935, 
#' #     1846, 1795, 1731, 863, 850, 825, 819, 816, 1348, 1342, 1246, 
#' #     1138, 1101, 391, 520, 585, 560, 659, 670, 750, 686, 634, 604, 
#' #     353, 340, 270, 246, 247, 4143)
#' #  Age <- 0:80
#' # 
#' #  smoothed_females <- smooth_age_5(Value = pop_female_counts,
#' #                                   Age = Age,
#' #                                   method = "Arriaga",
#' #                                   OAG = TRUE,
#' #                                   young.tail = "Original")
#' #  smoothed_males <- smooth_age_5(Value = pop_male_counts,
#' #                                   Age = Age,
#' #                                   method = "Arriaga",
#' #                                   OAG = TRUE,
#' #                                   young.tail = "Original")
#'
#'  # For adjusting using BPA for males, we need to specify
#'  # female = FALSE with Males and nLxMale.
#'  
#'  # This needs work still
#'  # bpa_male <-
#'  #   basepop_single(
#'  #     refDate = refDate,
#'  #     Males_single = pop_male_counts,
#'  #     Females_single = pop_female_counts,
#'  #     SRB = sex_ratio,
#'  #     nLxFemale = nLxFemale,
#'  #     nLxMale = nLxMale,
#'  #     nLxDatesIn = nLxDatesIn,
#'  #     AsfrMat = asfrmat,
#'  #     AsfrDatesIn = AsfrDatesIn
#'  #     )
#'
#'  # See adjustments?
#'  # pop_male_counts[1:10]
#'  # bpa_male[1:10]
#'
#'  # Adjusting the BPA for females requires less arguments
#'  # bpa_female <-
#'  #   basepop_single(
#'  #     refDate = refDate,
#'  #     Females_single = pop_female_counts,
#'  #     SmoothedFemales = smoothed_females,
#'  #     SRB = sex_ratio,
#'  #     nLxFemale = nLxFemale,
#'  #     nLxDatesIn = nLxDatesIn,
#'  #     AsfrMat = asfrmat,
#'  #     AsfrDatesIn = AsfrDatesIn
#'  #   )
#'
#'  # pop_female_counts[1:10]
#'  # bpa_female[1:10]
#'  # 
#'  # # For adjustment using BPE, we use exactly the same definitions as above
#'  # # but remove SmoothedFemales.
#'  # bpe_male <-
#'  #   basepop_single(
#'  #     refDate = refDate,
#'  #     Males_single = pop_male_counts,
#'  #     Females_single = pop_female_counts,
#'  #     SRB = sex_ratio,
#'  #     nLxFemale = nLxFemale,
#'  #     nLxMale = nLxMale,
#'  #     nLxDatesIn = nLxDatesIn,
#'  #     AsfrMat = asfrmat,
#'  #     AsfrDatesIn = AsfrDatesIn,
#'  #     female = FALSE
#'  #     )
#'
#'  # See adjustments?
#'  # pop_male_counts[1:10]
#'  # bpa_male[1:10]
#'  # bpe_male[1:10]
#'
#'  # Adjusting the BPA for females requires less arguments
#'  # bpe_female <-
#'  #   basepop_single(
#'  #     refDate = refDate,
#'  #     Females_single = pop_female_counts,
#'  #     SRB = sex_ratio,
#'  #     nLxFemale = nLxFemale,
#'  #     nLxDatesIn = nLxDatesIn,
#'  #     AsfrMat = asfrmat,
#'  #     AsfrDatesIn = AsfrDatesIn
#'  #   )
#'  # 
#'  # pop_female_counts[1:10]
#'  # bpa_female[1:10]
#'  # bpe_female[1:10]
#'
#'  }
#'
#' @references
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#'
basepop_five <- function(country = NULL,
                         refDate,
                         Age = NULL,
                         Females_five,
                         Males_five = NULL,
                         nLxFemale = NULL,
                         nLxMale = NULL,
                         nLxDatesIn = NULL,
                         AsfrMat = NULL,
                         AsfrDatesIn = NULL,
                         ...,
                         SRB = NULL,
                         radix = NULL,
                         verbose = TRUE) {

  options(basepop_verbose = verbose)
  on.exit(options(basepop_verbose = NULL))

  # Ensure census date is numeric.
  # "YYYY-MM-DD" input is acceptable
  refDate <- dec.date(refDate)
  
  if (!is.null(Age)){
    stopifnot(is_abridged(Age))
  } else {
    if (!is.null(names(Females_five))){
      Age <- names2age(Females_five)
    } else {
      if (verbose) {
        cat("Assuming age groups are in standard abridged intervals")
      }
      # last resort = assume in abrided ages!
      Age <- inferAgeIntAbr(Females_five)
    }
  }
  
  if (is.null(nLxDatesIn)) {
    nLxDatesIn <- c(abs(8 - refDate), refDate + 0.5)
    if (verbose) {
      cat(paste0("Assuming the two prior dates for the nLx matrix to be: ", paste0(nLxDatesIn, collapse = ", ")), sep = "\n")
    }
  }

  if (is.null(AsfrDatesIn)) {
    AsfrDatesIn <- abs(c(8, 0.5) - refDate)
    if (verbose) {
      cat(paste0("Assuming the two prior dates for the Asfr matrix to be: ", paste0(AsfrDatesIn, collapse = ", ")), sep = "\n")
    }
  }

  ## obtain nLx for males and females
  ## If these arguments have been specified, they return
  ## the same thing and don't download the data
  nLxFemale <-
    downloadnLx(
      nLx = nLxFemale,
      country = country,
      gender = "female",
      nLxDatesIn = nLxDatesIn
    )
  
  nLxMale <-
    downloadnLx(
      nLx = nLxMale,
      country = country,
      gender = "male",
      nLxDatesIn = nLxDatesIn
    )
  
  if (is.null(radix)) {
    # TR: not perfect, but it's a better guess. It would seem the radix
    # being pulled before was always 1, whereas the nLx columns was based on 100000
    radix <- lt_infer_radix_from_1L0(nLxMale[1,1])
    if (verbose) {
      cat(paste0("Setting radix to value of lx: ", radix, ". Can be overwritten with the `radix` argument"), sep = "\n")
    }
  }
  
  AsfrMat <-
    downloadAsfr(
      Asfrmat = AsfrMat,
      country = country,
      AsfrDatesIn = AsfrDatesIn
    )
  # get a vector of 3 SRB estimates matching the DatesOut dates.
  # if SRB was given as a vector of length 3 then we take it as-is
  # if only one value was given (or a vector of length not equal to 3), 
  # we repeat it 3 times and take the first 3 elements.
  # if it's NULL and we have the country in the DB then we look it up.
  # if it's NULL and we don't have the country then we assume 1.05,
  # because tradition.
  
  # TR saw no need for sapply()
  # DatesOut  <- sapply(c(0.5, 2.5, 7.5), function(x) refDate - x)
  DatesOut <- refDate - c(0.5, 2.5, 7.5)
  
  SRB <- downloadSRB(SRB, 
                     country, 
                     DatesOut)
  
  ## Check all arguments
  AllArgs <- as.list(environment())
  ArgsCheck(AllArgs)


  # Interpolate the gender specific nLx to the requested
  # dates out
  nLxf <- interp(
    nLxFemale,
    datesIn = nLxDatesIn,
    datesOut = DatesOut,
    ...
  )
  nLxm <- interp(
    nLxMale,
    datesIn = nLxDatesIn,
    datesOut = DatesOut,
   ...
  )

  # Interpolate the asfr to the requested dates.
  # This is gender agnostic.
  Asfr <- interp(
    AsfrMat,
    datesIn = AsfrDatesIn,
    datesOut = DatesOut,
    ...
  )

  # TR: Follows spreadsheet logic, can still be more elegant.
  # sometimes character indexing, sometimes position, but still
  ages_15_55         <- as.character(seq(15,55,by=5))
  ages_20_55         <- ages_15_55[-1]
  ages_15_50         <- ages_15_55[-9]
  ages_20_50         <- ages_15_55[-c(1,9)]
  ages_15_45         <- ages_15_55[-c(8,9)]
  ages_20_45         <- ages_15_55[-c(1,8,9)]
  ages_15_40         <- ages_15_55[-c(7,8,9)]
  
  FMiddleages        <- Females_five[ages_15_55]
  Ft_minus_5         <- FMiddleages[ages_20_55] * 
                          nLxf[ages_15_50, 2] / nLxf[ages_20_55, 2]
  names(Ft_minus_5)  <- ages_15_50
  
  Ft_minus_10        <- Ft_minus_5[ages_20_50] * 
                          nLxf[ages_15_45, 3] / nLxf[ages_20_50, 3]
  names(Ft_minus_10) <- ages_15_45
  
  # Now we take some averages to get to midpoints
  Ft_minus_.5        <- FMiddleages[ages_15_45] * .9 + Ft_minus_5[ages_15_45] * .1
  Ft_minus_2.5       <- FMiddleages[ages_15_45] * .5 + Ft_minus_5[ages_15_45] * .5
  Ft_minus_7.5       <- Ft_minus_5[ages_15_45] * .5 + Ft_minus_10[ages_15_45] * .5
  
  # 3 column matrix of sort-of-exposures for ages 15-45, matched to ASFR
  fExpos             <- cbind(Ft_minus_.5, Ft_minus_2.5, Ft_minus_7.5)

  # Calculate births
  Bt     <- colSums(fExpos * Asfr)


  #GenderCounts <- if (male) Males_five else Females_five
  Males_five_out     <- Males_five
  Females_five_out   <- Females_five
  ## Currently, this assumes that there can only be 3 dates.
  
  ## We only have 3 age groups to adjust and 3 dates
  PF <- 1 / (SRB + 1)
  
  # Age 0
  Females_five_out[1] <- Bt[1] * PF[1] * nLxf[1, 1] / radix
  Males_five_out[1]   <- Bt[1] * (1 - PF[1]) * nLxm[1, 1] / radix 
  
  # Age 1-4
  Females_five_out[2] <- Bt[2] * PF[2] * 5 *           
    sum(nLxf[1:2, 2]) / (radix * 5) - 
    Females_five_out[1] 
  
  Males_five_out[2]   <- Bt[2] * (1 - PF[2]) * 5 * 
    sum(nLxm[1:2, 2]) / (radix * 5) - 
    Males_five_out[1]
  
  # Age 5-9
  Females_five_out[3] <- Bt[3] * PF[3] * 5 * 
    sum(nLxf[1:2,3]) / (radix * 5) *
    nLxf[3,2] / sum(nLxf[1:2,2])
  
  Males_five_out[3]   <-  Bt[3] * (1 - PF[3]) * 5 * 
    sum(nLxm[1:2,3]) / (radix * 5) *
    nLxm[3,2] / sum(nLxm[1:2,2])
  
  # return the important things
  list(
    Females_adjusted = Females_five_out,
    Males_adjusted = Males_five_out,
    Females_five = Females_five,
    Males_five = Males_five,
    nLxf = nLxf,
    nLxm = nLxm,
    Asfr = Asfr,
    Exposure_female = fExpos,
    Bt = Bt,
    SRB = SRB,
    Age = Age,
    Ft_minus_.5 =Ft_minus_.5,
    Ft_minus_5= Ft_minus_5,
    Ft_minus_10 = Ft_minus_10
    )
}


# #' @rdname basepop_five
# #' @aliases basepop_five
# #' @param Females_single A named numeric vector. Reported population by 1-year age # groups for \code{refDate} for females. The names of the vector should reflect the age # groups. See examples. The method assumes that the last age group is open (for example# , the population ends at '80+' and '100+')
# #' @param Males_single A named numeric vector. Reported population by 1-year age # groups format \code{refDate} for males. The names of the vector should reflect the # age groups. See examples. The method assumes that the last age group is open (for # example, the population ends at '80+' and '100+')
# #'
# #' @export
# #'
# basepop_single <- function(country = NULL,
#                            refDate,
#                            Females_single,
#                            nLxFemale = NULL,
#                            nLxDatesIn = NULL,
#                            AsfrMat = NULL,
#                            AsfrDatesIn = NULL,
#                            ...,
#                            female = TRUE,
#                            SmoothedFemales = NULL,
#                            Males_single = NULL,
#                            nLxMale = NULL,
#                            SRB = 1.05,
#                            radix = NULL,
#                            verbose = TRUE) {
# 
#   stopifnot(
#     !is.null(names(Females_single)),
#     is_single(as.numeric(names(Females_single)))
#   )
# 
#   Females_abridged <- single2abridged(Females_single)
#   males_present <- !is.null(Males_single)
# 
#   if (males_present) {
#     stopifnot(
#       !is.null(names(Males_single)),
#       is_single(as.numeric(names(Males_single)))
#     )
# 
#     Males_abridged <- single2abridged(Males_single)
#     gender_single <- Males_single
#   }  else {
#     Males_abridged <- Males_single
#     gender_single <- Females_single
#   }
# 
#   res <-
#     basepop_five(
#       country = country,
#       refDate = refDate,
#       Females_five = Females_abridged,
#       nLxFemale = nLxFemale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       ... = ...,
#       female = female,
#       SmoothedFemales = SmoothedFemales,
#       Males_five = Males_abridged,
#       nLxMale = nLxMale,
#       SRB = SRB,
#       radix = radix
#     )
# 
#   # Since diff always returns a vector of length `length(x) - 1`,
#   # the 1 in the end is to reflct the the open ages for 80+ or 100+
#   AgeBins1 <- c(diff(as.integer(names(gender_single))), 1)
#   AgeBins2 <- c(diff(as.integer(names(res))), 1)
# 
#   rescaled_res <-
#     rescaleAgeGroups(
#       Value1 = gender_single,
#       AgeInt1 = AgeBins1,
#       Value2 = res,
#       AgeInt2 = AgeBins2,
#       splitfun = graduate_uniform
#     )
# 
#   round(rescaled_res, 3)
# }

#

# TR: modified to assume males and females always given
ArgsCheck <- function(ArgList) {
  with(ArgList, {
    stopifnot(
      is.numeric(Females_five),
      is.numeric(Males_five),
      is.numeric(SRB),
      is.matrix(nLxFemale),
      is.matrix(nLxMale),
      is.matrix(AsfrMat),
      #is.logical(female),
      ncol(nLxFemale) == length(nLxDatesIn),
      ncol(nLxMale) == length(nLxDatesIn)
      # TR no check on ASFRmat dates?
    )})
}


lt_infer_radix_from_1L0 <- function(L0){

  if (L0 > 1){
    radix_check <- L0 %>% as.integer() %>% log10()
    is_it_a_radix <- (radix_check - round(radix_check)) == 0
    
    if (!is_it_a_radix){
      pow <- L0 %>% round() %>% as.integer() %>% nchar()
      
      the_radix <- 10^pow
    } else {
      the_radix <- L0
    }
  } else {
    the_radix <- 1
  }
  the_radix
}

# TR: radix removed, as it seems lx was 1 but nLx was based on 1e5...
# will use indirect inference.
downloadnLx <- function(nLx, country, gender, nLxDatesIn) {
require(magrittr)
  verbose <- getOption("basepop_verbose", TRUE)
    if (!is.null(nLx)) {
      # TR: ensure colnames passed
      colnames(nLx) <- nLxDatesIn
      n <- nrow(nLx)
      Age <- c(0,1,seq(5,(n-2)*5,by=5))
      rownames(nLx) <- Age
      return(nLx)
    }
    
    if (is.null(nLx)){

    if (is.null(country)) stop("You need to provide a country to download the data for nLx")
      
      if (verbose) {
        cat(paste0("Downloading nLx data for ", country, ", years ", paste(nLxDatesIn,collapse=", "), ", gender ", gender), sep = "\n")
      }
  nLx <-
    lapply(nLxDatesIn, function(x) {
      fertestr::FetchLifeTableWpp2019(country, x, gender)$Lx
    }) %>% do.call("cbind",.)

  colnames(nLx) <- nLxDatesIn
  n <- nrow(nLx)
  Age <- c(0,1,seq(5,(n-2)*5,by=5))
  rownames(nLx) <- Age
  return(nLx)
   }
}

downloadAsfr <- function(Asfrmat, country, AsfrDatesIn) {

  verbose <- getOption("basepop_verbose", TRUE)

  if (!is.null(Asfrmat)) {
    # TR: can we assume colnames are AsfrDatesIn ?
    return(Asfrmat)
  }

  if (is.null(country)) stop("You need to provide a country to download the data for Asfrmat")

  tmp <-
    lapply(AsfrDatesIn, function(x) {

      if (verbose) {
        cat(paste0("Downloading Asfr data for ", country, ", year ", x), sep = "\n")
      }

      res <- fertestr::FetchFertilityWpp2019(country, x)["asfr"]
      names(res) <- NULL
      as.matrix(res)[2:nrow(res), , drop = FALSE]
    })

  Asfrmat <- do.call(cbind, tmp)
  colnames(Asfrmat) <- AsfrDatesIn
  Asfrmat
}


downloadSRB <- function(SRB, country, DatesOut){
  verbose <- getOption("basepop_verbose", TRUE)
  
  WPP2019_births <- DemoToolsData::WPP2019_births 
  # If not given and we have the country, then we use it
  if (is.null(SRB) & !is.null(country)){
    if (country %in% WPP2019_births$LocName){
    # TODO: really this should take a weighted average of SRB
    # over the period represented by each cetral date?
     
      SRB <- WPP2019_births %>% 
               dplyr::filter(LocName == country,
                             Year %in% floor(DatesOut)) %>% 
               dplyr::pull(SRB)
    }  else {
      if (verbose){
        cat(paste(country,"not available in WPP LocName list\n"))
      }
    }
    # otherwise will need to assume
  } 
  
  # if still not given then assume something
  if (is.null(SRB)){
    SRB <- rep(1.05,3)
    if (verbose){
      cat(paste(country,"not available in WPP LocName list\n"))
    }
  }
    
  # if given but not with 3 elements then repeat and cut as necessary
  if (is.numeric(SRB) & length(SRB) != 3){
    SRB <- rep(SRB,3)[1:3]
  }
  names(SRB) <-  DatesOut
  # return, potentially the same as input
  SRB
}





