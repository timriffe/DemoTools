#' The BPA and BPE methods for adjusting age groups under 10
#' @description Adjust population counts for the age groups 0, 1-4 and 5-9
#' @details
#' \code{basepop} can estimate both the BPA and BPE methods. If the user
#' specifies \code{SmoothedFemales}, \code{basepop} will return the BPA
#' method. If \code{SmoothedFemales} is left empty, \code{basepop} will
#' adjusting using the BPE method.
#'
#' Adjusting the the female population counts is the default. For this, all
#' arguments except \code{Males} and \code{nlxMale} are needed. However, for
#' adjusting the male population counts, the user needs to specify all arguments
#' as well as \code{Males}, \code{nlxMale} and \code{nlxFemale} and set
#' \code{sex == "m"}. 
#'
#' Currently, \code{basepop} only works with abridged five year age groups.
#'
#' @section BPA:
#'
#' Description:
#'
#' The method estimates a smoothed population ages 10 and over and adjusts the population under age 10 using the smoothed population and estimates of fertility and mortality.
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
#' (1) The method adjusts under-10 population to be consistent with fertility and mortality levels and adjusted adult female population.  
#'
#' Limitations:
#' 
#' (1) BPA assumes a linear change in fertility and mortality during the decade prior to the reference year.
#' 
#' (2) The procedure ignores migration, which can lead to misleading results.  There are two issues. First, age groups 0-4 and 5-9 are subject to migration, which would affect the comparability of estimated and reported populations in the base year.  Second, the estimated size of age groups 0-4 and 5-9 are calculated from numbers of women of reproductive age in the base year rejuvenated to points in the past.  With migration, rejuvenated number of women may exceed or be smaller than the number present, and giving birth to children, in the decade prior to the base year.
#' 
#' (3) BPA’s smoothing calculations may mask unusual, but real, variations in population age group size.  Smoothing irregularities in age structure not attributable to age misreporting will distort estimated births and survived children in the base year.
#'
#' Assumptions:
#'
#' (1) No significant international migration took place within the reference periods for the population, mortality, and fertility input.
#' 
#' (2) The data input as the "reported" population is not affected by underenumeration of persons in certain ages, nor by age misreporting.  
#'
#' @section BPE:
#'
#' Description:
#'
#' The method adjusts the population under age 10 using the reported population ages 10 and above and estimates of fertility and mortality.
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
#' (1) The method adjusts the under-10 population to be consistent with fertility and mortality levels and adult female population.  
#'
#' Limitations:
#'
#' (1) BPE assumes a linear change in fertility and mortality during the decade prior to the reference year.
#' 
#' (2) The procedure ignores migration, which can lead to misleading results.  There are two issues.  First, age groups 0-4 and 5-9 are subject to migration, which would affect the comparability of estimated and reported populations in the base year.  Second, the estimated size of age groups 0-4 and 5-9 are calculated from numbers of women of reproductive age in the base year rejuvenated to points in the past.  With migration, rejuvenated number of women may exceed or be smaller than the number present, and giving birth to children, in the decade prior to the base year.
#' 
#' (3) The method does not adjust for possible underenumeration and age misreporting errors in the over-10 “reported” population. If the reported population is subject to age-misreporting or age-sex-specific underenumeration, the over-10 population should be smoothed or otherwise corrected prior to use.
#' 
#' Assumptions:
#'
#' (1) No significant international migration took place within the reference periods for the population, mortality, and fertility input.
#' 
#' (2) The data input as the “reported” population is not affected by underenumeration of persons in certain ages, nor by age misreporting.  
#' 
#' @return A numeric vector of adjusted counts for the age groups 0, 1-4 and
#' 5-9 for either \code{Females} or \code{Males}, depending on whether the
#' user specified \code{"m"} or \code{"f"} in \code{Sex}. The remaining age
#' counts are left as is.
#' 
#' @param refYear An integer. The decimal reference year for which the reported population pertain
#' @param Females A numeric vector. Reported population by abridged 5-year lower bound age groups for females pertaining to \code{refYear}. This means that age 0 must be kept as a single age, ages 1-4 must grouped together as abridged age 1, and thereafter 5-year age groups
#' @param SRB A numeric. Sex ratio at birth (males / females)
#' @param nlxFemale A numeric matrix. The female nLx function of two life tables with ages in the rows and time in columns. The earlier date should be at least 7.5 years before the reference date of the "reported" population. The later date should be no earlier than one-half year before the reference date of the "reported" population
#' @param nlxDatesIn Vector of dates. The dates which pertain to the columns in \code{nlxFemale} and \code{nlxMale}. For details on the format of the dates, see the details of \code{\link{interp}}
#' @param asfrMat A numeric matrix. An age-period matrix of age specific fertility rates with age in rows, time in columns. 
#' @param AsfrDatesIn A vector of dates. The dates which pertain to the columns in \code{asfrMat}. For details on the format of the dates, see the details of \code{\link{interp}}
#' @param ... Arguments passed to \code{\link{interp}}. In particular, users might be interested in changing the interpolation method. By default, it's linearly interpolated
#' @param Sex A character string. Either 'f' or 'm', for adjusting either \code{Females} or \code{Males}. See the details section for the mandatory arguments for adjusting the male counts.
#' @param SmoothedFemales A numeric vector. This should be the result of passing the \code{Females} arguments through \code{smooth_age_5}. It is assumed that this argument is a numeric vector with the lower bound age group definitions as names. By default, this is what \code{smooth_age_5} returns. See the example section for an example.
#' @param Males A numeric vector. Reported population by abridged 5-year lower bound age groups for males pertaining to \code{refYear}. This means that age 0 must be kept as a single age, ages 1-4 must grouped together as abridged age 1, and thereafter 5-year age groups. This argument is only used when \code{Sex} is set to \code{"m"}
#' @param nlxMale A numeric matrix. The male nLx function of two life tables with ages in the rows and time in columns. The male nLx functions should be for age groups 0, 1 to 4, and 5 to 9. The dates which are represented in the columns are assumed to be the same as \code{nlxDatesIn}. This argument is only used when \code{Sex} is set to \code{"m"}
#' 
#' @export
#' @examples
#'
#'  \dontrun{
#'  # (1) RefYear
#'  refYear <- 1986.21
#'  
#'  # (2) Reported population by 5-year age groups and sex in the base year (Include unknowns).
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
#'  nlxDatesIn <- c(1977.31, 1986.50)
#'  
#'  nlxMale <- matrix(c(87732,
#'                      304435,
#'                      361064,
#'                      88451,
#'                      310605,
#'                      370362
#'                      ),
#'                    nrow = 3, ncol = 2)
#'  
#'  nlxFemale <- matrix(c(89842,
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
#'  # Sex = "m" with Males and nlxMale.
#'  bpa_male <-
#'    basepop(
#'      refYear = refYear,
#'      Males = pop_male_counts,
#'      Females = pop_female_counts,
#'      SmoothedFemales = smoothed_females,
#'      SRB = sex_ratio,
#'      nlxFemale = nlxFemale,
#'      nlxMale = nlxMale,
#'      nlxDatesIn = nlxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn,
#'      Sex = "m"
#'      )
#' 
#'  # See adjustments?
#'  pop_male_counts[1:3]
#'  bpa_male[1:3]
#' 
#'  # Adjusting the BPA for females requires less arguments
#'  
#'  bpa_female <-
#'    basepop(
#'      refYear = refYear,
#'      Females = pop_female_counts,
#'      SmoothedFemales = smoothed_females,
#'      SRB = sex_ratio,
#'      nlxFemale = nlxFemale,
#'      nlxDatesIn = nlxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn,
#'      Sex = "f"
#'    )
#' 
#'  pop_female_counts[1:3]
#'  bpa_female[1:3]
#' 
#'  # For adjustment using BPE, we use exactly the same definitions as above
#'  # but remove SmoothedFemales.
#'  bpe_male <-
#'    basepop(
#'      refYear = refYear,
#'      Males = pop_male_counts,
#'      Females = pop_female_counts,
#'      SRB = sex_ratio,
#'      nlxFemale = nlxFemale,
#'      nlxMale = nlxMale,
#'      nlxDatesIn = nlxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn,
#'      Sex = "m"
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
#'    basepop(
#'      refYear = refYear,
#'      Females = pop_female_counts,
#'      SRB = sex_ratio,
#'      nlxFemale = nlxFemale,
#'      nlxDatesIn = nlxDatesIn,
#'      asfrMat = asfrmat,
#'      AsfrDatesIn = AsfrDatesIn,
#'      Sex = "f"
#'    )
#' 
#'  pop_female_counts[1:3]
#'  bpa_female[1:3]
#'  bpe_female[1:3]
#'  }
#' 
#' @references
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertRef{PAS}{DemoTools}
basepop <- function(refYear,
                    Females,
                    SRB,
                    nlxFemale,
                    nlxDatesIn,
                    asfrMat,
                    AsfrDatesIn,
                    ...,
                    Sex = c("f", "m"),
                    SmoothedFemales = NULL,
                    Males = NULL,
                    nlxMale = NULL) {

  Sex <- match.arg(Sex)
  
  ## Check all arguments
  AllArgs <- as.list(environment())
  ArgsCheck(AllArgs)

  # We calculate basepop for males **ONLY** if the nlxMale
  # argument is not null
  nlxGender <- if (Sex == "m") nlxMale else nlxFemale
  male <- if (Sex == "m") TRUE else FALSE
  DatesOut <- sapply(c(0.5, 2.5, 7.5), function(x) refYear - x)

  # Interpolate the gender specific nlx to the requested
  # dates out
  nlx <- interp(
    nlxGender,
    datesIn = nlxDatesIn,
    datesOut = DatesOut,
    ...
  )
  # Turn the columns from the matrix into a list
  nlx <- lapply(as.data.frame(nlx), identity)

  # Interpolate only the female nlx to the requested dates.
  # This is independent of the nlx above, since to apply
  # basepop to either gender you need to have the female nlx
  nlxFemale <- interp(
    nlxFemale,
    datesIn = nlxDatesIn,
    datesOut = DatesOut,
    ...
  )
  # Turn the columns from the matrix into a list
  nlxFemale <- lapply(as.data.frame(nlxFemale), identity)

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
  
  EarliestDate <- which.max(names(nlxFemale))
  OlderNlxFemale <- nlxFemale[-EarliestDate]
  # Reorder the years such that the earlier year is first
  OlderNlxFemale <- OlderNlxFemale[sort(names(OlderNlxFemale), decreasing = TRUE)]

  # If SmoothedFemales is specified, we assume that it is returned
  # by smooth_age_5, so it should be named with the age groups.
  # This is checked in ArgList
  FemalePops <- if (!is.null(SmoothedFemales)) SmoothedFemales else Females

  # We use the smoothed vector names to only get certain age groups
  SmoothedFMiddleages <- FemalePops[as.character(seq(15, 55, by = 5))]

  # Currently, we assume that for calculating the estimated
  # female population for the DatesOut, the first OlderNlxFemale
  # date is multiplied by the SmoothedFMiddleages. However,
  # all dates after that one are mutltiplied by the previous
  # OlderNlxFemale in the order.
  for (i in seq_along(OlderNlxFemale)) {
    .x <- OlderNlxFemale[[i]]
    VecMult <- if (i == 1) SmoothedFMiddleages else OlderNlxFemale[[i-1]]
    IndLength <- if (i == 1) 0 else 1

    NlxDiv <- .x[4:(length(.x) - IndLength)]
    iter <- c(1, rep(seq(2, length(NlxDiv) - 1), each = 2), length(NlxDiv))
    NlxSeq <- lapply(seq(1, length(iter), by = 2), function(i) NlxDiv[iter[i:(i+1)]])

    # Here VecMult is either SmoothedFMiddleAges
    OlderNlxFemale[[i]] <- mapply(function(.y, .z) .y * .z[1] / .z[2],
                                  VecMult[-1],
                                  NlxSeq
                                  )
  }

  # We always assume that the FirstDate will be calculate with the first date
  # from the OlderNlxFemale, regardless of the number of dates.
  FirstDate <- list(0.2 * OlderNlxFemale[[1]] + 0.8 * SmoothedFMiddleages[-length(SmoothedFMiddleages)])
  names(FirstDate) <- max(names(nlxFemale))
  OlderNlxFemale <- c(
    FirstDate,
    OlderNlxFemale
  )
  
  # Currently, we assume that for estimating the population
  # for DatesOut, starting from the third year, we carry forward
  # the mean between the current year the previous one instead
  # of the SmoothedFMiddleages as fixed first year
  EstPop <- vector("list", length(OlderNlxFemale))
  names(EstPop) <- names(OlderNlxFemale)
  for (i in seq_along(OlderNlxFemale)) {
    .x <- OlderNlxFemale[[i]]
    MeanVec <- if (i > 2) OlderNlxFemale[[i-1]] else SmoothedFMiddleages

    # Loop along the 15-45 age groups
    EstPop[[i]] <- vapply(seq_along(seq(15, 45, by = 5)),
                           function(i) mean(c(MeanVec[i], .x[i])),
                           FUN.VALUE = numeric(1))
  }

  # This contains the estimate male/female population (depending on whether
  # the user specified nlxMale) for the DatesOut.
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

AdjustFemalePop <- function(x, asfr_pastyear, sex_ratio, male = TRUE) {
  births_age_female <- x * asfr_pastyear
  est_tot <- sum(births_age_female) * sex_ratio / (1 + sex_ratio)

  if (!male) est_tot <- sum(births_age_female) - est_tot

  est_tot
}

ArgsCheck <- function(ArgList) {
  with(ArgList, {
    stopifnot(is.numeric(Females),
              length(SRB) == 1,            
              is.numeric(SRB),
              is.matrix(nlxFemale),
              is.matrix(asfrMat),
              ncol(nlxFemale) == length(nlxDatesIn)
              )

    if ((Sex == "m")) {
      stopifnot(is.numeric(Males),
                is.matrix(nlxMale),
                ncol(nlxMale) == length(nlxDatesIn),
                !is.null(Males),
                !is.null(nlxMale))
    }
    
    is.named <- function(x) !is.null(names(x))
    if (!is.null(SmoothedFemales)) stopifnot(is.named(SmoothedFemales))

  })
}
