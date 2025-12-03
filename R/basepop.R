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
#' default. For this, only the \code{location}, \code{refDate} and
#' \code{Females_five} are needed. All other arguments are downloaded
#' or set to sensible defaults. For adjusting the male population
#' counts, the user needs to specify the \code{Males_five} population
#' counts and set \code{female = FALSE}.
#'
# # For \code{basepop_single}, the same procedure applies. The only difference
# # is that the vector \code{Males_five} is named \code{Males_single} and accepts
# # a vector of single ages rather than five year abridged age groups. Similarly,
# # the vector for females is \code{Females_single} rather than
# # \code{Females_five} and accepts single age groups.
#'
#' Currently, \code{basepop_five} works only with five year abridged age groups
# #' while \code{basepop_single} works only with single year age groups.
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
#' @return `basepop_five` returns a list with the following elements:
#'  *
#'  * `Females_adjusted` numeric vector of adjusted population counts for females. Age groups 0, 1-4, and 5-9 are adjusted, while ages 10 and higher are unchanged.
#'  * `Males_adjusted` numeric vector of adjusted population counts for males. Age groups 0, 1-4, and 5-9 are adjusted, while ages 10 and higher are unchanged.
#'  * `Females_five` numeric vector of female population counts given as input.
#'  * `Males_five` numeric vector of male population counts given as input.
#'  * `nLxf` numeric matrix of female `nLx`, abridged ages in rows and (potentially interpolated) time in columns. Potentially downloaded.
#'  * `nLxm` numeric matrix of male `nLx`, abridged ages in rows and (potentially interpolated) time in columns. Potentially downloaded.
#'  * `Asfr` numeric matrix of age specific fertility in 5-year age groups ages 15-19 until 45-49 in rows, and (potentially interpolated) time in columns. Potentially downloaded.
#'  * `Exposure_female` numeric matrix of approximated age-specific exposure in 5-year age groups ages 15-19 until 45-49 in rows, and (potentially interpolated) time in columns.
#'  * `Bt` births at three time points prior to census corresponding to the midpoints of the cohorts entering ages 0, 1-4, and 5-9.
#'  * `SRB` sex ratio at birth at three time points prior to census corresponding to the midpoints of the cohorts entering ages 0, 1-4, and 5-9. Potentially downloaded.
#'  * `Age` age groups of the input population counts.
#'
# #' `basepop_single` is used, the return value is a numeric vector with
# #' **single year age groups** where the counts between 0 and 10 are adjusted.
#'
#' @param location UN Pop Division `LocName` or `LocID`
#' @param refDate The reference year for which the reported population pertain
#' (these are the population counts in `Females_five` and
#' \code{Males_five}). Can either be a decimal date, a `Date` class.
#' If \code{nLxDatesIn} or \code{AsfrDatesIn} are not supplied and the
#' corresponding \code{nLxFemale/Male}/\code{AsfrMat} is not supplied,
#' \code{refDate} must be at a minimum 1962.5. This is because we can only
#' fetch WPP data from 1955 onwards, and these minimum date is assumed to be
#' 7.5 years before \code{refDate}, meaning 1955.
#'
#' @param Age integer vector of lower bounds of abridged age groups given in `Females_five` and `Males_five`.
#'
#' @param Females_five A named numeric vector with the population counts for
#' five-year abridged age groups for females in `refDate`. The names of the
#' vector should reflect the age groups. See the example section for some
#' examples.
#'
#' @param nLxFemale A numeric matrix. The female nLx function of two abridged life tables
#' with ages in the rows and time in columns. The earlier date should be at least
#' 7.5 years before the reference date of the "reported" population. The later
#' date should be no earlier than one-half year before the reference date of
#' the "reported" population. If not provided, it's automatically downloaded if
#' `location`, `refDate` and the equivalent population counts
#' `*_five` are provided.
#'
#' @param nLxDatesIn A vector of numeric years (for example, 1986). The dates
#' which pertain to the columns in `nLxFemale` and `nLxMale`. If not
#' provided, the function automatically determines two dates which are 8 years
#' before `refDate` and 0.5 years after `refDate`.
#'
#' @param AsfrMat A numeric matrix. An age-period matrix of age specific
#' fertility rates with age in rows, time in columns. If not provided, the
#' function automatically downloads the ASFR matrix based on the dates in
#' `AsfrDatesIn`.
#'
#' @param AsfrDatesIn A vector of numeric years (for example, 1986). These are
#' the dates which pertain to the columns in `AsfrMat`. If not provided,
#' the function automatically determines two dates which are 8 years before
#' `refDate` and 0.5 before `refDate`.
#'
#' @param ... Arguments passed to `\link{interp}`. In particular, users
#' might be interested in changing the interpolation method for the `nLx*`
#' matrices and the `Asfr` matrix. By default, it's linearly interpolated.
#'
#' @param Males_five A named numeric vector with the population counts for
#' five-year abridged age groups for males in `refDate`. The names of
#' the vector should reflect the age groups. See the example section for
#' some examples.
#'
#' @param nLxMale A numeric matrix. The male nLx function of two abridged life tables
#' with ages in the rows and time in columns. The dates which are represented
#' in the columns are assumed to be the same as `nLxDatesIn`. This
#' argument is only used when `female` is set to `FALSE` and
#' `Males_five` is provided. If `Males_five` is provided and
#' `female` set to `FALSE`, the `nLx` for males is
#' automatically downloaded for the dates in `nLxDatesIn`.
#'
#' @param SRB A numeric. Sex ratio at birth (males / females). Default is set
#' to 1.046. Only a maximum of three values permitted.
#'
#' @param SRBDatesIn A vector of numeric years (for example, 1986). Only a maximum
#' number of three dates allowed. These are
#' the dates which pertain to the values in `SRB`. If not provided,
#' the function automatically determines three dates which are 7.5 years,
#' 2.5 and 0.5 years before `refDate`.
#'
#' @param radix starting point to use in the adjustment of the three first age
#' groups. Default is NULL. If not provided, it is inferred based on the scale of age `1L0`.
#'
#' @param verbose when downloading new data, should the function print details
#' about the download at each step? Defaults to `TRUE`. We recommend the
#' user to set this to `TRUE` at all times because the function needs to
#' make decisions (such as picking the dates for the Asfr and nLx) that the user
#' should be aware of.
#'
#' @export
#' @examples
#'
#'  \dontrun{
#'refDate  <- 1985
#'location <- "Brazil"
#'data("popF",  package = "wpp2024")
#'data("popM",  package = "wpp2024")
#'
#'Age <- parse_number(unique(popF$age))
#'
#'# so technically this function is for abridged data
#'# it does not work if first group is in 5 years
#'pop_female_counts <- popF %>%
#'  dplyr::filter(name == location) %>% 
#'  pull("1985")
#'
#'pop_male_counts   <- popM %>%
#'  dplyr::filter(name == location) %>% 
#'  pull("1985")
#'
#'names(pop_female_counts) <- names(pop_female_counts) <-  Age
#'
#'# graduate
#'Females_five <- graduate_pclm(pop_female_counts, Age = Age) %>% 
#'  single2abridged()
#'Males_five   <- graduate_pclm(pop_male_counts, Age = Age) %>% 
#'  single2abridged()
#'
#'Age <- names2age(Males_five)
#'
#'# Automatically downloads the nLx, ASFR, and SRB data
#'bpe <- basepop_five(
#'  location     = location,
#'  refDate      = refDate,
#'  Females_five = Females_five,
#'  Males_five   = Males_five,
#'  Age          = Age,
#'  verbose      = TRUE
#')
#'
#'  }
#'
#' @references
#' \insertRef{arriaga1994population}{DemoTools}
#' \insertRef{PAS}{DemoTools}
#'
# basepop_five <- function(location = NULL,
#                          refDate,
#                          Age = NULL,
#                          Females_five,
#                          Males_five = NULL,
#                          nLxFemale = NULL,
#                          nLxMale = NULL,
#                          nLxDatesIn = NULL,
#                          AsfrMat = NULL,
#                          AsfrDatesIn = NULL,
#                          ...,
#                          SRB = NULL,
#                          SRBDatesIn = NULL,
#                          radix = NULL,
#                          verbose = TRUE) {
# 
#   options(basepop_verbose = verbose)
#   on.exit(options(basepop_verbose = NULL))
# 
#   # Ensure census date is numeric.
#   # "YYYY-MM-DD" input is acceptable
#   refDate <- dec.date(refDate)
# 
#   if (!is.null(Age)){
#     stopifnot(is_abridged(Age))
#     stopifnot(length(Age) == length(Females_five))
#   } else {
#     if (!is.null(names(Females_five))){
#       Age <- names2age(Females_five)
#     } else {
#       if (verbose) {
#         cat("Assuming age groups are in standard abridged intervals")
#       }
#       # last resort = assume in abrided ages!
#       Age <- inferAgeIntAbr(Females_five)
#     }
#   }
# 
#   if (is.null(nLxDatesIn)) {
#     # re PJ issue #183 suggested default
#     nLxDatesIn <- refDate - c(0.5, 7.5)
#     #nLxDatesIn <- c(abs(8 - refDate), refDate + 0.5)
#     if (verbose) {
#       cat(paste0("Assuming the two prior dates for the nLx matrix to be: ", paste0(nLxDatesIn, collapse = ", ")), sep = "\n")
#     }
#   }
# 
#   if (is.null(AsfrDatesIn)) {
#     # re PJ issue #183 suggested default
#     AsfrDatesIn <- refDate - c(0.5, 7.5)
#     #AsfrDatesIn <- abs(c(8, 0.5) - refDate)
#     if (verbose) {
#       cat(paste0("Assuming the two prior dates for the Asfr matrix to be: ", paste0(AsfrDatesIn, collapse = ", ")), sep = "\n")
#     }
#   }
# 
#   # ensure vectors named, for purposes of selection
#   names(Females_five) <- Age
#   names(Males_five)   <- Age
#   ## obtain nLx for males and females
#   ## If these arguments have been specified, they return
#   ## the same thing and don't download the data
#   nLxFemale <-
#     downloadnLx(
#       nLx = nLxFemale,
#       location = location,
#       gender = "female",
#       nLxDatesIn = nLxDatesIn
#     )
# 
#   nLxMale <-
#     downloadnLx(
#       nLx = nLxMale,
#       location = location,
#       gender = "male",
#       nLxDatesIn = nLxDatesIn
#     )
# 
#   if (is.null(radix)) {
#     # TR: not perfect, but it's a better guess. It would seem the radix
#     # being pulled before was always 1, whereas the nLx columns was based on 100000
#     radix <- lt_infer_radix_from_1L0(nLxMale[1,1])
#     if (verbose) {
#       cat(paste0("Setting radix to value of lx: ", radix, ". Can be overwritten with the `radix` argument"), sep = "\n")
#     }
#   }
# 
#   AsfrMat <-
#     downloadAsfr(
#       Asfrmat = AsfrMat,
#       location = location,
#       AsfrDatesIn = AsfrDatesIn
#     )
# 
#   DatesOut <- refDate - c(0.5, 2.5, 7.5)
#   SRBDatesIn <- if (!is.null(SRBDatesIn)) SRBDatesIn else DatesOut
# 
#   SRB <- downloadSRB(SRB,
#                      location,
#                      DatesOut = SRBDatesIn,
#                      verbose = verbose)
# 
#   ## Check all arguments
#   AllArgs <- as.list(environment())
#   ArgsCheck(AllArgs)
# 
#   lower_bound <- abs(min(nLxDatesIn) - min(DatesOut))
#   upper_bound <- abs(max(nLxDatesIn) - max(DatesOut))
# 
#   if (lower_bound > 5 || upper_bound > 5) {
#     stop("nLxDatesIn implies an extrapolation of > 5 years to achieve the needed reference dates")
#   }
# 
#   # Interpolate the gender specific nLx to the requested
#   # dates out
#   nLxf <- interp(
#     nLxFemale,
#     datesIn = nLxDatesIn,
#     datesOut = DatesOut,
#     ...
#   )
# 
#   nLxm <- interp(
#     nLxMale,
#     datesIn = nLxDatesIn,
#     datesOut = DatesOut,
#    ...
#   )
# 
#   lower_bound <- abs(min(AsfrDatesIn) - min(DatesOut))
#   upper_bound <- abs(max(AsfrDatesIn) - max(DatesOut))
# 
#   if (lower_bound > 5 || upper_bound > 5) {
#     stop("AsfrDatesIn implies an extrapolation of > 5 years to achieve the needed reference dates")
#   }
# 
#   # Interpolate the asfr to the requested dates.
#   # This is gender agnostic.
#   Asfr <- interp(
#     AsfrMat,
#     datesIn = AsfrDatesIn,
#     datesOut = DatesOut,
#     ...
#   )
# 
#   # TR: Follows spreadsheet logic, can still be more elegant.
#   # sometimes character indexing, sometimes position, but still
#   ages_15_55         <- as.character(seq(15,55,by=5))
#   ages_20_55         <- ages_15_55[-1]
#   ages_15_50         <- ages_15_55[-9]
#   ages_20_50         <- ages_15_55[-c(1,9)]
#   ages_15_45         <- ages_15_55[-c(8,9)]
#   ages_20_45         <- ages_15_55[-c(1,8,9)]
#   ages_15_40         <- ages_15_55[-c(7,8,9)]
# 
#   FMiddleages        <- Females_five[ages_15_55]
#   Ft_minus_5         <- FMiddleages[ages_20_55] *
#                           nLxf[ages_15_50, 2] / nLxf[ages_20_55, 2]
#   names(Ft_minus_5)  <- ages_15_50
# 
#   Ft_minus_10        <- Ft_minus_5[ages_20_50] *
#                           nLxf[ages_15_45, 3] / nLxf[ages_20_50, 3]
#   names(Ft_minus_10) <- ages_15_45
# 
#   # Now we take some averages to get to midpoints
#   Ft_minus_.5        <- FMiddleages[ages_15_45] * .9 + Ft_minus_5[ages_15_45] * .1
#   Ft_minus_2.5       <- FMiddleages[ages_15_45] * .5 + Ft_minus_5[ages_15_45] * .5
#   Ft_minus_7.5       <- Ft_minus_5[ages_15_45] * .5 + Ft_minus_10[ages_15_45] * .5
# 
#   # 3 column matrix of sort-of-exposures for ages 15-45, matched to ASFR
#   fExpos             <- cbind(Ft_minus_.5, Ft_minus_2.5, Ft_minus_7.5)
# 
#   # Calculate births
#   Bt     <- colSums(fExpos * Asfr)
# 
# 
#   #GenderCounts <- if (male) Males_five else Females_five
#   Males_five_out     <- Males_five
#   Females_five_out   <- Females_five
#   ## Currently, this assumes that there can only be 3 dates.
# 
#   ## We only have 3 age groups to adjust and 3 dates
#   PF <- 1 / (SRB + 1)
# 
#   # Age 0
#   Females_five_out[1] <- Bt[1] * PF[1] * nLxf[1, 1] / radix
#   Males_five_out[1]   <- Bt[1] * (1 - PF[1]) * nLxm[1, 1] / radix
# 
#   # Age 1-4
#   Females_five_out[2] <- Bt[2] * PF[2] * 5 *
#     sum(nLxf[1:2, 2]) / (radix * 5) -
#     Females_five_out[1]
# 
#   Males_five_out[2]   <- Bt[2] * (1 - PF[2]) * 5 *
#     sum(nLxm[1:2, 2]) / (radix * 5) -
#     Males_five_out[1]
# 
#   # Age 5-9
#   Females_five_out[3] <- Bt[3] * PF[3] * 5 *
#     sum(nLxf[1:2,3]) / (radix * 5) *
#     nLxf[3,2] / sum(nLxf[1:2,2])
# 
#   Males_five_out[3]   <-  Bt[3] * (1 - PF[3]) * 5 *
#     sum(nLxm[1:2,3]) / (radix * 5) *
#     nLxm[3,2] / sum(nLxm[1:2,2])
# 
#   # return the important things
#   list(
#     Females_adjusted = Females_five_out,
#     Males_adjusted = Males_five_out,
#     Females_five = Females_five,
#     Males_five = Males_five,
#     nLxf = nLxf,
#     nLxm = nLxm,
#     Asfr = Asfr,
#     Exposure_female = fExpos,
#     Bt = Bt,
#     SRB = SRB,
#     Age = Age
#     )
# }

basepop_five <- function(location     = NULL,
                         refDate      = NULL,
                         Age          = NULL,
                         Females_five = NULL,
                         Males_five   = NULL,
                         nLxFemale    = NULL,
                         nLxMale      = NULL,
                         nLxDatesIn   = NULL,
                         AsfrMat      = NULL,
                         AsfrDatesIn  = NULL,
                         SRB          = NULL,
                         SRBDatesIn   = NULL,
                         radix        = NULL,
                         verbose      = TRUE,
                         ...) {
  
  options(basepop_verbose = verbose)
  on.exit(options(basepop_verbose = NULL))
  refDate <- dec.date(refDate)
  
  # if population are missing stop
  
  if(is.null(Females_five) | is.null(Males_five)) { 
    
    stop("Male or female population was not provided.")
    
  }
  
  
  if (length(Females_five) != length(Males_five)) {
    stop("Females population length differs from male population length.")
  }
  
  
  
  if(!is.null(Age)) {
    
    stopifnot(is_abridged(Age))
    stopifnot(length(Age) == length(Females_five))
    
  } else {
    
    if(!is.null(names(Females_five))) {
      
      Age <- names2age(Females_five)
      
    } else {
      
      if(verbose) {
        
        message("Assuming age groups are in standard abridged intervals")
        
      }
      
      Age <- inferAgeIntAbr(Females_five)
      
    }
  }
  
  if(is.null(nLxDatesIn)) {
    
    nLxDatesIn <- refDate - c(0.5, 7.5)
    
    if(verbose) {
      
      message(paste0(
        "Assuming the two prior dates for the nLx matrix to be: ",
        paste0(nLxDatesIn, collapse = ", ")
      ),
      sep = "\n")
      
    }
  }
  
  if(is.null(AsfrDatesIn)) {
    
    AsfrDatesIn <- refDate - c(0.5, 7.5)
    
    if(verbose){
      
      message(paste0(
        "Assuming the two prior dates for the Asfr matrix to be: ",
        paste0(AsfrDatesIn, collapse = ", ")
      ),
      sep = "\n")
      
    }
  }
  
  names(Females_five) <- Age
  names(Males_five)   <- Age
  
  if(is.null(radix)) { 
    
    radix       <- 1
    check_radix <- FALSE
    
    message("Since radix was not provided we assume your radix to be equal to 1")
    
  } else { 
    
    check_radix <- TRUE
    
  }
  
  if(is.null(nLxFemale)) { 
    # download nLx data
    nLxFemale <- downloadnLx(
      nLx        = nLxFemale,
      Age        = Age,
      location   = location,
      gender     = "female",
      nLxDatesIn = sort(nLxDatesIn),
      output     = "abridged",
      radix      = radix,
      verbose    = verbose
    )
  }
  
  if(is.null(nLxMale)) { 
    
    nLxMale <- downloadnLx(
      nLx        = nLxMale,
      Age        = Age,
      location   = location,
      gender     = "male",
      nLxDatesIn = sort(nLxDatesIn),
      output     = "abridged",
      radix      = radix,
      verbose    = verbose
    )
    
  }
  
  
  if(check_radix) {
    
    radix1 <- radix
    radix  <- lt_rule_l0_1L0(nLxMale[1, 1])
    
    if(radix1 != radix) { 
      
      stop("Your provided radix value do not match with the radix obtained from your nLx table. Please enter correct radix or nLx data.") 
      
    }
  }
  
  if(is.null(AsfrMat)) { 
    # asfr now
    AsfrMat <- downloadASFR(Asfrmat     = AsfrMat,
                            location    = location,
                            AsfrDatesIn = sort(AsfrDatesIn),
                            Age         = NULL,
                            output      = "5-year",
                            verbose     = verbose)
    
  }
  
  DatesOut   <- refDate - c(0.5, 2.5, 7.5)
  
  # SRB now
  SRBDatesIn <- if(!is.null(SRBDatesIn)) {
    
    SRBDatesIn
    
  } else { 
    
    DatesOut
    
  }
  
  if(is.null(SRB)) {
    SRB <- downloadSRB(SRB       = SRB,
                       location  = location,
                       DatesOut  = sort(SRBDatesIn),
                       verbose   = verbose)
  }
  # check argumnts if everything is ok
  AllArgs <- as.list(environment())
  ArgsCheck(AllArgs)
  
  # check bounds
  lower_bound <- abs(min(nLxDatesIn) - min(DatesOut))
  upper_bound <- abs(max(nLxDatesIn) - max(DatesOut))
  
  if(lower_bound > 5 || upper_bound > 5) {
    
    stop(
      "nLxDatesIn implies an extrapolation of > 5 years to achieve the needed reference dates"
    )
    
  }
  
  # same for asfr
  lower_bound <- abs(min(AsfrDatesIn) - min(DatesOut))
  upper_bound <- abs(max(AsfrDatesIn) - max(DatesOut))
  
  if(lower_bound > 5 || upper_bound > 5) {
    stop(
      "AsfrDatesIn implies an extrapolation of > 5 years to achieve the needed reference dates"
    )
  }
  
  # interpolate nLx and ASFR
  nLxf <- interp(nLxFemale,
                 datesIn  = nLxDatesIn, 
                 datesOut = DatesOut)
  nLxm <- interp(nLxMale,
                 datesIn  = nLxDatesIn, 
                 datesOut = DatesOut)
  Asfr <- interp(AsfrMat, 
                 datesIn  = AsfrDatesIn, 
                 datesOut = DatesOut)
  
  # now here a little update
  # we need ASFR match fExpos for colSums(fExpos * Asfr)
  # thus age adjustment
  Asfr <- Asfr[c(which(names2age(Asfr) %in% seq(from = 15, to = 45, by = 5))), ]
  # Asfr                <- Asfr[-c(1, ncol(Asfr)), ]
  ages_15_55          <- as.character(seq(15, 55, by = 5))
  ages_20_55          <- ages_15_55[-1]
  ages_15_50          <- ages_15_55[-9]
  ages_20_50          <- ages_15_55[-c(1, 9)]
  ages_15_45          <- ages_15_55[-c(8, 9)]
  ages_20_45          <- ages_15_55[-c(1, 8, 9)]
  ages_15_40          <- ages_15_55[-c(7, 8, 9)]
  FMiddleages         <- Females_five[ages_15_55]
  Ft_minus_5          <- FMiddleages[ages_20_55] * nLxf[ages_15_50, 2] / nLxf[ages_20_55, 2]
  names(Ft_minus_5)   <- ages_15_50
  Ft_minus_10         <- Ft_minus_5[ages_20_50]  * nLxf[ages_15_45, 3] / nLxf[ages_20_50, 3]
  names(Ft_minus_10)  <- ages_15_45
  Ft_minus_.5         <- FMiddleages[ages_15_45] * 0.9 + Ft_minus_5[ages_15_45]  * 0.1
  Ft_minus_2.5        <- FMiddleages[ages_15_45] * 0.5 + Ft_minus_5[ages_15_45]  * 0.5
  Ft_minus_7.5        <- Ft_minus_5[ages_15_45]  * 0.5 + Ft_minus_10[ages_15_45] * 0.5
  fExpos              <- cbind(Ft_minus_.5, Ft_minus_2.5, Ft_minus_7.5)
  Bt                  <- colSums(fExpos * Asfr)
  Males_five_out      <- Males_five
  Females_five_out    <- Females_five
  PF                  <- 1 / (SRB + 1)
  Females_five_out[1] <- Bt[1] * (PF[1])     * nLxf[1, 1] / radix
  Males_five_out[1]   <- Bt[1] * (1 - PF[1]) * nLxm[1, 1] / radix
  Females_five_out[2] <- Bt[2] * (PF[2] * 5) * sum(nLxf[1:2, 2])     / (radix * 5) - Females_five_out[1]
  Males_five_out[2]   <- Bt[2] * (1 - PF[2]) * 5 * sum(nLxm[1:2, 2]) / (radix * 5) - Males_five_out[1]
  Females_five_out[3] <- Bt[3] * (PF[3] * 5) * sum(nLxf[1:2, 3])     / (radix * 5) * nLxf[3, 2] / sum(nLxf[1:2, 2])
  Males_five_out[3]   <- Bt[3] * (1 - PF[3]) * 5 * sum(nLxm[1:2, 3]) / (radix * 5) * nLxm[3, 2] / sum(nLxm[1:2, 2])
  
  return(
    list(
      Females_adjusted = Females_five_out,
      Males_adjusted   = Males_five_out,
      Females_five     = Females_five,
      Males_five       = Males_five,
      nLxf             = nLxf,
      nLxm             = nLxm,
      Asfr             = Asfr,
      Exposure_female  = fExpos,
      Bt               = Bt,
      SRB              = SRB,
      Age              = Age
    ))
}






# Authors: Tim, Rustam
#' Construct Single-Year Base Population Using Reverse Survival
#'
#' Constructs a single-year, age-specific population for both sexes using the reverse
#' survival method. This combines life table data, fertility rates, and sex ratio
#' at birth (SRB) to estimate age-specific population counts, optionally using
#' user-supplied demographic data or WPP inputs.
#'
#' @param location Optional character string naming the location (used only for messages).
#' @param refDate Numeric reference date (e.g., census year). Required.
#' @param Age Optional numeric vector of single ages. Must match population vector length if supplied.
#' @param country_code Numeric ISO country code. Required for WPP data retrieval.
#' @param Females_single Optional numeric vector of female single-age population counts.
#' @param Males_single Optional numeric vector of male single-age population counts.
#' @param nLxFemale Optional life table \code{nLx} matrix for females. If supplied, used to derive period and cohort survival ratios.
#' @param nLxMale Optional life table \code{nLx} matrix for males. If supplied, used to derive period and cohort survival ratios.
#' @param nLxDatesIn Optional numeric vector of years corresponding to life table data.
#' @param AsfrMat Optional matrix of age-specific fertility rates (ASFR). If supplied, reshaped internally; if not, WPP ASFR will be downloaded.
#' @param AsfrDatesIn Optional numeric vector of years corresponding to ASFR data.
#' @param SRB Optional numeric vector or single value giving sex ratio at birth (males per female); if not supplied, WPP data are used.
#' @param SRBDatesIn Optional numeric vector of years corresponding to SRB input.
#' @param radix Numeric starting radix for life table calculations. Defaults to 1 if not inferred.
#' @param verbose Logical; if \code{TRUE}, prints progress messages. Default = \code{TRUE}.
#' @param ... Additional arguments passed to internal DemoTools functions.
#'
#' @details
#' The function performs the following sequence:
#' \enumerate{
#'   \item Validates input arguments and loads the most recent available \code{wpp} dataset.
#'   \item Retrieves population, mortality (\code{mx}), fertility (ASFR), and SRB data from WPP if not supplied by the user.
#'   \item Converts \code{nLx} inputs (if provided) into period life tables and reconstructs cohort survival ratios needed for reverse survival.
#'   \item Computes reverse-survival multipliers (\code{SxRev}) and cumulative inflation factors.
#'   \item Estimates exposures for reproductive ages and computes annual births.
#'   \item Splits births into male and female counts using the sex ratio at birth.
#'   \item Reverse-survives births to reconstruct ages 0–9, and reverse-survives older ages using cohort-specific inflation factors.
#'   \item Returns reconstructed population for males and females.
#' }
#'
#' The reverse survival method is useful when historical population data are unavailable
#' or incomplete but life table and fertility information are known.
#'
#' @return A named list with two tibbles:
#' \describe{
#'   \item{pop_hat_m}{Reconstructed male population by single year of age and cohort.}
#'   \item{pop_hat_f}{Reconstructed female population by single year of age and cohort.}
#' }
#' @importFrom dplyr select mutate arrange group_by ungroup right_join left_join summarize
#' @importFrom tidyr pivot_longer pivot_wider unnest
#' @importFrom tibble tibble rownames_to_column as_tibble_col
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' # Example: Estimate base population for Germany, year 2000
#' basepop <- basepop_single(
#'   country_code = 276,
#'   refDate = 2000
#' )
#'
#' basepop$pop_hat_m  # Male reconstructed population
#' basepop$pop_hat_f  # Female reconstructed population
#'
#' @export

basepop_single <- function(refDate        = NULL,
                           Age            = NULL,
                           country_code   = NULL,
                           Females_single = NULL,
                           Males_single   = NULL,
                           nLxFemale      = NULL,
                           nLxMale        = NULL,
                           nLxDatesIn     = NULL,
                           AsfrMat        = NULL,
                           AsfrDatesIn    = NULL,
                           SRB            = NULL,
                           SRBDatesIn     = NULL,
                           radix          = NULL,
                           verbose        = TRUE,
                           ...) {
  
  # --- Global setup ---------------------------------------------------------
  # options from DemoTools version
  options(basepop_verbose = verbose)
  on.exit(options(basepop_verbose = NULL))
  
  if(!is.numeric(country_code)) { 
    
    stop("Please provide a numeric code of the country you working with.")
    
  }
  
  
  if(length(refDate) > 1) { 
    
    stop("You should pick one reference date refDate")
    
  }
  
  
  # Jan 1 2000 female pop;
  # note:  1999 means dec 31, 1999, so we treat as jan 1, 2000.
  # we need these dates to filter period
  # refDate is treated as Jan 1 of the following year for consistency with WPP
  refDate       <- dec.date(refDate) - 1
  refDate_start <- refDate - 9
  
  # --- Basic validation -----------------------------------------------------
  # we need minimum of date and country to run a function
  if(is.null(refDate) | is.null(country_code)) {
    
    stop("At least refDate and country_code should be provided for function to work")
    
  }
  
  # here we attach the latest wpp package available
  installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
  
  if(length(installed_wpp) == 0) {
    
    stop("No wpp package installed.")
    
  }
  
  latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
  library(latest_wpp, character.only = TRUE, quietly = TRUE)
  
  if(parse_number(latest_wpp) < 2022) { 
    
    stop("No single ages are availabe in wpp versions earlier that wpp2022. Consider updating the wpp package or change to five year solution, basepop_five.")
    
  }
  
  # --- Load population and mortality data ----------------------------------
  # if user did not provide population we can calculate it from wpp
  if(is.null(Females_single) | is.null(Males_single)) {
    
    data("popAge1dt", package = latest_wpp)
    
  }
  
  # this data is for future calculations
  data("mx1dt", package = latest_wpp)
  
  # --- Prepare user-supplied population vectors ----------------------------
  # Convert numeric vectors to tidy tibbles with cohort tagging
  # if user provided Population vector we turn it into tibble
  if(!is.null(Females_single)) {
    
    Females_single <- tibble("pop"  = Females_single,
                             "year" = refDate + 1) %>% 
      mutate("age"    = row_number() - 1,
             "cohort" = year - age - 1)
    
    
  }
  
  # same for Males
  if(!is.null(Males_single)) {
    
    Males_single <- tibble("pop"    = Males_single,
                           "year"   = refDate + 1) %>% 
      mutate("age"    = row_number() - 1,
             "cohort" = year - age - 1)
    
  }
  
  # --- Retrieve WPP population if user did not provide one -----------------
  # if user did not provide Population
  if(is.null(Females_single)) {
    
    Females_single <- popAge1dt %>%
      dplyr::filter(.data$country_code == !!country_code, 
             .data$year == refDate) %>%
      select("year", "age", "pop" = "popF") %>%
      mutate("year"   = .data$year + 1,
             "cohort" = .data$year - .data$age - 1)
    
  }
  
  # same operation for males
  if(is.null(Males_single)) {
    
    Males_single <- popAge1dt %>% 
      dplyr::filter(.data$country_code == !!country_code, 
             .data$year == refDate) %>% 
      select("year", "age", "pop" = "popM") %>% 
      mutate("year"   = .data$year + 1,
             "cohort" = .data$year - .data$age - 1)
    
  }
  
  # --- Age validation and setup --------------------------------------------
  # Age setup
  if(!is.null(Age)) {
    
    stopifnot(is_single(Age))
    stopifnot(length(Age) == nrow(Females_single))
    
  } else {
    
    Age <- as.integer(sort(unique(Females_single$age)))
    
  }
  
  # --- Determine life table radix ------------------------------------------
  if(is.null(radix) & !is.null(nLxFemale)) {
    
    radix <- lt_infer_radix_from_1L0(nLxFemale[1, 1])
    
    if(verbose) {
      
      cat(paste0("Setting radix to value of lx: ", radix, 
                 ". Can be overwritten with the `radix` argument"),
          sep = "\n")
      
    }
    
  }
  
  if(is.null(radix) & is.null(nLxFemale)) { 
    
    radix <- 1
    cat(paste0("Setting radix to value of lx: ", radix, 
               ". Can be overwritten with the `radix` argument"),
        sep = "\n")
    
  }
  
  # --- Compute reverse survival ratios (SxRev) ------------------------------
  # These represent the probability of surviving backward in time.
  # Female and male SxRev matrices are created either from nLx input or WPP mx data.
  # Inflation factor is cumulative survival ratio (used for exposure reconstruction).
  
  if(!is.null(nLxFemale)) {
    
    mxF <- downloadnLx(nLx = nLxFemale) %>% 
      as.data.frame() %>% 
      { if (!"age" %in% names(.)) rownames_to_column(., "age") else . } %>%
      mutate("age" = parse_number(.data$age)) %>% 
      pivot_longer(-c("age"),
                   names_to  = "year",
                   values_to = "Lxp",
                   names_transform  = list("year" = as.numeric)) %>% 
      mutate("n"  = age2int(.data$age, OAvalue = Inf),
             "ax" = ifelse(.data$age == 0,
                           lt_rule_L0_1a0(L0_target = .data$Lxp[.data$age == 0],
                                          l0        = radix,
                                          sex       = "f"), 0.5), .by = "year") %>% 
      unnest("ax") %>% 
      mutate("qx"      = (.data$Lxp - lead(.data$Lxp)) / .data$Lxp,
             "lx"      = lt_id_q_l(nqx   = .data$qx, 
                                   radix = radix),
             "dx"      = .data$lx * .data$qx,
             "Lx"      = .data$lx - .data$dx * (1 - .data$ax),
             "mx"      = lt_id_qa_m(nqx    = .data$qx, 
                                    nax    = .data$ax, 
                                    AgeInt = .data$n),
             "cohort"  = .data$year - .data$age - 1,
             "qx"      = ifelse(.data$age == max(.data$age), 1, .data$qx),
             "mx"      = ifelse(.data$age == max(.data$age), 1 / 
                                  .data$ax[.data$age == max(.data$age)], .data$mx)) %>% 
      arrange(.data$cohort, .data$age)  %>% 
      mutate(
        "lx"    = lt_id_q_l(nqx   = .data$qx, 
                            radix = radix),
        "dx"    = .data$lx * .data$qx,
        "Lx"    = .data$lx - .data$dx * (1 - .data$ax),
        "SxRev" = .data$Lx / lead(.data$Lx, default = 1),
        .by = "cohort"
      )  %>% 
      arrange(.data$year, .data$age)  %>% 
      mutate(
        "lxp"   = lt_id_q_l(nqx   = .data$qx, 
                            radix = radix),
        "dxp"   = .data$lxp * .data$qx,
        "Lxp"   = .data$lxp - .data$dxp * (1 - .data$ax),
        "SxRev" = ifelse(year == max(.data$year), .data$Lxp / 
                           lead(.data$Lxp), .data$SxRev), 
        .by = "year")  %>% 
      dplyr::filter(.data$age < 100)  %>% 
      select("cohort", "year", "age", "SxRev")  %>% 
      arrange(.data$cohort, -c(.data$age)) %>% 
      # inflation for reverse-surviving
      mutate("inflation_factor" = cumprod(.data$SxRev), .by = "cohort")
    
  }
  
  if(!is.null(nLxMale)) {
    
    # tst <- downloadnLx(
    #   nLx         = NULL,
    #   location    = country_code,
    #   gender      = "m",
    #   nLxDatesIn  = refDate_start:refDate,
    #   Age         = NULL,
    #   method      = "linear",
    #   output      = "single",
    #   radix       = radix)
    
    mxM <- downloadnLx(nLx = nLxMale) %>% 
      as.data.frame() %>% 
      { if (!"age" %in% names(.)) rownames_to_column(., "age") else . } %>%
      mutate("age" = parse_number(.data$age)) %>% 
      pivot_longer(-c("age"),
                   names_to  = "year",
                   values_to = "Lxp",
                   names_transform  = list(year = as.numeric)) %>% 
      mutate("n"  = age2int(.data$age, OAvalue = Inf),
             "ax" = ifelse(.data$age == 0,
                           lt_rule_L0_1a0(L0_target = .data$Lxp[.data$age == 0],
                                          l0        = radix,
                                          sex       = "m"), 0.5), .by = "year") %>% 
      unnest("ax") %>% 
      mutate("qx"      = (.data$Lxp - lead(.data$Lxp)) / .data$Lxp,
             "lx"      = lt_id_q_l(nqx   = .data$qx, 
                                   radix = radix),
             "dx"      = .data$lx * .data$qx,
             "Lx"      = .data$lx - .data$dx * (1 - .data$ax),
             "mx"      = lt_id_qa_m(nqx    = .data$qx, 
                                    nax    = .data$ax, 
                                    AgeInt = .data$n),
             "cohort"  = .data$year - .data$age - 1,
             "qx"      = ifelse(.data$age == max(.data$age), 1, .data$qx),
             "mx"      = ifelse(.data$age == max(.data$age), 1 / 
                                  .data$ax[.data$age == max(.data$age)], .data$mx)) %>% 
      arrange(.data$cohort, .data$age)  %>%
      mutate(
        "lx"    = lt_id_q_l(nqx   = .data$qx, 
                            radix = radix),
        "dx"    = .data$lx * .data$qx,
        "Lx"    = .data$lx - .data$dx * (1 - .data$ax),
        "SxRev" = .data$Lx / lead(.data$Lx, default = 1),
        .by = "cohort"
      ) %>% 
      arrange(.data$year, .data$age) %>% 
      mutate(
        "lxp"   = lt_id_q_l(nqx   = .data$qx, 
                            radix = radix),
        "dxp"   = .data$lxp * .data$qx,
        "Lxp"   = .data$lxp - .data$dxp * (1 - .data$ax),
        "SxRev" = ifelse(year == max(.data$year), .data$Lxp / 
                           lead(.data$Lxp), .data$SxRev), 
        .by = "year") %>% 
      dplyr::filter(.data$age < 100) %>% 
      select("cohort", "year", "age", "SxRev") %>% 
      arrange(.data$cohort, -c(.data$age)) %>% 
      # inflation for reverse-surviving
      mutate("inflation_factor" = cumprod(.data$SxRev), .by = "cohort")
    
  }
  
  # if not provided by user
  if(is.null(nLxMale)) {
    
    mxM <- mx1dt  %>% 
      dplyr::filter(.data$country_code == !!country_code,
             between(.data$year, refDate_start, refDate)) %>% 
      as_tibble() %>% 
      select("year", "age", "mx" = "mxM") %>% 
      # could try to warp to PC shape here,
      # but uncertain infants. Maybe using
      # an a0 assumption it'd be doable.
      # need cohorts to structure reverse survival
      mutate(
        "cohort" = .data$year - .data$age - 1,
        "age_int" = 1, 
        "ax" = if_else(.data$age == 0, lt_rule_1a0_ak(M0  = .data$mx, 
                                                      Sex = "m"), 0.5),
        "qx" = lt_id_ma_q(nMx = .data$mx, nax = .data$ax, AgeInt = .data$age_int)
      ) %>% 
      arrange(.data$cohort, .data$age) %>% 
      mutate(
        "lx" = lt_id_q_l(nqx   = .data$qx, 
                         radix = radix),
        "dx" = .data$lx * .data$qx,
        "Lx" = .data$lx - .data$dx * (1 - .data$ax),
        "SxRev" = .data$Lx / lead(.data$Lx, default = 1), .by = "cohort"
      ) %>% 
      arrange(.data$year, .data$age) %>% 
      mutate(
        "lxp" = lt_id_q_l(nqx   = .data$qx, 
                          radix = radix), 
        "dxp" = .data$lxp * .data$qx,
        "Lxp" = .data$lxp - .data$dxp * (1 - .data$ax),
        "SxRev" = ifelse(.data$year == max(.data$year), .data$Lxp / lead(.data$Lxp), .data$SxRev), 
        .by = "year") %>% 
      dplyr::filter(.data$age < 100) %>% 
      select("cohort", "year", "age", "SxRev") %>% 
      arrange(.data$cohort, -c(.data$age)) %>% 
      # inflation for reverse-surviving
      mutate("inflation_factor" = cumprod(.data$SxRev), .by = "cohort")
    
  }
  
  # if not provided by user
  if(is.null(nLxFemale)) {
    
    mxF <- mx1dt  %>% 
      dplyr::filter(.data$country_code == !!country_code,
             between(.data$year, refDate_start, refDate))  %>% 
      select("year", "age", "mx" = "mxF")  %>% 
      # could try to warp to PC shape here,
      # but uncertain infants. Maybe using
      # an a0 assumption it'd be doable.
      # need cohorts to structure reverse survival
      mutate(
        "cohort" = .data$year - .data$age - 1,
        "age_int" = 1,
        "ax" = if_else(.data$age == 0, lt_rule_1a0_ak(M0  = .data$mx, 
                                                      Sex = "f"), 0.5),
        "qx" = lt_id_ma_q(nMx = .data$mx, nax = .data$ax, AgeInt = .data$age_int)
      )  %>% 
      arrange(.data$cohort, .data$age)  %>% 
      mutate(
        "lx" = lt_id_q_l(nqx   = .data$qx, 
                         radix = radix),
        "dx" = .data$lx * .data$qx,
        "Lx" = .data$lx - .data$dx * (1 - .data$ax),
        "SxRev" = .data$Lx / lead(.data$Lx, default = 1), .by = "cohort"
      )  %>% 
      arrange(.data$year, .data$age)  %>% 
      mutate(
        "lxp" = lt_id_q_l(nqx   = .data$qx, 
                          radix = radix), 
        "dxp" = .data$lxp * .data$qx,
        "Lxp" = .data$lxp - .data$dxp * (1 - .data$ax),
        "SxRev" = ifelse(.data$year == max(.data$year), .data$Lxp / lead(.data$Lxp), .data$SxRev), 
        .by = "year") %>% 
      dplyr::filter(.data$age < 100)  %>% 
      select("cohort", "year", "age", "SxRev")  %>% 
      arrange(.data$cohort, -c(.data$age))  %>% 
      # inflation for reverse-surviving
      mutate("inflation_factor" = cumprod(.data$SxRev), .by = "cohort")
  }
  
  # --- Compute exposures ---------------------------------------------------
  # Exposures are estimated as mean of adjacent populations (mid-year approx).
  # Restricted to ages 15–49 for fertility-related computations.
  
  # (Exposure calculations for expF and expM unchanged)
  
  expF <- Females_single  %>% 
    select("cohort", "pop") %>% 
    right_join(mxF, by = c("cohort")) %>% 
    mutate("pop_hat" = .data$pop * .data$inflation_factor) %>%
    select("year", "age", "pop" = "pop_hat") %>%
    bind_rows(select(Females_single, "year", "age", "pop")) %>%
    dplyr::filter(between(.data$age, 10, 55)) %>% # are we sure about this age range??
    arrange(.data$age, .data$year) %>%
    mutate("pop_1p1"  = lead(.data$pop),
           "exposure" = (.data$pop + .data$pop_1p1) / 2) %>%
    dplyr::filter(.data$age < 55)
  
  # --- Fertility rates (ASFR) ----------------------------------------------
  # If user supplies AsfrMat, reshape to long format; otherwise, download from WPP.
  if(!is.null(AsfrMat)) {
    
    AsfrMat <- downloadASFR(Asfrmat = AsfrMat) %>% 
      as.data.frame() %>%
      # only if age is not present in data
      { if (!"age" %in% names(.)) rownames_to_column(., "age") else . } %>%
      mutate(age = readr::parse_number(age)) %>% 
      pivot_longer(-any_of(c("country_code", "name", "age")),
                   names_to  = "year",
                   values_to = "asfr") %>%
      mutate("year" = as.numeric(.data$year))
    
  }
  
  if(is.null(AsfrMat)) {
    
    AsfrMat <- downloadASFR(Asfrmat     = NULL,
                            location    = country_code,
                            AsfrDatesIn = refDate_start:refDate,
                            Age         = NULL,
                            method      = "linear",
                            output      = "single") %>%
      as.data.frame() %>% 
      rownames_to_column("age") %>%
      mutate("age" = parse_number(.data$age)) %>% 
      pivot_longer(-c("age"),
                   names_to  = "year",
                   values_to = "asfr",
                   names_transform = list("year" = as.numeric))
    
  }
  
  # --- Birth computations --------------------------------------------------
  # Estimate annual births from ASFR and exposures, then compute age since birth.
  Bt <- left_join(expF, AsfrMat , by = c("year", "age")) %>% 
    dplyr::filter(.data$year < refDate + 1) %>% 
    mutate("Bx" = .data$asfr * .data$exposure)  %>% 
    summarize(B = sum(Bx), .by = "year") %>% 
    mutate("age"  = refDate - .data$year)
  
  # --- Sex ratio at birth (SRB) --------------------------------------------
  # If missing, download from WPP; otherwise, align with cohort years.
  if(!is.null(SRB)) {
    
    SRB <- downloadSRB(SRB = SRB) %>%
      as_tibble_col(column_name = "SRB") %>%
      mutate("cohort" = refDate_start:refDate)  
  }
  
  if(is.null(SRB)) {
    
    SRB <- downloadSRB(SRB      = NULL,
                       location = country_code,
                       DatesOut = refDate_start:refDate,
                       verbose  = TRUE) %>%
      as_tibble() %>%
      mutate("cohort" = refDate_start:refDate) %>%
      rename("SRB" = "value")
    
  }
  
  # Separate births into M and F
  Bt <- Bt %>% 
    left_join(SRB, by = c("year" = "cohort")) %>% 
    mutate("Bm" = .data$B * .data$SRB  / (1 + .data$SRB),
           "Bf" = .data$B - .data$Bm)
  
  # --- Reverse-survival of infants -----------------------------------------
  # Use infant mortality (Lx values) to reverse-survive births and estimate
  # age 0–9 populations for females and males.
  # final part
  pop_hat_f <- mx1dt %>% 
    dplyr::filter(.data$country_code == !!country_code) %>% 
    select("year", "age", "mx" = "mxF") %>% 
    # could try to warp to PC shape here,
    # but uncertain infants. Maybe using
    # an a0 assumption it'd be doable.
    # need cohorts to structure reverse survival
    mutate("cohort" = .data$year - .data$age - 1,
           "age_int" = 1,
           "ax" = ifelse(.data$age == 0,
                         lt_rule_1a0_ak(M0  = .data$mx, 
                                        Sex = "f"),
                         0.5),
           "qx" = lt_id_ma_q(nMx    = .data$mx, 
                             nax    = .data$ax, 
                             AgeInt = .data$age_int)) %>% 
    dplyr::filter(between(.data$cohort, refDate_start, refDate),
           between(.data$year,   refDate_start, (refDate + 1)),
           .data$age < 10) %>% 
    arrange(.data$cohort, .data$age) %>% 
    mutate("lx" = lt_id_q_l(nqx   = .data$qx, 
                            radix = radix),
           "dx" = .data$lx * .data$qx,
           "Lx" = .data$lx - .data$dx * (1 - .data$ax), .by = "cohort") %>% 
    dplyr::filter(.data$year == max(.data$year), .by = "cohort") %>%  # ??????
    select("age", "Lx", "cohort") %>% 
    left_join(Bt, by = c("age")) %>% 
    left_join(SRB, by = c("cohort")) %>% 
    mutate("pop_hat_f" = .data$Lx * .data$Bf) %>%
    select("age", "pop_hat_f", "cohort")
  
  # for males
  pop_hat_m <- mx1dt %>% 
    dplyr::filter(.data$country_code == !!country_code) %>% 
    select("year", "age", "mx" = "mxM") %>% 
    # could try to warp to PC shape here,
    # but uncertain infants. Maybe using
    # an a0 assumption it'd be doable.
    # need cohorts to structure reverse survival
    mutate("cohort" = .data$year - .data$age - 1,
           "age_int" = 1,
           "ax" = ifelse(.data$age == 0,
                         lt_rule_1a0_ak(M0  = .data$mx, 
                                        Sex = "m"),
                         0.5),
           "qx" = lt_id_ma_q(nMx    = .data$mx, 
                             nax    = .data$ax, 
                             AgeInt = .data$age_int)) %>% 
    dplyr::filter(between(.data$cohort, refDate_start, refDate),
           between(.data$year,   refDate_start, (refDate + 1)),
           .data$age < 10) %>% 
    arrange(.data$cohort, .data$age) %>% 
    mutate("lx" = lt_id_q_l(nqx   = .data$qx, 
                            radix = radix),
           "dx" = .data$lx * .data$qx,
           "Lx" = .data$lx - .data$dx * (1 - .data$ax), .by = "cohort") %>% 
    dplyr::filter(.data$year == max(.data$year), .by = "cohort") %>%  # ??????
    select("age", "Lx", "cohort") %>% 
    left_join(Bt, by = c("age")) %>% 
    left_join(SRB, by = c("cohort")) %>% 
    mutate("pop_hat_m" = .data$Lx * .data$Bm) %>%
    select("age", "pop_hat_m", "cohort")
  
  return(list(pop_hat_m = pop_hat_m,
              pop_hat_f = pop_hat_f,
              SRB       = SRB,
              AsfrMat   = AsfrMat,
              Bt        = Bt,
              nLxFemale = nLxFemale,
              nLxMale   = nLxMale,
              expF      = expF
  ))
  
}


# #' @rdname basepop_five
# #' @aliases basepop_five
# #' @param Females_single A named numeric vector. Reported population by 1-year age # groups for \code{refDate} for females. The names of the vector should reflect the age # groups. See examples. The method assumes that the last age group is open (for example# , the population ends at '80+' and '100+')
# #' @param Males_single A named numeric vector. Reported population by 1-year age # groups format \code{refDate} for males. The names of the vector should reflect the # age groups. See examples. The method assumes that the last age group is open (for # example, the population ends at '80+' and '100+')
# #'
# #' @export
# #'
# basepop_single <- function(location = NULL,
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
#       location = location,
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

# # TR: modified to assume males and females always given
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
# 
# 
# lt_infer_radix_from_1L0 <- function(L0){
# 
#   if (L0 > 1){
#     radix_check <- L0 %>% as.integer() %>% log10()
#     is_it_a_radix <- (radix_check - round(radix_check)) == 0
# 
#     if (!is_it_a_radix){
#       pow <- L0 %>% round() %>% as.integer() %>% nchar()
# 
#       the_radix <- 10^pow
#     } else {
#       the_radix <- L0
#     }
#   } else {
#     the_radix <- 1
#   }
#   the_radix
# }
