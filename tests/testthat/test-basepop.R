context("test-basepop")
test_that("Setup: prepare Brazil population and input age groups", {
  
  data("popF", package = "wpp2024")
  data("popM", package = "wpp2024")
  
  location <- "Brazil"
  refDate  <- 1985
  
  Age_raw <- parse_number(unique(popF$age))
  
  pop_female_counts <- popF %>%
    dplyr::filter(name == location) %>%
    dplyr::pull("1985")
  
  pop_male_counts <- popM %>%
    dplyr::filter(name == location) %>%
    dplyr::pull("1985")
  
  names(pop_female_counts) <- Age_raw
  names(pop_male_counts)   <- Age_raw
  
  Females_five <<- graduate_pclm(pop_female_counts, Age_raw) %>% single2abridged()
  Males_five   <<- graduate_pclm(pop_male_counts, Age_raw) %>% single2abridged()
  Age_five     <<- names2age(Males_five)
  
  refDate_global <<- refDate
  location_global <<- location
  
  expect_true(all(Age_five >= 0))
  expect_true(length(Females_five) == length(Males_five))
})

# ---------------------------------------------------------
# 1. Minimal automatic use-case: only Females/Males_five
# ---------------------------------------------------------
test_that("basepop_five runs with only population inputs", {
  
  res <- basepop_five(
    location     = location_global,
    refDate      = refDate_global,
    Females_five = Females_five,
    Males_five   = Males_five,
    Age          = Age_five,
    verbose      = FALSE
  )
  
  expect_type(res, "list")
  expect_true("Females_adjusted" %in% names(res))
  expect_true("Males_adjusted" %in% names(res))
  
  # anything beyond age 10 is not adjusted
  expect_equal(
    res$Males_adjusted[4:length(res$Males_adjusted)], 
    res$Males_five[4:length(res$Males_five)])
  
})

# ---------------------------------------------------------
# 2. Missing Age but radix provided
# ---------------------------------------------------------
test_that("basepop_five works when Age = NULL but radix is provided", {
  
  res <- basepop_five(
    location     = location_global,
    refDate      = refDate_global,
    Females_five = Females_five,
    Males_five   = Males_five,
    Age          = NULL,
    radix        = 1,
    verbose      = FALSE
  )
  
  expect_type(res, "list")
  
})

# ---------------------------------------------------------
# 3. Provided nLxFemale only
# ---------------------------------------------------------
test_that("basepop_five works with supplied nLxFemale only", {
  
  nLxFemale <- downloadnLx(
    nLx         = NULL,
    Age         = Age_five,
    location    = location_global,
    nLxDatesIn  = refDate_global - c(0.5, 7.5),
    gender      = "female",
    output      = "abridged",
    radix       = 1,
    verbose     = FALSE
  )
  
  res <- basepop_five(
    location     = location_global,
    refDate      = refDate_global,
    Females_five = Females_five,
    Males_five   = Males_five,
    nLxFemale    = nLxFemale,
    nLxMale      = NULL,
    Age          = Age_five,
    verbose      = FALSE,
    radix        = 1
  )
  
  expect_type(res, "list")
})

# ---------------------------------------------------------
# 4. Provided nLxFemale AND nLxMale
# ---------------------------------------------------------
test_that("basepop_five works with supplied nLxFemale and nLxMale", {
  
  nLxFemale <- downloadnLx(
    nLx         = NULL,
    Age         = Age_five,
    location    = location_global,
    nLxDatesIn  = refDate_global - c(0.5, 7.5),
    gender      = "female",
    output      = "abridged",
    radix       = 1,
    verbose     = FALSE
  )
  
  nLxMale <- downloadnLx(
    nLx         = NULL,
    Age         = Age_five,
    location    = location_global,
    nLxDatesIn  = refDate_global - c(0.5, 7.5),
    gender      = "male",
    output      = "abridged",
    radix       = 1,
    verbose     = FALSE
  )
  
  res <- basepop_five(
    location     = location_global,
    refDate      = refDate_global,
    Females_five = Females_five,
    Males_five   = Males_five,
    nLxFemale    = nLxFemale,
    nLxMale      = nLxMale,
    Age          = Age_five,
    verbose      = FALSE
  )
  
  expect_type(res, "list")
})

# ---------------------------------------------------------
# 5. Provided nLxFemale + nLxMale + nLxDatesIn
# ---------------------------------------------------------
test_that("basepop_five works with nLxFemale, nLxMale and nLxDatesIn", {
  
  nLxFemale <- downloadnLx(
    nLx         = NULL,
    Age         = Age_five,
    location    = location_global,
    nLxDatesIn  = refDate_global - c(0.5, 7.5),
    gender      = "female",
    output      = "abridged",
    radix       = 1,
    verbose     = FALSE
  )
  
  nLxMale <- downloadnLx(
    nLx         = NULL,
    Age         = Age_five,
    location    = location_global,
    nLxDatesIn  = refDate_global - c(0.5, 7.5),
    gender      = "male",
    output      = "abridged",
    radix       = 1,
    verbose     = FALSE
  )
  
  res <- basepop_five(
    location     = location_global,
    refDate      = refDate_global,
    Females_five = Females_five,
    Males_five   = Males_five,
    nLxFemale    = nLxFemale,
    nLxMale      = nLxMale,
    nLxDatesIn   = refDate_global - c(0.5, 7.5),
    Age          = Age_five,
    verbose      = FALSE
  )
  
  expect_type(res, "list")
})

# ---------------------------------------------------------
# 6. Supplied nLxFemale + nLxMale but Age = NULL
# ---------------------------------------------------------
test_that("basepop_five works when nLx provided but Age = NULL", {
  
  nLxFemale <- downloadnLx(
    nLx         = NULL,
    Age         = Age_five,
    location    = location_global,
    nLxDatesIn  = refDate_global - c(0.5, 7.5),
    gender      = "female",
    output      = "abridged",
    radix       = 1,
    verbose     = FALSE
  )
  
  nLxMale <- downloadnLx(
    nLx         = NULL,
    Age         = Age_five,
    location    = location_global,
    nLxDatesIn  = c(1977.31, 1986.50),
    gender      = "male",
    output      = "abridged",
    radix       = 1,
    verbose     = FALSE
  )
  
  res <- basepop_five(
    location     = location_global,
    refDate      = refDate_global,
    Females_five = Females_five,
    Males_five   = Males_five,
    nLxFemale    = nLxFemale,
    nLxMale      = nLxMale,
    Age          = NULL,
    verbose      = FALSE
  )
  
  expect_type(res, "list")
})

# ---------------------------------------------------------
# 7. Provided ASFR matrix
# ---------------------------------------------------------
test_that("basepop_five works with supplied ASFR matrix", {
  
  AsfrMat <- downloadAsfr(
    Asfrmat     = NULL,
    location    = location_global,
    AsfrDatesIn = refDate_global - c(0.5, 7.5),
    Age         = NULL,
    output      = "5-year",
    verbose     = FALSE
  )
  
  nLxFemale <- downloadnLx(
    nLx         = NULL,
    Age         = Age_five,
    location    = location_global,
    nLxDatesIn  = refDate_global - c(0.5, 7.5),
    gender      = "female",
    output      = "abridged",
    radix       = 1,
    verbose     = FALSE
  )
  
  nLxMale <- downloadnLx(
    nLx         = NULL,
    Age         = Age_five,
    location    = location_global,
    nLxDatesIn  = refDate_global - c(0.5, 7.5),
    gender      = "male",
    output      = "abridged",
    radix       = 1,
    verbose     = FALSE
  )
  
  res <- basepop_five(
    location     = location_global,
    refDate      = refDate_global,
    Females_five = Females_five,
    Males_five   = Males_five,
    nLxFemale    = nLxFemale,
    nLxMale      = nLxMale,
    AsfrMat      = AsfrMat,
    Age          = NULL,
    verbose      = FALSE
  )
  
  expect_type(res, "list")
})

# ---------------------------------------------------------
# 8. Provided SRB
# ---------------------------------------------------------
test_that("basepop_five works with supplied SRB", {
  
  SRB <- downloadSRB(
    SRB      = NULL,
    location = location_global,
    DatesOut = refDate_global - c(0.5, 7.5),
    verbose  = FALSE
  )
  
  res <- basepop_five(
    location     = location_global,
    refDate      = refDate_global,
    Females_five = Females_five,
    Males_five   = Males_five,
    SRB          = SRB,
    Age          = Age_five,
    verbose      = FALSE
  )
  
  expect_type(res, "list")
})


# ---------------------------------------------------------------
# Helper population to use for error tests
# ---------------------------------------------------------------

data("popF", package = "wpp2024")
data("popM", package = "wpp2024")

location <- "Brazil"
refDate  <- 1985

Age_raw <- parse_number(unique(popF$age))

pop_female_counts <- dplyr::filter(popF, name == location)$"1985"
pop_male_counts   <- dplyr::filter(popM, name == location)$"1985"

names(pop_female_counts) <- Age_raw
names(pop_male_counts)   <- Age_raw

Females_five <- single2abridged(graduate_pclm(pop_female_counts, Age_raw)) 
Males_five   <- single2abridged(graduate_pclm(pop_male_counts, Age_raw))
Age_five     <- names2age(Females_five)

# ---------------------------------------------------------------
# ERROR TEST 1 — missing population
# ---------------------------------------------------------------
test_that("Error when Females_five or Males_five missing", {
  expect_error(
    basepop_five(
      location     = location,
      refDate      = refDate,
      Females_five = NULL,
      Males_five   = Males_five
    ),
    "Male or female population was not provided"
  )
  
  expect_error(
    basepop_five(
      location     = location,
      refDate      = refDate,
      Females_five = Females_five,
      Males_five   = NULL
    ),
    "Male or female population was not provided"
  )
})

# ---------------------------------------------------------------
# ERROR TEST 2 — Age provided but not abridged
# ---------------------------------------------------------------
test_that("Error when Age provided but not abridged", {
  bad_age <- 0:30  # not abridged
  
  expect_error(
    basepop_five(
      location     = location,
      refDate      = refDate,
      Age          = bad_age,
      Females_five = Females_five,
      Males_five   = Males_five
    ),
    regexp = "is_abridged"
  )
})

# ---------------------------------------------------------------
# ERROR TEST 3 — Age length mismatches population length
# ---------------------------------------------------------------
test_that("Error when Age and Females_five lengths differ", {
  short_age <- Age_five[-1]  # delete element
  
  expect_error(
    basepop_five(
      location     = location,
      refDate      = refDate,
      Age          = short_age,
      Females_five = Females_five,
      Males_five   = Males_five
    )
  )
})

# ---------------------------------------------------------------
# ERROR TEST 4 — radix does not match downloaded nLx radix
# ---------------------------------------------------------------
test_that("Error when user-supplied radix mismatches implied nLx radix", {
  
  # download real nLxFemale
  nLxFemale <- downloadnLx(
    nLx        = NULL,
    Age        = Age_five,
    location   = location,
    gender     = "female",
    nLxDatesIn = refDate - c(0.5, 7.5),
    output     = "abridged",
    radix      = 1,
    verbose    = FALSE
  )
  
  # trick: give wrong radix
  expect_error(
    basepop_five(
      location     = location,
      refDate      = refDate,
      Females_five = Females_five,
      Males_five   = Males_five,
      Age          = Age_five,
      radix        = 999,    # impossible
      nLxFemale    = nLxFemale,
      nLxMale      = nLxFemale
    ),
    regexp = "Your provided radix value do not match"
  )
})

# ---------------------------------------------------------------
# ERROR TEST 5 — nLxDatesIn requires >5y extrapolation
# ---------------------------------------------------------------
test_that("Error when nLxDatesIn would require >5 years extrapolation", {
  expect_error(
    basepop_five(
      location     = location,
      refDate      = refDate,
      Females_five = Females_five,
      Males_five   = Males_five,
      Age          = Age_five,
      nLxDatesIn   = c(1800, 1805),  # far away → extrapolation > 5
      verbose      = FALSE
    ),
    regexp = "extrapolation of > 5 years"
  )
})

# ---------------------------------------------------------------
# ERROR TEST 6 — AsfrDatesIn requires >5y extrapolation
# ---------------------------------------------------------------
test_that("Error when AsfrDatesIn would require >5 years extrapolation", {
  expect_error(
    basepop_five(
      location     = location,
      refDate      = refDate,
      Females_five = Females_five,
      Males_five   = Males_five,
      Age          = Age_five,
      AsfrDatesIn  = c(1800, 1805),  # also far away
      verbose      = FALSE
    ),
    regexp = "AsfrDatesIn implies an extrapolation of > 5 years"
  )
})

# ---------------------------------------------------------------
# ERROR TEST 7 — Negative population counts
# ---------------------------------------------------------------
test_that("Error when population contains negative values", {
  Fem_neg <- Females_five
  Fem_neg[3] <- -10
  
  expect_error(
    basepop_five(
      location     = location,
      refDate      = refDate,
      Females_five = Fem_neg,
      Males_five   = Males_five,
      Age          = Age_five
    )
  )
})

# ---------------------------------------------------------------
# ERROR TEST 8 — mismatch Females_five and Males_five lengths
# ---------------------------------------------------------------
test_that("Error when male/female population differ in length", {
  Fem_short <- Females_five[-1]
  
  expect_error(
    basepop_five(
      location     = location,
      refDate      = refDate,
      Females_five = Fem_short,
      Males_five   = Males_five,
      Age          = Age_five[-1]
    )
  )
})


# BASEPOP SINGLE ----------------------------------- 
refDate  <- 2000
location <- "Germany"
data("popF",  package = "wpp2024")
data("popM",  package = "wpp2024")

Age <- parse_number(unique(popF$age))

pop_female_counts <- popF %>%
  dplyr::filter(name == location) %>% 
  dplyr::pull("1985")

pop_male_counts   <- popM %>%
  dplyr::filter(name == location) %>% 
  dplyr::pull("1985")

names(pop_female_counts) <- names(pop_female_counts) <-  Age

# graduate
Females <- graduate_pclm(pop_female_counts, Age = Age)
Males   <- graduate_pclm(pop_male_counts,  Age = Age)
Age     <- 0:100


# test 1 - errors when insufficinet data is provided
test_that("Error when country_code is not numeric or ref date is not provided or more than 1 reference date is requires", {
  expect_error(basepop_single(country_code = "Germany",
                              refDate      = 1985))
  expect_error(basepop_single(country_code = "276",
                              refDate      = 1985))
  expect_error(basepop_single(country_code = 276,
                              refDate      = NULL))
  expect_error(
    basepop_single(
      country_code = 276,
      refDate = c(2000, 2001)
    ))
  
})

# test 2 - function works with minimal data required
test_that("Minimal data setting", {
  # germany
  res <- basepop_single(
    country_code = 276,
    refDate = 2000
  )
  
  expect_type(res, "list")
  expect_true(nrow(res$pop_hat_m) == 10)
  
})

# test 3 - radix is provided
test_that("radix setup and more than one refDate", {
  # germany
  res <- basepop_single(
    country_code = 276,
    refDate = 2000,
    radix = 1
  )
  expect_true(nrow(res$pop_hat_m) == 10)
  
})


# test 4 - Females and Males single are provided
test_that("Females and Males single are provided AND age is provided OR age is wrong", {
  # germany
  res <- basepop_single(
    Females_single = Females,
    Males_single   = Males,
    country_code   = 276,
    refDate        = 2000,
    radix          = 1
  )
  expect_true(nrow(res$pop_hat_m) == 10)
  expect_true(nrow(na.omit(res$pop_hat_m)) == 10)
  
  
  res1 <- basepop_single(
    Females_single = Females,
    Males_single   = Males,
    country_code   = 276,
    refDate        = 2000,
    radix          = 1,
    Age            = 0:100
  )
  
  expect_true(nrow(res1$pop_hat_m) == 10)
  expect_true(nrow(na.omit(res1$pop_hat_m)) == 10)
  
  # wrong age shold return error
  expect_error(basepop_single(
    Females_single = Females,
    Males_single   = Males,
    country_code   = 276,
    refDate        = 2000,
    radix          = 1,
    Age            = 0:80
  ))
  
})

# test 5 - nLxFemale and nLxMale AsfrMat, SRB are provided one by one
test_that("nLxFemale and nLxMale AsfrMat, SRB are provided", {
  
  # germany
  nLxFemale <- downloadnLx(location = 276,
                           gender = "f",
                           nLxDatesIn = 1990:1999,
                           output = "single",
                           radix = 1)
  nLxMale <- downloadnLx(location = 276,
                         gender = "m",
                         nLxDatesIn = 1990:1999,
                         output = "single",
                         radix = 1)
  AsfrMat <- downloadAsfr(location = 276,
                          AsfrDatesIn = 1990:1999,
                          output = "single")
  SRB <- downloadSRB(location = 276,
                     DatesOut = 1990:1999)
  
  # start with only nLx available
  res <- basepop_single(
    Females_single = Females,
    Males_single   = Males,
    country_code   = 276,
    refDate        = 2000,
    radix          = 1,
    nLxFemale      = nLxFemale,
    nLxMale        = nLxMale
  )
  
  expect_true(nrow(na.omit(res$pop_hat_m)) == 10)
  expect_true(nrow(res$pop_hat_m) == 10)
  
  # now also ASFR and Age
  res1 <- basepop_single(
    Females_single = Females,
    Males_single   = Males,
    country_code   = 276,
    refDate        = 2000,
    radix          = 1,
    nLxFemale      = nLxFemale,
    nLxMale        = nLxMale,
    AsfrMat        = AsfrMat,
    Age            = 0:100
  )
  
  expect_true(nrow(na.omit(res1$pop_hat_m)) == 10)
  expect_true(nrow(res1$pop_hat_m) == 10)
  
  # now also SRB
  res2 <- basepop_single(
    Females_single = Females,
    Males_single   = Males,
    country_code   = 276,
    refDate        = 2000,
    radix          = 1,
    nLxFemale      = nLxFemale,
    nLxMale        = nLxMale,
    AsfrMat        = AsfrMat,
    Age            = 0:100,
    SRB            = SRB
  )
  
  expect_true(nrow(na.omit(res2$pop_hat_m)) == 10)
  expect_true(nrow(res2$pop_hat_m) == 10)
  
  # now everything
  res3 <- basepop_single(
    Females_single = Females,
    Males_single   = Males,
    country_code   = 276,
    refDate        = 2000,
    radix          = 1,
    nLxFemale      = nLxFemale,
    nLxMale        = nLxMale,
    AsfrMat        = AsfrMat,
    Age            = 0:100,
    SRB            = SRB,
    nLxDatesIn     = 1990:1999,
    AsfrDatesIn    = 1990:1999,
    SRBDatesIn     = 1990:1999
  )
  
  expect_true(nrow(na.omit(res3$pop_hat_m)) == 10)
  expect_true(nrow(res3$pop_hat_m) == 10)
  
})



# # (1) RefYear
# refDate <- 1986.21
# library(testthat)
# # (2) Reported population by 5-year age groups and sex in the base year (Include unknowns).
# pop_male_counts <- c(`0` = 11684,`1` = 46738,`5` = 55639,`10` = 37514,
#                      `15` = 29398, `20` = 27187,`25` = 27770,`30` = 20920,
#                      `35` = 16973,`40` = 14999, `45` = 11330,`50` = 10415,
#                      `55` = 6164,`60` = 7330,`65` = 3882, `70` = 3882,
#                      `75` = 1840,`80` = 4200
#                      )
# 
# pop_female_counts <- c(`0` = 11673,`1` = 46693,`5` = 55812,`10` = 35268,
#                        `15` = 33672, `20` = 31352,`25` = 33038,`30` = 24029,
#                        `35` = 16120,`40` = 14679, `45` = 8831,`50` = 9289,
#                        `55` = 4172,`60` = 6174,`65` = 2715, `70` = 3344,
#                        `75` = 1455, `80` = 4143)
# 
# # (4) Sex ratio at birth (m/f)
# sex_ratio <- 1.0300
# 
# # (6) The male and female nLx functions for ages under 1 year, 1 to 4 years, and 5 to 9
# # years, pertaining to an earlier and later date
# nLxDatesIn <- c(1977.31, 1986.50)
# 
# nLxMale <- matrix(c(87732,
#                     304435,
#                     361064,
#                     88451,
#                     310605,
#                     370362
#                     ),
#                   nrow = 3, ncol = 2)
# # includes age 10 patch
# nLxFemale <- matrix(c(89842,314521,372681,666666,353053,340650,326588,311481,
#                       295396,278646,261260,241395,217419,90478,320755,
#                       382531,666666,364776, 353538,340687, 326701, 311573,
#                       295501, 278494, 258748,234587),
#                     nrow = 13,
#                     ncol = 2)
# rownames(nLxFemale) <- c(0,1,seq(5,55,by=5))
# # (7) A set of age-specific fertility rates pertaining to an earlier and later
# # date
# 
# AsfrMat <- matrix(c(0.2000,0.3000,0.3000, 0.2500, 0.2000,
#                     0.1500, 0.0500,0.1500,0.2000,0.2750,
#                     0.2250, 0.1750, 0.1250, 0.0500),
#                   nrow = 7,
#                   ncol = 2)
# rownames(AsfrMat) <- seq(15,45,by=5)
# AsfrDatesIn <- c(1977.81, 1985.71)
# 
# smoothed_females <- smooth_age_5(Value = pop_female_counts,
#                                  Age = as.numeric(names(pop_female_counts)),
#                                  method = "Arriaga",
#                                  OAG = TRUE,
#                                  young.tail = "Original")
# smoothed_females <- c(pop_female_counts[1:2], smoothed_females[-1])
# 
# smoothed_males <- smooth_age_5(Value = pop_male_counts,
#                                  Age = as.numeric(names(pop_male_counts)),
#                                  method = "Arriaga",
#                                  OAG = TRUE,
#                                  young.tail = "Original")
# smoothed_males <- c(pop_male_counts[1:2], smoothed_males[-1])
# 
# ## This is the only number that messes up the whole calculation.
# ## smooth_age_5 returns the same result as the PASS excel sheet
# ## except for the age groups 10-15 and 15-19. Here we only use
# ## age group 15-19. If we plug in manually the correct value,
# ## we get all results match exactly, otherwise there are
# ## some differences.
# smoothed_females[4] <- 34721
# # smoothed_males will be off in the same cell, but inconsequential for children
# test_that("basepop_five - bpa matches the expected result from PASS", {
# 
#   bpa <-
#     basepop_five(
#       refDate = 1986.21,
#       Males_five = smoothed_males,
#       Females_five = smoothed_females,
#       SRB = sex_ratio,
#       Age = NULL,
#       nLxFemale = nLxFemale,
#       nLxMale    = nLxMale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       method = "linear",
#       radix = 100000
#     )
# 
#   # matches BASEPOP1altA.xlsx
#   #expect_true(all(round(bpa$Males_adjusted[1:3], 0) == c(13559, 47444, 54397)))
#   expect_true(all(round(bpa$Males_adjusted[1:3], 0) == c(13406, 47184, 54397)))
#   #expect_true(all(round(bpa_female[1:3], 0) == c(13467, 47576, 54554)))
#   expect_true(all(round(bpa$Females_adjusted[1:3], 0) == c(13315, 47315, 54554)))
# })
# 
# test_that("basepop_five - bpe matches the expected result from PASS", {
# 
#   bpe <-
#     basepop_five(
#       refDate = refDate,
#       Males_five = pop_male_counts,
#       Females_five = pop_female_counts,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxMale = nLxMale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       method = "linear",
#       radix = 100000
#     )
# 
#   expect_true(all(round(bpe$Males_adjusted[1:3], 0) == c(13679, 47967, 55721)))
#   expect_true(all(round(bpe$Females_adjusted[1:3], 0) == c(13587, 48101, 55882)))
# })
# 
# 
# test_that("basepop_single fails if provided five year age groups", {
#
#   female_single <-
#   c(
#     `0` = 11673,
#     `1` = 11474,
#     `2` = 11670,
#     `3` = 11934,
#     `4` = 11614,
#     `5` = 10603,
#     `6` = 11144,
#     `7` = 11179,
#     `8` = 11269,
#     `9` = 11617,
#     `10` = 6772,
#     `11` = 6948,
#     `12` = 7030,
#     `13` = 7211,
#     `14` = 7306,
#     `15` = 6531,
#     `16` = 6443,
#     `17` = 6535,
#     `18` = 6951,
#     `19` = 7213,
#     `20` = 6096,
#     `21` = 6234,
#     `22` = 6327,
#     `23` = 6410,
#     `24` = 6285,
#     `25` = 6464,
#     `26` = 6492,
#     `27` = 6549,
#     `28` = 6739,
#     `29` = 6795,
#     `30` = 5013,
#     `31` = 4888,
#     `32` = 4735,
#     `33` = 4747,
#     `34` = 4646,
#     `35` = 3040,
#     `36` = 3068,
#     `37` = 3107,
#     `38` = 3246,
#     `39` = 3658,
#     `40` = 2650,
#     `41` = 2788,
#     `42` = 2977,
#     `43` = 3108,
#     `44` = 3156,
#     `45` = 1756,
#     `46` = 1784,
#     `47` = 1802,
#     `48` = 1764,
#     `49` = 1724,
#     `50` = 1982,
#     `51` = 1935,
#     `52` = 1846,
#     `53` = 1795,
#     `54` = 1731,
#     `55` = 863,
#     `56` = 850,
#     `57` = 825,
#     `58` = 819,
#     `59` = 816,
#     `60` = 1348,
#     `61` = 1342,
#     `62` = 1246,
#     `63` = 1138,
#     `64` = 1101,
#     `65` = 391,
#     `66` = 520,
#     `67` = 585,
#     `68` = 560,
#     `69` = 659,
#     `70` = 670,
#     `71` = 750,
#     `72` = 686,
#     `73` = 634,
#     `74` = 604,
#     `75` = 353,
#     `76` = 340,
#     `77` = 270,
#     `78` = 246,
#     `79` = 247,
#     `80` = 4143
#   )
#
#   # Correction for males
#   # To test that males are checked for single ages, we first
#   # correctly define the female pop counts as single ages
#   expect_error(
#     basepop_single(
#       refDate = refDate,
#       Males_single = pop_male_counts,
#       Females_single = female_single,
#       SmoothedFemales = smoothed_females,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxMale = nLxMale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       method = "linear"
#     ),
#     "is_single(as.numeric(names(Males_single))) is not TRUE",
#     fixed = TRUE
#   )
#
#   # check that pop male counts are named
#   expect_error(
#     basepop_single(
#       refDate = refDate,
#       Males_single = setNames(pop_male_counts, NULL),
#       Females_single = female_single,
#       SmoothedFemales = smoothed_females,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxMale = nLxMale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       method = "linear"
#     ),
#     "!is.null(names(Males_single)) is not TRUE",
#     fixed = TRUE
#   )
#
#   # Check female counts are single ages
#   expect_error(
#     basepop_single(
#       refDate = refDate,
#       Females_single = pop_female_counts,
#       SmoothedFemales = smoothed_females,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       method = "linear"
#     ),
#     "is_single(as.numeric(names(Females_single))) is not TRUE",
#     fixed = TRUE
#   )
#
#   # Check that pop female counts are named
#   expect_error(
#     basepop_single(
#       refDate = refDate,
#       Females_single = setNames(pop_female_counts, NULL),
#       SmoothedFemales = smoothed_females,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       method = "linear"
#     ),
#     "!is.null(names(Females_single)) is not TRUE",
#     fixed = TRUE
#   )
#
#   # Works as well for smoothing
#   expect_error(
#     basepop_single(
#       refDate = refDate,
#       Males_single = pop_male_counts,
#       Females_single = pop_female_counts,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxMale = nLxMale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       method = "linear"
#     ),
#     "is_single(as.numeric(names(Females_single))) is not TRUE",
#     fixed = TRUE
#   )
# })
# 
# pop_male_counts <-
#   c(`0` = 11684,`1` = 11473,`2` = 11647,`3` = 11939,`4` = 11680,
#     `5` = 10600,`6` = 11100,`7` = 11157,`8` = 11238,`9` = 11544,
#     `10` = 7216,`11` = 7407,`12` = 7461,`13` = 7656,`14` = 7774,
#     `15` = 5709,`16` = 5629,`17` = 5745,`18` = 6056,`19` = 6259,
#     `20` = 5303,`21` = 5423,`22` = 5497,`23` = 5547,`24` = 5417,
#     `25` = 5441,`26` = 5466,`27` = 5500,`28` = 5668,`29` = 5694,
#     `30` = 4365,`31` = 4252,`32` = 4122,`33` = 4142,`34` = 4039,
#     `35` = 3210,`36` = 3222,`37` = 3258,`38` = 3413,`39` = 3871,
#     `40` = 2684,`41` = 2844,`42` = 3052,`43` = 3182,`44` = 3237,
#     `45` = 2263,`46` = 2298,`47` = 2318,`48` = 2257,`49` = 2194,
#     `50` = 2231,`51` = 2172,`52` = 2072,`53` = 2008,`54` = 1932,
#     `55` = 1301,`56` = 1262,`57` = 1213,`58` = 1197,`59` = 1191,
#     `60` = 1601,`61` = 1593,`62` = 1490,`63` = 1348,`64` = 1299,
#     `65` = 568,`66` = 745,`67` = 843,`68` = 801,`69` = 925,
#     `70` = 806,`71` = 883,`72` = 796,`73` = 725,`74` = 672,
#     `75` = 470,`76` = 441,`77` = 340,`78` = 300,`79` = 289,
#     `80` = 4200
#   )
# 
# pop_female_counts <-
#   c(`0` = 11673,`1` = 11474,`2` = 11670,`3` = 11934,`4` = 11614,
#     `5` = 10603,`6` = 11144,`7` = 11179,`8` = 11269,`9` = 11617,
#     `10` = 6772,`11` = 6948,`12` = 7030,`13` = 7211,`14` = 7306,
#     `15` = 6531,`16` = 6443,`17` = 6535,`18` = 6951,`19` = 7213,
#     `20` = 6096,`21` = 6234,`22` = 6327,`23` = 6410,`24` = 6285,
#     `25` = 6464,`26` = 6492,`27` = 6549,`28` = 6739,`29` = 6795,
#     `30` = 5013,`31` = 4888,`32` = 4735,`33` = 4747,`34` = 4646,
#     `35` = 3040,`36` = 3068,`37` = 3107,`38` = 3246,`39` = 3658,
#     `40` = 2650,`41` = 2788,`42` = 2977,`43` = 3108,`44` = 3156,
#     `45` = 1756,`46` = 1784,`47` = 1802,`48` = 1764,`49` = 1724,
#     `50` = 1982,`51` = 1935,`52` = 1846,`53` = 1795,`54` = 1731,
#     `55` = 863,`56` = 850,`57` = 825,`58` = 819,`59` = 816,
#     `60` = 1348,`61` = 1342,`62` = 1246,`63` = 1138,`64` = 1101,
#     `65` = 391,`66` = 520,`67` = 585,`68` = 560,`69` = 659,
#     `70` = 670,`71` = 750,`72` = 686,`73` = 634,`74` = 604,
#     `75` = 353,`76` = 340,`77` = 270,`78` = 246,`79` = 247,
#     `80` = 4143
#   )
# 
# smoothed_females <- smooth_age_5(Value = pop_female_counts,
#                                  Age = as.numeric(names(pop_female_counts)),
#                                  method = "Arriaga",
#                                  OAG = TRUE,
#                                  young.tail = "Original")

# test_that("basepop_single does calculation for males when providing Males_single", {
#
#   # Since the default is for females, I saved the correct calculations for males
#   # and just check that the current implementation matches it
#   res <-
#     basepop_single(
#       refDate = refDate,
#       Males_single = pop_male_counts,
#       Females_single = pop_female_counts,
#       SmoothedFemales = smoothed_females,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxMale = nLxMale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       method = "linear",
#       radix = 100000
#     )
#
#   correct_res_males <-
#     c(
#       `0` = 13315,
#       `1` = 11614,
#       `2` = 11791,
#       `3` = 12086,
#       `4` = 11824,
#       `5` = 10393,
#       `6` = 10883,
#       `7` = 10939,
#       `8` = 11019,
#       `9` = 11319
#     )
#
#   expect_equivalent(round(res[1:10], 0), correct_res_males)
# })

# refDate <- 1986
# location <- "Spain"
# res <- fertestr::FetchPopWpp2019(location, refDate, ages = 0:100, sex = "female")
# pop_female_counts <- single2abridged(setNames(res$pop, res$ages))
# 
# res <- fertestr::FetchPopWpp2019(location, refDate, ages = 0:100, sex = "male")
# pop_male_counts <- single2abridged(setNames(res$pop, res$ages))
# 
# # Download asfr matrix to test that it can download the nLx only
# invisible(
#   capture.output(
#     nLxFemale <-
#       downloadnLx(
#         NULL,
#         location,
#         gender = "female",
#         c(1978, 1985.5)
#       )
#   )
# )
# 
# # Download asfr matrix to test that it can download the nLx only
# invisible(capture.output(AsfrMat <- downloadAsfr(NULL, location, c(1978, 1985.5))))
# 
# test_that("basepop_five can download data for nLx", {
# 
#   # Female
#   output <-
#     capture.output(
#       basepop_five(
#         location = location,
#         refDate = refDate,
#         Females_five = pop_female_counts,
#         Males_five = pop_male_counts,
#         AsfrMat = AsfrMat
#       )
#     )
# 
#   # Test that the output print doesn't print its
#   # downloading the Asfr
#   expect_false(any(grepl("^Downloading Asfr", output)))
# 
#   # Test that it's downloading nLx for the two years and female
#   expect_true(sum(grepl("^Downloading nLx", output)) == 2)
# 
# })
# 
# test_that("basepop_five can download data for asfr", {
# 
#   output <-
#     capture.output(
#       basepop_five(
#         location = location,
#         refDate = refDate,
#         nLxFemale = nLxFemale,
#         nLxMale = nLxMale,
#         Females_five = pop_female_counts,
#         Males_five = pop_male_counts
#       )
#     )
# 
#   # Test that the output print doesn't print it's
#   # downloading the nLx
#   expect_false(any(grepl("^Downloading nLx", output)))
# 
#   # For males
#   # There's no need to test this separately between gender because
#   # the Asfr matrix is downloaded regardless of gender. However,
#   # since the code is becoming convoluted, I'm testing just in case
# 
#   # output <-
#   #   capture.output(
#   #     basepop_five(
#   #       location = location,
#   #       refDate = refDate,
#   #       nLxFemale = nLxFemale,
#   #       Females_five = pop_female_counts,
#   #       Males_five = pop_male_counts,
#   #       female = TRUE
#   #     )
#   #   )
#   #
#   # expect_false(any(grepl("^Downloading nLx", output)))
# })
# 
# test_that("basepop_five infers radix if not provided", {
# 
#   output <-
#     capture.output(
#       basepop_five(
#         location = location,
#         refDate = refDate,
#         nLxFemale = nLxFemale,
#         nLxMale = nLxMale,
#         Females_five = pop_female_counts,
#         Males_five = pop_male_counts
#       )
#     )
# 
#   expect_true(sum(grepl("^Setting radix", output)) == 1)
# })
# 
# test_that("basepop raises error when downloads needed but no location is specified", {
# 
#   expect_error(
#     basepop_five(
#       refDate = refDate,
#       Females_five = pop_female_counts,
#       Males_five = pop_male_counts,
#       verbose = FALSE
#     ),
#     "You need to provide a location to download the data for nLx"
#   )
# 
# 
#   expect_error(
#     basepop_five(
#       refDate = refDate,
#       AsfrMat = AsfrMat,
#       Females_five = pop_female_counts,
#       Males_five = pop_male_counts,
#       radix = 1,
#       verbose = FALSE
#     ),
#     "You need to provide a location to download the data for nLx"
#   )
# 
#   expect_error(
#     basepop_five(
#       refDate = refDate,
#       nLxFemale = nLxFemale,
#       nLxMale = nLxMale,
#       Females_five = pop_female_counts,
#       Males_five = pop_male_counts,
#       radix = 1,
#       verbose = FALSE
#     ),
#     "You need to provide a location to download the data for Asfrmat"
#   )
# 
#   # If provided all correct arguments, it download the data
#   # successfully
#   expect_success({
#     res <-
#       basepop_five(
#         location = "Spain",
#         refDate = refDate,
#         AsfrMat = AsfrMat,
#         Females_five = pop_female_counts,
#         Males_five = pop_male_counts,
#         radix = 1,
#         verbose = FALSE
#       )
# 
#     expect_type(res$Females_adjusted, "double")
#   })
# 
# })
# 
# test_that("basepop_five can download from dates provided", {
# 
#   output <-
#     capture.output(
#       basepop_five(
#         location = location,
#         refDate = refDate,
#         nLxDatesIn = c(1978, 1986.5),
#         AsfrDatesIn = c(1978, 1985.5),
#         Females_five = pop_female_counts,
#         Males_five = pop_male_counts
#       )
#     )
# 
#   # Test that the output print doesn't print its
#   # calculating the datesin for both nLx and Asfr.
#   expect_false(any(grepl("^Assuming the two", output)))
# 
#   output <-
#     capture.output(
#       basepop_five(
#         location = location,
#         refDate = refDate,
#         nLxDatesIn = c(1978, 1986.5),
#         Females_five = pop_female_counts,
#         Males_five = pop_male_counts
#       )
#     )
# 
#   # Test that it calculates the dates for the asfr
#   expect_true(sum(grepl("^Assuming the two", output)) == 1)
# 
#   output <-
#     capture.output(
#       basepop_five(
#         location = location,
#         refDate = refDate,
#         Females_five = pop_female_counts,
#         Males_five = pop_male_counts
#       )
#     )
# 
#   # Tests that it calculates the dates for nLx and asfr
#   expect_true(sum(grepl("^Assuming the two", output)) == 2)
# })
# 
# test_that("basepop works with up to year 1955", {
# 
#   # This is where the test actually happens since
#   # internally we subtract 7.5 from the refDate to download
#   # the nLx and Asfr data. So the minimum year will be 1955.
#   res <-
#     basepop_five(
#       location = "Spain",
#       refDate= 1962.5,
#       Males_five = smoothed_males,
#       Females_five = smoothed_females,
#       SRB = sex_ratio,
#       method = "linear",
#       radix = 100000
#   )
# 
#   expect_true("1955" %in% colnames(res$Asfr))
#   expect_true("1955" %in% colnames(res$nLxm))
#   expect_true("1955" %in% colnames(res$nLxf))
# 
# })
# 
# test_that("basepop works well with SRBDatesIn", {
#   # Most of the tests for SRBDatesIn are actually done when
#   # testing downloadSRB. These tests just make sure that
#   # basepop can handle NULL/Non-NULL SRBDatesIn dates. Everything
#   # else is forward to downloadSRB.
# 
#   # Works when SRBDatesIn is NULL, meaning that it convert them
#   # to refDate - c(0.5, 2.5, 7.5)
# 
#   expect_success({
#     res <-
#       basepop_five(
#         location = "Spain",
#         refDate= 1962.5,
#         Males_five = smoothed_males,
#         Females_five = smoothed_females,
#         SRB = sex_ratio,
#         SRBDatesIn = NULL,
#         method = "linear",
#         radix = 100000
#       )
# 
#     expect_type(res, "list")
#   })
# 
#   # Works when SRBDatesIn is an actual date.
#   expect_success({
#     res <-
#       basepop_five(
#         location = "Spain",
#         refDate= 1962.5,
#         Males_five = smoothed_males,
#         Females_five = smoothed_females,
#         SRB = sex_ratio,
#         SRBDatesIn = 1960,
#         method = "linear",
#         radix = 100000
#       )
# 
#     expect_type(res, "list")
#   })
# 
# })
# 
# # IW: no need capping at 1955
# # test_that("basepop caps nLxDatesIn to 1955 when provided a date below that", {
# #
# #   tmp_nlx <- c(1954, 1960)
# #   expect_output(
# #     tmp <-
# #       basepop_five(
# #         location = "Spain",
# #         refDate = 1960,
# #         Males_five = smoothed_males,
# #         Females_five = smoothed_females,
# #         SRB = sex_ratio,
# #         nLxDatesIn = tmp_nlx,
# #         AsfrDatesIn = c(1955, 1960),
# #         method = "linear",
# #         radix = 100000
# #       ),
# #     regexp = "nLxDate\\(s\\) 1954 is/are below 1955\\. Capping at 1955",
# #     all = FALSE
# #   )
# #
# #   tmp_asfr <- c(1954, 1960)
# #   expect_output(
# #     tmp <-
# #       basepop_five(
# #         location = "Spain",
# #         refDate = 1960,
# #         Males_five = smoothed_males,
# #         Females_five = smoothed_females,
# #         SRB = sex_ratio,
# #         nLxDatesIn = c(1955, 1960),
# #         AsfrDatesIn = tmp_asfr,
# #         method = "linear",
# #         radix = 100000
# #       ),
# #     regexp = "AsfrDate\\(s\\) 1954 is/are below 1955\\. Capping at 1955",
# #     all = FALSE
# #   )
# #
# # })
# 
# test_that("basepop fails when it implies an extrapolation of > 5 years", {
# 
#   ## For nLxDatesIn ##
# 
#   ## By setting refDate to 1974, the difference between 1974 - 7.5 and the
#   ## minimum of nLxDatesIn is greater than five.
# 
#   expect_error(
#     basepop_five(
#       refDate = 1974,
#       Males_five = smoothed_males,
#       Females_five = smoothed_females,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxMale = nLxMale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       radix = 100000
#     ),
#     regexp = "nLxDatesIn implies an extrapolation of > 5 years to achieve the needed reference dates",
#     fixed = TRUE
#   )
# 
#   ## By setting refDate to 1995, the difference between 1995 - 0.5 and the
#   ## maximum of nLxDatesIn is greater than five.
# 
#   expect_error(
#     basepop_five(
#       refDate = 1995,
#       Males_five = smoothed_males,
#       Females_five = smoothed_females,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxMale = nLxMale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = AsfrDatesIn,
#       radix = 100000
#     ),
#     regexp = "nLxDatesIn implies an extrapolation of > 5 years to achieve the needed reference dates",
#     fixed = TRUE
#   )
# 
#   ## For AsfrDatesIn
#   ## Here we just provide AsfrDatesIn which we are much higher than refDate - 7.5
#   expect_error(
#     basepop_five(
#       refDate = 1986,
#       Males_five = smoothed_males,
#       Females_five = smoothed_females,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxMale = nLxMale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = c(1925, 1930),
#       radix = 100000
#     ),
#     regexp = "AsfrDatesIn implies an extrapolation of > 5 years to achieve the needed reference dates",
#     fixed = TRUE
#   )
# 
#   ## Here we just provide AsfrDatesIn which we are much higher than refDate - 0.5
#   expect_error(
#     basepop_five(
#       refDate = 1986,
#       Males_five = smoothed_males,
#       Females_five = smoothed_females,
#       SRB = sex_ratio,
#       nLxFemale = nLxFemale,
#       nLxMale = nLxMale,
#       nLxDatesIn = nLxDatesIn,
#       AsfrMat = AsfrMat,
#       AsfrDatesIn = c(2020, 2025),
#       radix = 100000
#     ),
#     regexp = "AsfrDatesIn implies an extrapolation of > 5 years to achieve the needed reference dates",
#     fixed = TRUE
#   )
# 
# })
# 
# srb_checker <- function(x, ordered_name = FALSE) {
#   expect_length(x, 3)
#   expect_named(x)
#   expect_type(x, "double")
# 
#   # Check the names of x are returned ordered
#   # This only makes sense when we download the data from WPP because
#   # otherwise we respect the order provided by the user.
#   if (ordered_name) {
#     expect_true(all(names(x) == as.character(sort(as.numeric(names(x))))))
#   }
# }
# 
# test_that("downloadSRB works as expected", {
#   # Should return same estimate three times
#   srb_checker(downloadSRB(c(1.05), DatesOut = c(1999, 1920, 1930)))
# 
#   # Should return three values with the last being 1.05
#   srb_checker(downloadSRB(c(1.05, 1.07), DatesOut = c(1999, 1920, 1930)))
# 
#   # Should return the same thing
#   srb_checker(downloadSRB(c(1.05, 1.07, 1.08), DatesOut = c(1999, 1920, 1930)))
# 
#   # Should error
#   expect_error(
#     downloadSRB(1:4),
#     regexp = "SRB can only accept three dates at maximum",
#     fixed = TRUE
#   )
# 
#   # Should return three SRBs
#   srb_checker(downloadSRB(SRB = NULL, location = "Spain", DatesOut = 1955:1957))
# 
#   # Assumes SRB
#  outp <- capture_output_lines(
#    downloadSRB(SRB = NULL,
#                location = "Whatever",
#                DatesOut = 1955:1957)
#    )
# 
#  expect_true(all(
#    c("Whatever not available in DemoToolsData::WPP2019_births",
#      "Assuming SRB to be 1.047 ") %in% outp) )
# 
#   # Should fail because of number of years
#   expect_error(
#     downloadSRB(SRB = NULL, location = "Whatever", DatesOut = 1955:1958),
#     regexp = "SRB can only accept three dates at maximum",
#     fixed = TRUE
#   )
# 
#   # Should impute the first two years with the last
#   srb_checker(
#     downloadSRB(SRB = NULL, location = "Germany", DatesOut = 1948:1950),
#     ordered_name = TRUE
#   )
# 
#   # Should impute all values
#   srb_checker(
#     downloadSRB(SRB = NULL, location = "Germany", DatesOut = 1947:1949),
#     ordered_name = TRUE
#   )
# })
# 
# # Interpolation between pivot wpp years, also into the period 1950-1955
# tfr_pj <- data.frame(
#            year = c(1950.0,1950.5,1951.5,1952.5,1953.0,
#                     1953.5,1954.5,1955.0,1955.5,
#                     1956.5,1957.5,1958.0,1958.5,
#                     1959.5,1960.0,1960.5,1961.5,
#                     1962.5,1963.0,1963.5,1964.5,
#                     1965.0,1965.5,1966.5,1967.5,
#                     1968.0,1968.5,1969.5,1970.0,
#                     1970.5,1971.5,1972.5,1973.0),
#            tfr = c(7.2986,7.3290,7.3898,7.4506,7.4810,7.5114,7.5722,7.6026,7.6330,
#                    7.6938,7.7546,7.7850,7.8130,7.8690,7.8970,7.9250,7.9810,8.0370,
#                    8.0650,8.0695,8.0785,8.0830,8.0875,8.0965,8.1055,8.1100,8.0980,
#                    8.0740,8.0620,8.0500,8.0260,8.0020,7.9900))
# 
# test_that("Replicate Peter Johnson´s excel for extrapolatebeyond 1955",{
#   expect_equal(tfr_pj$tfr,
#                downloadAsfr(Asfrmat = NULL, location = "Kenya",
#                             AsfrDatesIn = tfr_pj$year) %>% colSums() %>%
#                             as.numeric() * 5
#                )
# })
# 
# test_that("Receive a message if asked dates are not in 1950-2025 interval",{
#   expect_output(downloadAsfr(Asfrmat = NULL, location = "Kenya",
#                              AsfrDatesIn = 1900),
#                 regexp = "Careful, extrapolating beyond range 1950-2025")
#   expect_output(downloadnLx(nLx = NULL, location = "Kenya",
#                             nLxDatesIn = 1900, gender="both"),
#                 regexp = "Careful, extrapolating beyond range 1950-2025")
# })
