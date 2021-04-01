# Author: ilya

context("test-agesmth")


# smooth_age_5_cf -----------------------------------------------------

test_that("smooth_age_5_cf works", {
  Ages <- seq(0, 80, by = 5)

  CFmales <- smooth_age_5_cf(pop5m_pasex, Ages, OAG = TRUE)

  CFtest <- c(
    NA, NA, 346290, 287083, 285855, 261082, 237937,
    202809, 162973, 125720, 88730, 67352, 55187, 40657, NA, NA, NA
  )

  # test
  expect_equivalent(CFmales, CFtest, tolerance = 1)
  expect_equal(CFmales %>% sum(), CFtest %>% sum(), tolerance = 1)
  expect_true(all(CFmales > 0, na.rm = T))
})



# smooth_age_5_kkn ----------------------------------------------------------------

test_that("smooth_age_5_kkn works", {
  Ages <- seq(0, 80, by = 5)

  KKNtest <- c(
    NA, NA, 354871, 278502, 285508, 261429, 236513,
    204233, 162138, 126555, 90094, 65988, 54803, 41041, NA, NA, NA
  )

  KKNmales <- smooth_age_5_kkn(pop5m_pasex, Ages, OAG = TRUE)

  # test
  expect_equivalent(KKNmales, KKNtest, tolerance = 1)
  expect_equal(KKNmales %>% sum(), KKNtest %>% sum(), tolerance = 1)
  expect_true(all(KKNmales > 0, na.rm = T))
})


# smooth_age_5_arriaga ----------------------------------------------------------------

test_that("smooth_age_5_arriaga works", {
  Ages <- seq(0, 80, by = 5)

  Amales <- smooth_age_5_arriaga(Value = pop5m_pasex, Age = Ages, OAG = TRUE)

  #' # PAS spreadsheet result:
  Atest <- c(
    662761, 495126, 345744, 287629, 285919, 261018, 237469, 203277,
    161733, 126960, 88586, 67496, 54587, 41257, 28790, 17189, 34729
  )

  # test
  expect_equivalent(Amales, Atest, tolerance = 1)
  expect_equal(Amales %>% sum(), Atest %>% sum(), tolerance = 1)
  expect_true(all(Amales > 0, na.rm = T))
})


# smooth_age_5_strong ----------------------------------------------------------------

test_that("smooth_age_5_strong works", {
  Ages <- seq(0, 80, by = 5)

  Stest <- c(
    646617, 511270, 386889, 317345, 273736, 240058, 218645, 188297,
    153931, 124347, 93254, 71858, 53594, 39721, 27887, 18092, 34729
  )

  Smales <- smooth_age_5_strong(pop5m_pasex, Ages, OAG = TRUE)

  # test
  expect_equivalent(Smales, Stest, tolerance = 1)
  expect_equal(Smales %>% sum(), Stest %>% sum(), tolerance = 1)
  expect_true(all(Smales > 0, na.rm = T))
})

test_that("smooth_age_5_mav works", {

  # Without tails
  Ages <- seq(0, 80, by = 5)

  Mtest <- c(
    NA, 505239, 382964, 300570, 274160, 263151, 239782, 202228,
    162308, 128489, 92946, 73183, 51717, 41880, 26119, NA, NA
  )

  Mmales <- smooth_age_5_mav(pop5m_pasex, Ages, OAG = TRUE, tail = FALSE)

  # test
  expect_equivalent(Mmales, Mtest, tolerance = 1)
  expect_equal(Mmales %>% sum(), Mtest %>% sum(), tolerance = 1)
  expect_true(all(Mmales > 0, na.rm = T))

  # With tails
  Mtest_tail <- c(
    642367, 507810, 382964, 300570, 274160, 263151, 239782, 202228,
    162308, 128489, 92946, 73183, 51717, 41880, 27038, 29796, 34729
  )

  Mmales_tail <- smooth_age_5_mav(pop5m_pasex, Ages, OAG = TRUE, tail = TRUE)

  # test
  expect_equivalent(Mmales_tail, Mtest_tail, tolerance = 1)
  expect_equal(Mmales_tail %>% sum(), Mtest_tail %>% sum(), tolerance = 1)
  expect_true(all(Mmales_tail > 0, na.rm = T))

  # Check tail/no_tail are the same except
  # for the tails:
  expect_equivalent(Mmales[3:15], Mmales_tail[3:15], tolerance = 1)
  # Exclue first/last pair of values because cascading alters both
  # numbers.

})



test_that("smooth_age_5_feeney works", {
  Pop <- c(
    2337, 3873, 3882, 3952, 4056, 3685, 3687, 3683, 3611, 3175,
    3457, 2379, 3023, 2375, 2316, 2586, 2014, 2123, 2584, 1475,
    3006, 1299, 1236, 1052, 992, 3550, 1334, 1314, 1337, 942,
    3951, 1128, 1108, 727, 610, 3919, 1221, 868, 979, 637,
    3409, 887, 687, 533, 313, 2488, 677, 426, 524, 333,
    2259, 551, 363, 290, 226, 1153, 379, 217, 223, 152,
    1500, 319, 175, 143, 89, 670, 149, 96, 97, 69,
    696, 170, 60, 38, 23, 745
  )

  Ages <- c(0:75)

  result <- smooth_age_5_feeney(Pop, Ages, maxit = 200, OAG = TRUE)
  # inlcude original unstated age (15)
  result <- rescale_vector(result, sum(result) + 15)
  tab2_answer <- c(
    18004, 17351, 14018, 10927, 8837, 8145, 7823, 7029,
    5748, 4326, 3289, 2415, 1794, 1197, 982, 741
  )
  names(tab2_answer) <- seq(0, 75, 5)


  testthat:::expect_equal(
    result,
    tab2_answer,
    tolerance = .001
  ) # TR: this on relative scale??
})
#' Age <- c(0,1,seq(5,90,by=5))
#' # defaults
#' zz <- zigzag_smth(dth5_zigzag, Age, OAG = TRUE, ageMin = 40, ageMax = 90)
