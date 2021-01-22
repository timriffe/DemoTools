check_form <- function(x) {
  expect_is(x, "matrix")
  expect_true(nrow(x) == 101)
  expect_true(all(!is.na(x)))
  #expect_true(ncol(x) == 8)
  expect_true(length(colnames(x)) != 0)
}

births <- c(719511L, 760934L, 772973L, 749554L, 760831L,
            828772L, 880543L, 905380L, 919639L)

test_that("interp_coh works without midyear", {

  res <-
    interp_coh(
      country = "Russian Federation",
      sex = "male",
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      age1 = 0:100,
      births = births
    )

  check_form(res)
  expect_true(ncol(res) == 8)
})

test_that("interp_coh works with midyear", {

  res <-
    interp_coh(
      country = "Russian Federation",
      sex = "male",
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      age1 = 0:100,
      births = births,
      midyear = TRUE
    )

  check_form(res)
  # TODO: Tim said that we get an extra midyear on the right side
  # but there was a mistake where dates_out inside interp_coh was
  # returning a logical vector when it should've been a set of dates.
  # When I fixed that, it actually returned 8, so leaving it as 8
  # (before it was a 9).
  expect_true(ncol(res) == 8)
})


# Examples for interpolating between two Russian censuses


test_that("interp_coh works well with age1", {

  # 1) births given as vector
  #    mortality pulled from WPP2019  (graduated as needed)
  res1 <- interp_coh(
    country = "Russian Federation",
    sex = "male",
    c1 = pop1m_rus2002,
    c2 = pop1m_rus2010,
    date1 = "2002-10-16",
    date2 = "2010-10-25",
    age1 = 0:100,
    births = c(719511L, 760934L, 772973L, 749554L,
               760831L, 828772L, 880543L, 905380L, 919639L)
  )

  # Same, but age args totally inferred.
  res2 <- interp_coh(
    country = "Russian Federation",
    sex = "male",
    c1 = pop1m_rus2002,
    c2 = pop1m_rus2010,
    date1 = "2002-10-16",
    date2 = "2010-10-25",
    births = c(719511L, 760934L, 772973L, 749554L,
               760831L, 828772L, 880543L, 905380L, 919639L))

  expect_equal(res1, res2)
})


test_that("Births are pulled from post-processed WPP2019", {
  # 2) births pulled from post-processing of WPP2019;
  #    mortality from WPP2019 (graduated as needed)

  expect_output(
    interp_coh(
      country = "Russian Federation",
      sex = "male",
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      age1 = 0:100
    ),
    regexp = "Births fetched from WPP for: Russian Federation male population, years 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 " #nolintr
  )
})

test_that("interp_coh works well with different time points", {
  # 3) mortality (abridged, 2 and 3 time points) and fertility given: 
  mortdate1 <- 2003
  mortdate2 <- 2006
  mortdate3 <- 2010
  age_lx <- c(0,1,seq(5,100,by=5))
  lx1 <- fertestr::FetchLifeTableWpp2019(
                     locations = "Russian Federation",
                     year = mortdate1,
                     sex = "male")$lx

  lx2 <- fertestr::FetchLifeTableWpp2019(
                     locations = "Russian Federation",
                     year = mortdate2, sex = "male")$lx

  lx3 <- fertestr::FetchLifeTableWpp2019(
                     locations = "Russian Federation",
                     year = mortdate3, sex = "male")$lx

  lxmat2 <- cbind(lx1,lx3)
  lxmat3 <- cbind(lx1,lx2,lx3)

  # with 2 mort timepoints
  res1 <- interp_coh(
   c1 = pop1m_rus2002,
   c2 = pop1m_rus2010,
   date1 = "2002-10-16",
   date2 = "2010-10-25",
   lxMat = lxmat2,
   dates_lx = c(mortdate1,mortdate3),
   age_lx = age_lx,
   births = c(719511L, 760934L, 772973L, 749554L,
              760831L, 828772L, 880543L, 905380L, 919639L),
   years_births = 2002:2010)

  check_form(res1)

  # with 3 mort timepoints
  res2 <- interp_coh(
   c1 = pop1m_rus2002,
   c2 = pop1m_rus2010,
   date1 = "2002-10-16",
   date2 = "2010-10-25",
   lxMat = lxmat3,
   dates_lx = c(mortdate1,mortdate2,mortdate3),
   age_lx = age_lx,
   births = c(719511L, 760934L, 772973L, 749554L, 
              760831L, 828772L, 880543L, 905380L, 919639L),
   years_births = 2002:2010)

  check_form(res2)

  # Same as previous but with extra birth year specified (engage birth year filtering)
  res3 <- interp_coh(
   c1 = pop1m_rus2002,
   c2 = pop1m_rus2010,
   date1 = "2002-10-16",
   date2 = "2010-10-25",
   lxMat = lxmat3,
   dates_lx = c(mortdate1,mortdate2,mortdate3),
   age_lx = age_lx,
   births = c(719511L, 760934L, 772973L, 749554L, 
              760831L, 828772L, 880543L, 905380L, 919639L,1e6),
   years_births = 2002:2011)

  check_form(res3)

})

test_that("Test for stationary population using interp_coh", {
  # Test for a stationary population: Success = each year basically the same.
  LT    <- fertestr::FetchLifeTableWpp2019(
                       locations = "Russian Federation",
                       year = 2003, sex = "male"
                     )

  LT1   <- lt_abridged2single(ndx = LT$dx, nLx = LT$Lx, Age = LT$x, OAnew = 110)

  # We could get close by just taking lx from LT1,
  # But we can get even closer by converting this
  # lifetable to a PC one using the same approximations
  # used inside interp_coh.

  px    <- 1 - LT1$nqx
  pxsq  <- px ^ .5
  N     <- length(pxsq)
  # That is, the first element is just a lower tri surv prob,
  # which wrongly assumes that the survival probs in the upper
  # and lower infant triangles are equal.
  pxcp  <- c(pxsq,1) * c(1, pxsq)
  qxcp  <- 1 - pxcp

  # left and right-side stationary populalations, where radix
  # 1e5 is the horizontal birth line. Woot.
  c1    <- lt_single_qx(qxcp, OAnew = 110)$lx[-1][1:101]
  c2    <- c1

  lxMat <- cbind(LT1$lx,LT1$lx,LT1$lx)

  Pxt <- interp_coh(
    c1 = c1,
    c2 = c2,
    date1 = "2002-01-01",
    date2 = "2010-01-01",
    lxMat = lxMat,
    # linear interp, would be same w 2 or 10 lx columns.
    dates_lx = c(2002,2005,2008),
    age_lx = 0:110,
    births = rep(1e5,9), # stationary birth series
    years_births = 2002:2010)

  # here's the test:
  # now that's what I call a deterministic stationary population.
  # :-)
  expect_true(all(abs(diff(t(Pxt))) < 1e9))
})


test_that("interp_coh errors if not given correctly", {
  # 1) births given (no years_birth), but not right length
  expect_error(
    interp_coh(
      country = "Russian Federation",
      sex = "male",
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      age1 = 0:100,
      # Here we provide births with one year less
      births = births[-length(births)]
    )
  )

  # 2) births given, correct length, but not right years
  expect_error(
    interp_coh(
      country = "Russian Federation",
      sex = "male",
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      age1 = 0:100,
      # Correct births
      births = births,
      # Incorrect years of birth (should be 2010)
      years_births = 2002:2009
    )
  )
})

# We should error if
test_that("interp_coh fails when lxmat is not correct", {
  # Downloads some data
  mortdate1 <- 2003
  mortdate2 <- 2006
  mortdate3 <- 2010
  age_lx <- c(0,1,seq(5,100,by=5))
  lx1 <- fertestr::FetchLifeTableWpp2019(
                     locations = "Russian Federation",
                     year = mortdate1,
                     sex = "male")$lx

  lx2 <- fertestr::FetchLifeTableWpp2019(
                     locations = "Russian Federation",
                     year = mortdate2, sex = "male")$lx

  lx3 <- fertestr::FetchLifeTableWpp2019(
                     locations = "Russian Federation",
                     year = mortdate3, sex = "male")$lx

  lxmat <- cbind(lx1,lx2,lx3)

  # 3.1) lxMat given, but only one column
  expect_error(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      lxMat = lxmat[, 1, drop = FALSE],
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L, 
                 760831L, 828772L, 880543L, 905380L, 919639L),
      years_births = 2002:2010),
    regexp = "lxMat should have at least two or more dates as columns. lxMat contains only one column" #nolintr
  )

  ## 3.2) lxMat give, but the date range in it doesn't overlap
  ## with the date range of date1 to date2 (i.e. 100% extrapolation implied)
  expect_output(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2000-10-16",
      date2 = "2014-10-25",
      # Make up some very long dates
      lxMat = lxmat[, 1:2],
      dates_lx = c(2007, 2008),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L,
                 760831L, 828772L, 880543L, 905380L,
                 919639L, 719511L, 760934L, 772973L,
                 749554L, 760831L, 828772L),
      years_births = 2000:2014),
    regexp = "Range between `date1` and `date2` must overlap with `lx_dates` for at least 25% of the range or 6 years." #nolintr
  )

  # Full error when dates_lx are now within the date1 and date2 threshold.
  expect_error(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2000-10-16",
      date2 = "2014-10-25",
      # Make up some very long dates
      lxMat = lxmat[, 1:2],
      dates_lx = c(2020, 2021),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L,
                 760831L, 828772L, 880543L, 905380L,
                 919639L, 719511L, 760934L, 772973L,
                 749554L, 760831L, 828772L),
      years_births = 2000:2014),
    regexp = "All `dates_lx` must be within the range of `date1` and `date2`"
  )

})


# 4) age1 or age2 not single
test_that("Ages must be single in interp_coh", {

  # The error tests that they are the same length.
  # If ages are of anything other than single ages,
  # this will fail, capturing that the ages should
  # be single ages.
  expect_error(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      # Supply ages in five year age groups
      age1 = seq(0, 100, by = 5),
      lxMat = lxmat,
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L,
                 760831L, 828772L, 880543L, 905380L, 919639L),
      years_births = 2002:2010),
    regexp = "length(age1) == length(c1) is not TRUE",
    fixed = TRUE
  )

  expect_error(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      # Supply ages for second age group in five year age groups
      age2 = seq(0, 100, by = 5),
      lxMat = lxmat,
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L,
                 760831L, 828772L, 880543L, 905380L, 919639L),
      years_births = 2002:2010),
    regexp = "length(age2) == length(c2) is not TRUE",
    fixed = TRUE
  )

  expect_error(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      # Both ages supplied
      age1 = seq(0, 100, by = 5),
      age2 = seq(0, 100, by = 5),
      lxMat = lxmat,
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L,
                 760831L, 828772L, 880543L, 905380L, 919639L),
      years_births = 2002:2010),
    regexp = "length(age1) == length(c1) is not TRUE",
    fixed = TRUE
  )
})

test_that("interp_coh fails if arguments not supplied to download data ", {

  # 5) no births given, and no country/sex given
  expect_error(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      lxMat = lxmat,
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      age_lx = age_lx,
      years_births = 2002:2010),
    regexp = "births not specified, please specify country and sex",
    fixed = TRUE
  )

  # 6) no lxMat given, and no country/sex given
  expect_error(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L,
                 760831L, 828772L, 880543L, 905380L, 919639L),
      years_births = 2002:2010),
    regexp = "lxMat not specified, please specify country and sex",
    fixed = TRUE
  )
})

test_that("c1, c2 and lxmat should not have negatives", {

  # 7) c1, c2, lxMat, or births have negatives

  c1_neg <- pop1m_rus2002
  c1_neg[1] <- -c1_neg[1]

  expect_error(
    interp_coh(
      c1 = c1_neg,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      lxMat = lxmat,
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L, 760831L,
                 828772L, 880543L, 905380L, 919639L),
      years_births = 2002:2010
    ),
    regexp = "No negative values allowed in `c1`"
  )

  c2_neg <- pop1m_rus2010
  c2_neg[1] <- -c2_neg[1]

  expect_error(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = c2_neg,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      lxMat = lxmat,
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L, 760831L,
                 828772L, 880543L, 905380L, 919639L),
      years_births = 2002:2010
    ),
    regexp = "No negative values allowed in `c2`"
  )

  lxmat_neg <- lxmat
  lxmat_neg[2, 1] <- -lxmat_neg[2, 1]

  expect_error(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      lxMat = lxmat_neg,
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L, 760831L,
                 828772L, 880543L, 905380L, 919639L),
      years_births = 2002:2010
    ),
    regexp = "No negative values allowed in `lxMat`"
  )

})


test_that("interp_coh shows appropriate warnings when verbose = TRUE", {

  # 1) age1 and age2 not same range
  # TODO
  ## expect_error(
  ##   interp_coh(
  ##     c1 = pop1m_rus2002,
  ##     c2 = pop1m_rus2010,
  ##     date1 = "2002-10-16",
  ##     date2 = "2010-10-25",
  ##     lxMat = lxmat,
  ##     dates_lx = c(mortdate1,mortdate2,mortdate3),
  ##     age_lx = age_lx,
  ##     births = c(719511L, 760934L, 772973L, 749554L, 760831L,
  ##                828772L, 880543L, 905380L, 919639L),
  ##     years_births = 2002:2010,
  ##     verbose = TRUE
  ##   ),
  ##   regexp = "No negative values allowed in `lxMat`"
  ## )

  # 2) date2 - date1 > 15
  expect_output(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      # Here I set the year to 2020
      date2 = "2017-10-25",
      lxMat = lxmat,
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      age_lx = age_lx,
      # Add fake births/years_births so that they exceed more
      # than 15 years
      births = c(719511L, 760934L, 772973L, 749554L, 760831L,
                 828772L, 880543L, 905380L, 919639L, 919639L,
                 760831L, 880543L, 719511L, 760934L, 772973L,
                 749554L),
      years_births = 2002:2017,
      verbose = TRUE
    ),
    regexp = "FYI, there are 15.02466 years between c1 and c2\nBe wary.",
    fixed = TRUE
  )

  # 3) if the shortest distance from dates_lx to date1 or date2 is greater than 7
  # TODO
  expect_output(
    interp_coh(
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2000-10-16",
      date2 = "2017-10-25",
      lxMat = lxmat[, 1:2],
      dates_lx = c(2008, 2009),
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L, 760831L,
                 828772L, 880543L, 905380L, 919639L, 719511L,
                 760934L, 772973L, 749554L, 760831L, 828772L,
                 749554L, 760831L, 828772L),
      years_births = 2000:2017,
      verbose = TRUE
    ),
    regexp = "The shortest distance from `dates_lx` ( 2008 ) to `date1/date2`( 2000.79 ) is greater than 7 years. Be wary.",
    fixed = TRUE
  )

  # 4) any negatives detected in output (to be imputed with 0s)
  # TODO

})

test_that("interp_coh throws download messages when verbose = TRUE", {

  # 1) lx is downloaded
  expect_output(
    interp_coh(
      country = "Russian Federation",
      sex = "both",
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L, 760831L,
                 828772L, 880543L, 905380L, 919639L),
      years_births = 2002:2010,
      verbose = TRUE
    ),
    regexp = "lxMat not provided. Downloading lxMat for Russian Federation, gender: `both`, for years between 2002.8 and 2010.8"
  )

  # 2) births are downloaded
  expect_output(
    interp_coh(
      country = "Russian Federation",
      sex = "both",
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      lxMat = lxmat,
      dates_lx = c(mortdate1,mortdate2,mortdate3),
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      age_lx = age_lx,
      verbose = TRUE
    ),
    regexp = "Births fetched from WPP for: Russian Federation both population, years 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010"
  )

  # 3) dates_lx or years_births are being assumed anything
  expect_output(
    interp_coh(
      country = "Russian Federation",
      sex = "both",
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      lxMat = lxmat,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      age_lx = age_lx,
      births = c(719511L, 760934L, 772973L, 749554L, 760831L,
                 828772L, 880543L, 905380L, 919639L),
      years_births = 2002:2010,
      verbose = TRUE
    ),
    regexp = "lxMat specified, but not dates_lx\nAssuming: 2002.78904109589, 2006.80136986301, 2010.81369863014",
    fixed = TRUE
  )
})
