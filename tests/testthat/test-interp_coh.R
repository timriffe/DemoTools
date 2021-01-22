check_form <- function(x) {
  expect_is(x, "matrix")
  expect_true(nrow(x) == 101)
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
  expect_true(ncol(res) == 9) # we get an extra midyear on the right side
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

})

