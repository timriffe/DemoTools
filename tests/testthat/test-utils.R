# Author: ilya
###############################################################################
context("test-utils")
suppressWarnings(RNGversion("3.5.0"))

test_that("shift.vector works", {
    vec <- 1:5
    expect_length(shift.vector(vec), 5)
    expect_type(shift.vector(vec, shift = 2, fill = "new"), "character")
    expect_type(shift.vector(vec, shift = 2, fill = .1), "double")
})

x <- 1:10
#' stopifnot(all(ma(x,2) == x, na.rm = TRUE))
#' stopifnot(all(ma(x,3) == x, na.rm = TRUE))
#' stopifnot(all(ma(x,4) == x, na.rm = TRUE))
#' stopifnot(all(abs(ma(x,5) - x) < 1e-10, na.rm = TRUE))
test_that("ma works",{
    x  <- 1:10
    
    d2 <- ma(x,2) - x
    d3 <- ma(x,3) - x
    d4 <- ma(x,4) - x
    d5 <- ma(x,5) - x
    # pop <- sample.int(10, 5, replace = T)
    # all should give same for linear data.
    expect_true(all(abs(d2) < 1e-10, na.rm = TRUE))
    expect_true(all(abs(d3) < 1e-10, na.rm = TRUE))
    expect_true(all(abs(d4) < 1e-10, na.rm = TRUE))
    expect_true(all(abs(d5) < 1e-10, na.rm = TRUE))
    
    # can anticipate NAs in tails
    expect_true(sum(is.na(d2)) == 2)
    expect_true(sum(is.na(d3)) == 2)
    expect_true(sum(is.na(d4)) == 4)
    expect_true(sum(is.na(d5)) == 4)
    # expect_length(na.omit(ma(pop, 2)), 4)
    # expect_length(na.omit(ma(pop, 3)), 3)
})


test_that("rescale_vector works",{
    set.seed(911)
    x  <- runif(10)

    xx <- rescale_vector(x, 100)
    expect_equal(sum(xx), 100)
    # check if proportionality is preserved
    expect_equal(x[1]/x[2], xx[1]/xx[2])
})


test_that("dec.date works",{
  expect_equal(dec.date("2018-10-01"), 2018.75, tolerance = 1e-3)
  date_obj <- structure(17805, class = "Date")
  expect_equal(dec.date(date_obj), 2018.75, tolerance = 1e-3)
})

test_that("retx works",{
    expect_equal(ratx(seq(-10, 10, 3)), 
                 c(.7, .57, .25, -2, 2.5, 1.6), 
                 tolerance = 1e-2)
})



test_that("splitUniform works",{
    MalePop <- c(9544406, 7471790, 11590109, 11881844, 11872503, 12968350,
                 11993151, 10033918, 14312222, 8111523, 15311047, 6861510,
                 13305117, 7454575, 9015381, 10325432, 9055588, 5519173)
    Ages <- seq(0, 85, by = 5)
    one <- graduate_uniform(MalePop, Age = Ages)
    expect_equivalent(tail(as.vector(table(one)), 2), c(5, 1))
})

