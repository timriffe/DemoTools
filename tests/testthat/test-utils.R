# Author: ilya
###############################################################################
context("test-utils")

test_that("shift.vector works", {
    vec <- 1:5
    expect_length(shift.vector(vec), 5)
    expect_length(shift.vector(vec, shift = 2), 7)
    expect_type(shift.vector(vec, shift = 2, fill = "new"), "character")
    expect_type(shift.vector(vec, shift = 2, fill = .1), "double")
})


test_that("ma works",{
    set.seed(911)
    pop <- sample.int(10, 5, replace = T)
    expect_equal(head(ma(pop, 2), 2), c(3, 6))
    expect_length(na.omit(ma(pop, 2)), 4)
    expect_length(na.omit(ma(pop, 3)), 3)
})


test_that("rescale.vector works",{
    set.seed(911)
    x <- runif(10)
    xx <- rescale.vector(x, 100)
    expect_equal(sum(xx), 100)
    # check if proportionality is preserved
    expect_equal(x[1]/x[2], xx[1]/xx[2])
})



test_that("is_LeapYear works",{
    expect_equal(sum(is_LeapYear(1895:1905)), 2)
    expect_equal(sum(is_LeapYear(1995:2005)), 3)
})



test_that("ypart works",{
    expect_equal(ypart(2001, 2, 14), 0.1150685, tolerance = 1e-8)
    expect_equal(ypart(2001, 6, 30), .5)
    expect_equal(ypart(2001, 7, 1), .5)
    expect_equal(ypart(2000, 7, 1), .5)
    expect_equal(ypart(2000, 7, 1, FALSE), .5)
    expect_equal(ypart(2000, 7, 1, FALSE), .5)
    expect_false(ypart(2001, 7, 1, FALSE) == .5)
    expect_equal(ypart(2002, 12, 31, detect.start.end = FALSE), 1)
    expect_equal(ypart(2002, 1, 1, detect.start.end = FALSE), 
                 0.002739726, tolerance = 1e-8)
    expect_equal(ypart(2002, 1, 1, detect.start.end = TRUE), 0)
})


test_that("dec.date works",{
    expect_equal(dec.date("2018-10-01"), 2018.751, tolerance = 1e-3)
    expect_equal(dec.date(structure(17805, class = "Date")), 2018.751, 
                 tolerance = 1e-3)
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
    one <- splitUniform(MalePop, Age = Ages)
    expect_equivalent(tail(as.vector(table(one)), 2), c(5, 1))
})

