# test for basic interpolation function `interp`

context("test-AGEINT")

test_that("basic interpolation function",{
  
  # Toy example: two ages, two dates, expected result by eye
  popmatr <- matrix(1:2,2,2,byrow = T)
  dates <- dec.date(c(2000.5,2010.5))
  
  # interp
  expect_true(
    all( 1.5 ==
           interp(popmat = popmatr, 
                  datesIn = dates,
                  datesOut = dec.date(2005.5))
    )
  )
  # extrap
  expect_true(
    all( 2.5 ==
           interp(popmat = popmatr, 
                  datesIn = dates,
                  datesOut = dec.date(2015.5), extrap = T)
    )
  )
  # constant inputs
  expect_true(
    all( 1 ==
           interp(popmat = matrix(1,2,2,byrow = T), 
                  datesIn = dates,
                  datesOut = dec.date(2005.5)),
         1 ==
           interp(popmat = matrix(1,2,2,byrow = T), 
                  datesIn = dates,
                  datesOut = dec.date(2015.5), extrap = T)
    )
  )
  
  # you set an output date beyond observed range, but without extrap=T, receive NA
  # this is an interp function by default, so the argument to change the deafult behaviour
  # must be explicit
  expect_true(
    all(
      is.na(
        interp(popmat = popmatr, 
               datesIn = dates,
               datesOut = dec.date(2015.5))
      )
    )
  )
  
  # no negative values
  expect_true(
    all(interp(popmat = popmatr, 
               datesIn = dates,
               datesOut = dec.date(1900.5), 
               extrap = T)>=0
    )
  )
  expect_output(interp(popmat = popmatr, 
                       datesIn = dates,
                       datesOut = dec.date(1900.5), 
                       extrap = T),
                regexp = "Negative values were turned 0. No accepted in population counts, fertility rates or life table functions.")
})
