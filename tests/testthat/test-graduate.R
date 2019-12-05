context("test-graduate")

# todo: graduate_pclm, graduate_sprague, graduate_mono, graduate_uniform


# graduate_grabill
test_that("graduate_grabill_expand works",{
  a5 <- as.integer(rownames(pop5_mat))
  expect_equal(graduate_grabill_expand(pop5_mat[,1], Age = a5, OAG = TRUE)[1:3, 1], 
               c(0.0826, 0.0903, 0.0932), tolerance = 1e-3)
  expect_equal(graduate_grabill_expand(pop5_mat[,1], Age = a5, OAG = TRUE)[101, 21], 1)
})


test_that("graduate_grabill works",{
  check <- c(9350425.61, 10050059.85)
  names(check) <- c(0,1)
  expect_equal(c(graduate_grabill(pop1m_ind, 0:100, TRUE)[1:2]), 
               check, 
               tolerance = 1e-2)
})

# graduate_beers
# data --------------------------------------------------------------------

# Johnson_2016_BEERSP.XLS

M <- c(184499, (752124 - 184499), 582662, 463534, 369976, 286946, 235867,
       199561, 172133, 151194, 131502, 113439, 95614,
       78777, 60157, 40960, 21318, 25451)
Age <- c(0, 1, seq(5, 80, by = 5))
Age0        <- 184499
johnson     <- graduate_beers(
  Value = M,
  Age = Age,
  OAG = TRUE,
  method = "ord",
  johnson = TRUE)
names(johnson) <- NULL
# copied from spreadsheet output,
output <- c(184499,158163,143895,135416,130151,
            126153,122028,117182,111566,105733,
            101088,96719,92551)


# tests -------------------------------------------------------------------

test_that("beers works", {
  # allow for rounding differences, so maximum absolute difference of 1
  expect_equal(
    johnson[1:length(output)], 
    output, 
    tolerance = 1
  )
})



