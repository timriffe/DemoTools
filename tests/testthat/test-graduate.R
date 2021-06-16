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

test_that("age ranges flexible for uniform, mono, and pclm", {
  Pop5    <- c(88962, 73756, 51990, 55395, 48562)
  Age5    <- seq(0,20,5)
  Age5alt <- seq(45,65,5)
  
  m0      <- graduate_mono( Pop5, Age5 )
  m45     <- graduate_mono( Pop5, Age5alt ) 
  p0      <- graduate_pclm( Pop5, Age5 )
  p45     <- graduate_pclm( Pop5, Age5alt ) 
  u0      <- graduate_uniform( Pop5, Age5 )
  u45     <- graduate_uniform( Pop5, Age5alt ) 
  
  expect_true(all(abs(m0 - m45) < 1e-9))
  expect_true(all(abs(p0 - p45) < 1e-4))
  expect_true(all(abs(u0 - u45) < 1e-9))
  
  ma0     <- names2age(m0)
  ma45    <- names2age(m45)
  pa0     <- names2age(p0)
  pa45    <- names2age(p45)
  ua0     <- names2age(u0)
  ua45    <- names2age(u45)
  
  expect_equal(ma0, 0:20, tolerance = 0)
  expect_equal(pa0, 0:20, tolerance = 0)
  expect_equal(ua0, 0:20, tolerance = 0)
  expect_equal(ma45, 45:65, tolerance = 0)
  expect_equal(pa45, 45:65, tolerance = 0)
  expect_equal(ua45, 45:65, tolerance = 0)
  
})


