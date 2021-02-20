# TODO
# -[ ] check for success in a few different input scenarios
# -[ ] check sums match between input and output
# -[ ] check output of proper length
# -[ ] check for console output when informative warnings should be generated
# -[ ] expect error if OAnew > max(Age_nLx)
# -[ ] expect errors in other reasonable situations
# -[ ] canonical test: if the census is actually a stationary population, identical to nLx itself, then we hope to return something proportional to it.
# -[ ] canonical test: if input census is an exact stable pop (use auxiliary function OPAG_nLx_warp_r)

###############################################################################



context("test-OPAG")

# data --------------------------------------------------------------------
Pop   <- 
Age   <- c(0:85)
OAnow <- c()
StPop <- c()
StAge <- c(86:100)
OAnew <- max(StAge)
nLx <- 
r <- 0.005
AgeInt <- c()

# Insert data to be compared

test_that("OPAG_simple works", {
  expect_equal(
    OPAG_simple(
      Pop = Pop,
      Age = Age,
      OAnow = max(Age),
      StPop = StPop,
      OAnew = max(StAge)
    ),
    0.444659710063221, # think about the type of data would return and make equal
    tolerance = 1e-12
  )}
  )

test_that("OPAG_nLx_warp_r works", {
  expect_equal(
    OPAG_nLx_warp_r(
      nLx = nLx,
      Age = Age,
      r = r,
      AgeInt = NULL,
      continuous = TRUE,
      method = "uniform"
    ),
    0.444659710063221, # think about the type of data would return and make equal
    tolerance = 1e-12
  )}
)


test_that("OPAG_fit_stable_standard works", {
  expect_equal(
    OPAG_fit_stable_standard(Pop_fit,
                             Age_fit,
                             AgeInt_fit,
                             nLx,
                             Age_nLx,
                             AgeInt_nLx,          
                             method = "uniform",
                             continuous = TRUE), 
               value,
               tolerance = 0,0001)
})
  
