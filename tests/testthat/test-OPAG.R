# TODO
# -[x ] check for success in a few different input scenarios
# -[ x] check sums match between input and output
# -[ x] check output of proper length
# -[ ] check for console output when informative warnings should be generated
# -[ ] expect error if OAnew > max(Age_nLx)
# -[ ] expect errors in other reasonable situations
# -[ ] canonical test: if the census is actually a stationary population, identical to nLx itself, then we hope to return something proportional to it.
# -[ ] canonical test: if input census is an exact stable pop (use auxiliary function OPAG_nLx_warp_r)

###############################################################################



context("test-OPAG")
# pop sweden 2000
pop_swe <-  c(88367.00, 89890.00, 91008.00, 95773.00, 103678.00,
112627.00, 117695.00, 123777.00, 126366.00, 127545.00,
120657.00, 117454.00, 110691.00, 108696.00, 105477.00,
101220.00, 99075.00, 100090.00, 100684.00, 103623.00,
102193.00, 99242.00, 102375.00, 105360.00, 111201.00,
117816.00, 117423.00, 120121.00, 121586.00, 118213.00,
117669.00, 124087.00, 131642.00, 133692.00, 134292.00,
134758.00, 125367.00, 120009.00, 116535.00, 114949.00,
116379.00, 116192.00, 117963.00, 118680.00, 117149.00,
114963.00, 118060.00, 117424.00, 116446.00, 121530.00,
126177.00, 130623.00, 131833.00, 133485.00, 132313.00,
130098.00, 121803.00, 111254.00, 98568.00, 93645.00,
94070.00, 89995.00,  85443.00, 82849.00, 79006.00,
76534.00, 74682.00,  76992.00, 76164.00,  76840.00,
73736.00, 75105.00,  72258.00, 72409.00,  72466.00,
71159.00, 70307.00,  67816.00, 69750.00,  69855.00,
54218.00, 50316.00,  47688.00, 43035.00,  38938.00,
36031.21, 31859.38,  27876.99, 23447.49,  19537.57,
16325.19, 12829.43,   9890.83, 7421.46,   5308.28,
3849.87, 2690.06, 1762.22, 1119.00, 673.00,
386.00, 227.00, 127.00, 79.00, 43.00,
12.00, 10.00, 3.00, 4.00, 0.00,
2.00)


pop_swe50_110 <- pop_swe[51:111]
pop_85 <- c(pop_swe[1:85], sum(pop_swe[86:111]))
pop_check <- pop_swe
names_pop_check <- as.character(c(0:110))
names(pop_check)<- names_pop_check

# data --------------------------------------------------------------------
Pop   <- pop_85
Age   <- c(0:85)
OAnow <- max(Age)
StPop <- pop_swe50_110
StAge <- c(50:110)
OAnew <- max(StAge)
nLx <- 
r <- 0.005
AgeInt <- c()

# Insert data to be compared

test_that("OPAG_simple works", {
  OPAG_res <- OPAG_simple(
    Pop = Pop,
    Age = Age,
    OAnow = max(Age),
    StPop = StPop,
    StAge = StAge,
    OAnew = max(StAge)
  )
  expect_equal(OPAG_res,
    pop_check, # think about the type of data would return and make equal
    tolerance = 0.00001
  )}
  )


# test sum

test_that("OPAG_simple's sum of input and output are equal", {
  OPAG_res <- OPAG_simple(
    Pop = Pop,
    Age = Age,
    OAnow = max(Age),
    StPop = StPop,
    StAge = StAge,
    OAnew = max(StAge)
  )
  expect_equal(sum(OPAG_res),
               sum(pop_check), # think about the type of data would return and make equal
               tolerance = 0.00001
  )}
)

#length

test_that("OPAG_simple's output has a proper length", {
  OPAG_res <- OPAG_simple(
    Pop = Pop,
    Age = Age,
    OAnow = max(Age),
    StPop = StPop,
    StAge = StAge,
    OAnew = max(StAge)
  )
  expect_equal(length(OPAG_res),
               length(pop_check), # think about the type of data would return and make equal
               tolerance = 0.00001
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
  
