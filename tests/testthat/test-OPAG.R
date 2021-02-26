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

## testing warnings
# -[ ] check for console output when informative warnings should be generated

test_that("OAnew > max(Age_nLx) error", {
  expect_error(OPAG_simple(
    Pop = Pop,
    Age = Age,
    OAnow = max(Age),
    StPop = StPop,
    StAge = StAge,
    OAnew = max(StAge) + 5
  ), "OAnew %in% StAge is not TRUE")
})

test_that("length(Pop) == length(Age) error", {
  expect_error(OPAG_simple(
    Pop = pop_swe,
    Age = Age,
    OAnow = max(Age),
    StPop = StPop,
    StAge = StAge,
    OAnew = max(StAge)
  ), "length(Pop) == length(Age) is not TRUE")
})

## Stationary population ---------------------------------------
## Data

Lx <- c(0.997206968803419, 3.98651416246245, 4.98081195476269, 
                  4.97832219202643, 4.97409262747162, 4.96800317067612, 4.96130094726659, 
                  4.95324568198643, 4.9414967090444, 4.92251152310738, 4.89050587164337, 
                  4.83787412607811, 4.75471527967275, 4.62736885461426, 4.43325496677287, 
                  4.13111130068528, 3.65288437783847, 2.9050982849507, 1.88475330956745, 
                  0.865610023545313, 0.237771057845129, 0.0331623958830273)

r <-0
age_Lx     <- c(0,1,seq(5,100,by=5))
names(Lx)  <- age_Lx

test_that("OPAG_nLx_warp_r works", {
  res_stationary <- OPAG_nLx_warp_r(
    nLx = c(Lx),
    Age = age_Lx,
    r = r,
    continuous = TRUE,
    method = "uniform"
  )
  c_St <- res_stationary/sum(res_stationary)
  c_Lx <- c(Lx/sum(Lx))

  res <- (c_St - c_Lx) %>% abs() %>% max()
  expect_true(res < 0.0001)
  }
)


# Pop_fit was generated this way
# Pop_fit <- OPAG_nLx_warp_r(
#   nLx = c(Lx),
#   Age = age_Lx,
#   r = 0.01,
#   continuous = TRUE,
#   method = "uniform"
# )




# test_that("OPAG_fit_stable_standard works", {
#     # OPAG_nLx_warp_r(
#     #   nLx = c(Lx),
#     #   Age = age_Lx,
#     #   r = 0.01,
#     #   continuous = TRUE,
#     #   method = "uniform"
#     # )
#     
#   PopCheckStable <- c(0.0178828625366182,  0.0697292229670445, 
#                       0.0832903513014252,  0.0791886290502868, 0.075262556870,
#                       0.0715043132981018,  0.0679252465454812, 0.064507587072,
#                       0.0612159669721477,  0.0580067082482588, 0.054818928190,
#                       0.0515841876671949,  0.0482249540993725, 0.044644370355,
#                       0.0406855853175052,  0.0360636791712615, 0.03033362662,
#                       0.0229474544613544,  0.0141616376235189, 0.0061868062463,
#                       0.00161654762698326, 0.000218777749686752)
#   names(PopCheckStable) <- age_Lx
#   Pop_in                <- PopCheckStable[1:18] 
#   Pop_in[18]            <- sum(PopCheckStable[18:22])
#   Pop_in                <- Pop_in * 5e5
#   Age_Pop_in            <- names2age(Pop_in)
#   
#   AgeInt_in <- inferAgeIntAbr(Age_Pop_in, OAG = TRUE, OAvalue = 1)
#   AgeInt_nLx <- inferAgeIntAbr(age_Lx, OAG = TRUE, OAvalue = 1)
#   
#   
#   Pop_fit <- OPAG(Pop_in,
#                   Age_Pop = Age_Pop_in,
#                   AgeInt_Pop = AgeInt_in,
#                   nLx = nLx,
#                   Age_nLx = age_Lx,
#                   AgeInt_nLx = AgeInt_nLx,
#                   Age_fit =  c(50,60,70),
#                   AgeInt_fit = c(10,10,10),
#                   Redistribute_from = 80,
#                   continuous = TRUE,
#                   method = "uniform")
#   PopSt_Out <- rescale_vector(Pop_fit$Pop_out)
#   expect_equal(
#     PopSt_Out, 
#     PopCheckStable,
#                tolerance = 0.0001)
# })




