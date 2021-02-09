# tests against spreadsheet "Li_2018_Limited_Lee-Carter-v4.xlsm" -------------------------------------------------------------

context("test-interp_lc_lim")

# summary
# A to E tests againts spreadsheet
# F testing args from main fun
# G mixing input: single/abr with output single/abr, and mixing input nMx and lx
# H passing lt arguments
# I text/messages/warnings

# tolerance: .05% difference in e_0
tolerance_admited <- .005 

# data included in spredsheet -------------------------------------------------------------

Age <- c(0,1,seq(5,100,5))
dates_in <- c(1980.671, 1991.668, 1996.586, 2000.586, 2010.71)
input <- data.frame(Date= c(rep(sort(rep(dates_in, length(Age))),2)),
                     Age = rep(Age,  2 * length(dates_in)),
                     Sex = c(rep("m", length(Age) * length(dates_in)),
                             rep("f", length(Age) * length(dates_in))),
                     nMx = c(0.058872, 0.002992, 0.000697, 0.000703, 0.001424, # males
                             0.002102, 0.002519, 0.003109, 0.004072, 0.005968, 0.00739, 0.010927, 
                             0.013366, 0.018798, 0.028653, 0.037304, 0.052714, 0.059629, 0.070922, 
                             0.093256, 0.135563, 0.217859, 0.030825, 0.001351, 0.000517, 0.000528, 
                             0.001434, 0.002436, 0.002954, 0.003341, 0.003971, 0.004966, 0.006267, 
                             0.009662, 0.012983, 0.019135, 0.024503, 0.032664, 0.047827, 0.057952, 
                             0.073104, 0.099948, 0.148105, 0.237862, 0.026198, 0.001109, 0.000457, 
                             0.000555, 0.001571, 0.002404, 0.003012, 0.003674, 0.004129, 0.005016, 
                             0.006223, 0.008328, 0.012217, 0.017762, 0.025149, 0.032561, 0.042365, 
                             0.057642, 0.080202, 0.116701, 0.177582, 0.282593, 0.018484, 0.000883, 
                             0.000382, 0.000455, 0.001646, 0.002304, 0.002467, 0.003097, 0.003724, 
                             0.004507, 0.005908, 0.00794, 0.010738, 0.016865, 0.022493, 0.032624, 
                             0.040211, 0.051478, 0.068234, 0.09696, 0.147703, 0.241212, 0.01295, 
                             0.00063, 0.000332, 0.000433, 0.001641, 0.002581, 0.002578, 0.002547, 
                             0.00289, 0.004012, 0.005381, 0.007316, 0.009889, 0.013273, 0.018334, 
                             0.028212, 0.03749, 0.052073, 0.073922, 0.109615, 0.169785, 0.274699,
                             0.045269, 0.002704, 0.000507, 0.00046, 0.000734,  # females
                             0.000895, 0.001126, 0.001495, 0.002197, 0.003143, 0.003983, 0.005939, 
                             0.007469, 0.01166, 0.018486, 0.026548, 0.042649, 0.050858, 0.063509, 
                             0.086965, 0.130587, 0.215029, 0.023838, 0.001154, 0.000358, 0.000318, 
                             0.000502, 0.000698, 0.000918, 0.001144, 0.001572, 0.002207, 0.003151, 
                             0.005038, 0.007183, 0.011023, 0.014718, 0.022267, 0.035953, 0.048153, 
                             0.066424, 0.097196, 0.150869, 0.248412, 0.020248, 0.000933, 0.00031, 
                             0.000339, 0.000525, 0.000652, 0.000901, 0.001251, 0.001599, 0.00223, 
                             0.00313, 0.004514, 0.007125, 0.01058, 0.015764, 0.021294, 0.032344, 
                             0.049166, 0.07543, 0.117877, 0.18764, 0.304247, 0.014603, 0.000768, 
                             0.000271, 0.000287, 0.000487, 0.000565, 0.000715, 0.001059, 0.001481, 
                             0.002049, 0.002936, 0.004201, 0.006039, 0.009984, 0.013853, 0.021179, 
                             0.02809, 0.042159, 0.064247, 0.100939, 0.163497, 0.273028, 0.010488, 
                             0.000521, 0.00025, 0.00029, 0.000453, 0.000581, 0.000725, 0.000901, 
                             0.001171, 0.001816, 0.002734, 0.003782, 0.005293, 0.007575, 0.011174, 
                             0.018559, 0.026524, 0.041711, 0.066135, 0.106604, 0.174691, 0.291021)
                     )

# A to E tests againts spreadsheet -------------------------------------------------------------

# using e_dagger as summary measure
e_dagger <- function(lx){-sum(lx/lx[1]*log(lx/lx[1]))}

# A - test with input nMx, allowing cross-over, and NOT reproducing e0 at given years
outputA <- data.frame(Date = seq(1950,2015,5),
                      Sex = c(rep("m",14),rep("f",14)),
                      e0 = c(55.25,58.41,61.17,63.58,65.68,67.49,69.07,70.46,71.69,72.78,73.75,74.63,75.43,76.15,60.37,63.90,
                            66.97,69.64,71.92,73.89,75.58,77.06,78.35,79.48,80.48,81.36,82.15,82.86)) %>% 
                      dplyr::arrange(Sex, Date)
outputA_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5)) %>% 
                      dplyr::filter(Age == 0) %>% 
                      dplyr::select(Date, Sex, e0=ex) %>% 
                      dplyr::arrange(Sex, Date)
test_that("lc w lim data works", {
  expect_equal(
    as.numeric(outputA_test[,"e0"]), 
    as.numeric(outputA[,"e0"]), 
    tolerance = tolerance_admited
  )
})

# B - test with input nMx, NOT allowing cross-over, and NOT reproducing e0 at given years
outputB <- data.frame(Date = seq(1950,2015,5),
                      Sex = c(rep("m",14),rep("f",14)),
                      e0 = c(54.42,57.48,60.24,62.72,64.95,66.94,69.07,70.46,71.69,72.78,73.75,74.63,75.43,76.15,63.15,65.93,
                             68.39,70.56,72.46,74.13,75.58,77.06,78.35,79.48,80.48,81.36,82.15,82.86))%>% 
                      dplyr::arrange(Sex, Date)
outputB_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5), 
                              prev_divergence = TRUE, error_sheet = TRUE) %>% 
                        dplyr::filter(Age == 0) %>% 
                        dplyr::select(Date, Sex, e0=ex) %>% 
  dplyr::arrange(Sex, Date)

test_that("lc w lim data and prev divergence works", {
  # allow for rounding differences, so maximum absolute difference of 1
  expect_equal(
    as.numeric(outputB_test[,"e0"]), 
    as.numeric(outputB[,"e0"]), 
    tolerance = tolerance_admited
  )
})

# C - test input nMx, allowing cross-over, and reproducing e0 at given years (e0_swe)
data("e0_swe")
outputC <- data.frame(Date = seq(1950,2015,5),
                      Sex = c(rep("m",14),rep("f",14)),
                      e0 = c(70.49,70.99,71.48,71.98,72.21,72.48,73.27,74.30,
                             75.51,76.77,77.89,78.98,79.91,80.73,73.09,74.28,
                             75.47,76.64,77.59,78.39,79.26,80.05,80.91,81.74,
                             82.38,83.10,83.75,84.31)) %>% 
  dplyr::arrange(Sex, Date)
outputC_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                              dates_e0 = unique(e0_swe$Date),
                              e0_Males = e0_swe$e0[e0_swe$Sex=="m"], 
                              e0_Females = e0_swe$e0[e0_swe$Sex=="f"]) %>% 
                          dplyr::filter(Age == 0) %>% 
                          dplyr::select(Date, Sex, e0=ex) %>% 
                          dplyr::arrange(Sex, Date)
test_that("lc w lim data and fitting e0 works", {
  expect_equal(
    as.numeric(outputC_test[,"e0"]), 
    as.numeric(outputC[,"e0"]), 
    tolerance = tolerance_admited
  )
})

# D - test with input nqx, allowing cross-over, and NOT reproducing e0 at given years
input_nqx <- split(input, list(input$Date, input$Sex), drop = F) %>% 
                lapply(function(X){
                        LT = lt_abridged(nMx = X[["nMx"]], 
                                         Age = X[["Age"]],
                                         Sex = unique(X[["Sex"]]))
                        LT$Date = X$Date
                        LT$Sex = X$Sex
                        LT}) %>% 
                do.call("rbind", .) %>% 
  dplyr::select(Date, Sex, Age, nqx)
# paste in spreadsheet: input_nqx %>% tidyr::spread(Date,nqx) %>% write.csv("testD.csv")
outputD <- data.frame(Date = seq(1950,2015,5),
                      Sex = c(rep("m",14),rep("f",14)),
                      e0 = c(55.28,58.43,61.18,63.60,
                             65.69,67.50,69.09,70.47,
                             71.70,72.79,73.77,74.65,
                             75.44,76.17,60.37,63.91,
                             66.98,69.65,71.94,73.91,
                             75.61,77.08,78.37,79.51,
                             80.50,81.39,82.18,82.89
                      )) %>% 
  dplyr::arrange(Sex, Date)
outputD_test <- interp_lc_lim(input = input_nqx, dates_out = seq(1953,2018,5)) %>% 
                    dplyr::filter(Age == 0) %>% 
                    dplyr::select(Date, Sex, e0=ex) %>% 
  dplyr::arrange(Sex, Date)

test_that("lc w lim data and nqx as input works", {
  expect_equal(
    as.numeric(outputD_test[,"e0"]), 
    as.numeric(outputD[,"e0"]), 
    tolerance = tolerance_admited
  )
})
  
# E - test with input lx, allowing cross-over, and NOT reproducing e0 at given years
input_lx <- split(input, list(input$Date, input$Sex), drop = F) %>% 
                    lapply(function(X){
                      LT = lt_abridged(nMx = X[["nMx"]], 
                                       Age = X[["Age"]],
                                       Sex = unique(X[["Sex"]]))
                      LT$Date = X$Date
                      LT$Sex = X$Sex
                      LT}) %>% 
                    do.call("rbind", .) %>% 
  dplyr::select(Date, Sex, Age, lx)
# paste in spreadsheet: input_lx %>% tidyr::spread(Date,lx) %>% xlsx::write.xlsx("testD.xlsx")
outputE <- data.frame(Date = seq(1950,2015,5),
                      Sex = c(rep("m",14),rep("f",14)),
                      e0 = c(55.28,58.43,61.18,63.60,
                             65.69,67.50,69.09,70.47,
                             71.70,72.79,73.77,74.65,
                             75.44,76.17,60.37,63.91,
                             66.98,69.65,71.94,73.91,
                             75.61,77.08,78.37,79.51,
                             80.50,81.39,82.18,82.89)) %>% 
  dplyr::arrange(Sex, Date)
outputE_test <- interp_lc_lim(input = input_lx, dates_out = seq(1953,2018,5)) %>% 
                      dplyr::filter(Age == 0) %>% 
                      dplyr::select(Date, Sex, e0=ex) %>% 
  dplyr::arrange(Sex, Date)
test_that("lc w lim data and nqx as input works", {
  expect_equal(
    as.numeric(outputE_test[,"e0"]), 
    as.numeric(outputE[,"e0"]), 
    tolerance = tolerance_admited
  )
})


# F - testing args ------------------------------------------------------------

# single ages out
outputF1_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                               Single = T)


# single out diff OAG
outputF2_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                               Single = T, extrapLaw = "makeham", OAnew = 100)

# bunch of args
outputF3_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                               prev_divergence = T,
                               dates_e0 = unique(e0_swe$Date),
                               e0_Males = e0_swe$e0[e0_swe$Sex=="m"], 
                               e0_Females = e0_swe$e0[e0_swe$Sex=="f"],
                               Single = T, verbose = F, SVD = T,
                               extrapLaw = "ggompertz", OAnew = 100)

test_that("lc w lim data and nqx as input works", {
  expect_length(unique(outputF1_test$Age), 130)
  expect_s3_class(outputF3_test, "data.frame")
  expect_length(unique(outputF2_test$Age), 101)
  expect_length(unique(outputF3_test$ex), 101 * 2 * length(seq(1953,2018,5)))
}

# G - mixing input --------------------------------------------------------

# some dates gives rates and some lx
input_mix1 <- rbind(input %>% 
                     dplyr::filter(Date %in% dates_in[1:2]) %>% 
                     mutate(lx = NA),
                   input_lx %>% 
                     dplyr::filter(Date %in% dates_in[3:5]) %>% 
                     mutate(nMx = NA)
                   )
outputG1_test <- interp_lc_lim(input = input_mix1, dates_out = seq(1953,2018,5))


# some single and abr ages
input_single <- split(input, list(input$Date, input$Sex), drop = F) %>% 
                  lapply(function(X){
                    LT = lt_abridged2single(nMx = X[["nMx"]], 
                                     Age = X[["Age"]],
                                     Sex = unique(X[["Sex"]]),
                                     OAnew = 100)
                    LT$Date = unique(X$Date)
                    LT$Sex = unique(X$Sex)
                    LT}) %>% 
                  do.call("rbind", .) %>% 
  dplyr::select(Date, Sex, Age, nMx)
input_mix2 <- rbind(input %>% 
                      dplyr::filter(Date %in% dates_in[1:2]),
                    input_single %>% 
                      dplyr::filter(Date %in% dates_in[3:5])
                    )
outputG2_test <- interp_lc_lim(input = input_mix2, dates_out = seq(1953,2018,5))
outputG3_test <- interp_lc_lim(input = input_mix2, dates_out = seq(1953,2018,5),
                               Single = T)


test_that("mixing inputs works", {
expect_s3_class(outputG1_test, "data.frame")
expect_true(all(outputG1_test$nMx > 0))
expect_length(unique(outputG2_test$Age), 22)
expect_length(unique(outputG3_test$Age), 130)
}


# H - lt args -------------------------------------------------------------

# various comb args from lt functions
outputH1_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                              a0rule = "cd")
outputH2_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                              axmethod = "un")
outputH3_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                              region = "n")
outputH4_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                              extrapLaw = "makeham")
outputH5_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                              extrapLaw = "makeham", OAnew = 95)
outputH6_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                               extrapLaw = "gompertz", OAnew = 100, 
                               extrapFrom = 60, extrapFit = seq(50,95,5), radix = 1)

test_that("pass all lt args works", {
expect_s3_class(outputH1_test, "data.frame")
expect_s3_class(outputH2_test, "data.frame")
expect_s3_class(outputH3_test, "data.frame")
expect_s3_class(outputH4_test, "data.frame")
expect_s3_class(outputH5_test, "data.frame")
expect_s3_class(outputH6_test, "data.frame")
}

# I - messages/warnings -------------------------------------------------------

# need n(dates)>2
outputI1_test <- interp_lc_lim(input = input %>% dplyr::filter(Date %in% dates_in[1:2]), 
                               dates_out = seq(1953,2018,5))


# choose e0_dates for you
outputI2_test <- interp_lc_lim(input = input, dates_out = seq(1953,2018,5),
                               # dates_e0 = unique(e0_swe$Date),
                               e0_Males = e0_swe$e0[e0_swe$Sex=="m"], 
                               e0_Females = e0_swe$e0[e0_swe$Sex=="f"])


test_that("mess and warns works", {
  expect_error(outputI1_test)
  expect_error(outputI2_test)
  # tell me you´ll fit with gompertz in case max(Age) is <90
  expect_output(interp_lc_lim(input = input %>% dplyr::filter(Age < 85),
                               dates_out = seq(1953,2018,5)),
                 regexp = "A Gompertz function was fitted for older ages for sex ")
  expect_output(interp_lc_lim(input = input, dates_out = seq(1953,2018,5)),
                regexp = "A Kannisto function was fitted for older ages for sex ")
  
  })


