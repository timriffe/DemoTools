context("test-movepop")

# create example data for tests
# -------------------------------------------- #
# example for 5 year data
male_pop <- c(48875, 164390, 173551, 130297, 101143, 73615, 60594, 55175,
              49530, 46562, 39028, 27837, 22110, 18066, 15340, 13318,
              12002, 6424)
female_pop <- c(47105, 159546, 168760, 119437, 92080, 70515, 58801, 53381,
                46757, 41164, 33811, 24121, 19315, 16319, 14058, 12302,
                11047, 5922)
male_mx <- c(0.12427, 0.01639, 0.00274, 0.00167, 0.00251, 0.00380, 0.00382,
             0.00442, 0.00506, 0.00663, 0.00872, 0.01240, 0.01783, 0.02700,
             0.04126, 0.06785, 0.11287, 0.21015)
female_mx <- c(0.11050, 0.01577, 0.00254, 0.00159, 0.00232, 0.00304, 0.00344,
               0.00370, 0.00418, 0.00492, 0.00592, 0.00831, 0.01182, 0.01942,
               0.03221, 0.05669, 0.09771, 0.19385)
asfr <- c(0.199, 0.478, 0.418, 0.321, 0.163, 0.071, 0.028)

names(female_mx) <- seq(0, length.out = 18, by = 5)
dxm    <- male_pop * male_mx
dxf    <- female_pop * female_mx
births <- asfr * female_pop[c(4:10)]
# -------------------------------------------- #
# example for single age
dxms <- graduate_pclm(Value  = dxm,
                      Age    = seq(0, length.out = 18, by = 5),
                      OAG    = TRUE)
dxfs <- graduate_pclm(Value  = dxf,
                     Age    = seq(0, length.out = 18, by = 5),
                     OAG    = TRUE)

male_pop_single <- graduate_pclm(Value  = male_pop,
                                 Age    = seq(0, length.out = 18, by = 5),
                                 OAG    = TRUE)
female_pop_single <- graduate_pclm(Value  = female_pop,
                                   Age    = seq(0, length.out = 18, by = 5),
                                   OAG    = TRUE)

male_mx_single <- graduate_pclm(Value  = dxm,
                                Age    = seq(0, length.out = 18, by = 5),
                                OAG    = TRUE,
                                offset = male_pop)
female_mx_single <- graduate_pclm(Value  = dxf,
                                  Age    = seq(0, length.out = 18, by = 5),
                                  OAG    = TRUE,
                                  offset = female_pop)
Age_single <- names(male_pop_single)

asfr_single <- graduate_pclm(Value  = births,
                             Age    = seq(15, length.out = 7, by = 5),
                             OAG    = TRUE,
                             offset = female_pop[c(4:10)],
                             OAnew = max(seq(15, length.out = 7, by = 5)) + 4)

# -------------------------------------------- #
# example for abridged
# pop
male_pop_abr   <- single2abridged(male_pop_single)
female_pop_abr <- single2abridged(female_pop_single)

# deaths
dxmabr <-  single2abridged(dxms)
dxfabr <-  single2abridged(dxfs)

# mx
male_mx_abr   <- dxmabr / male_pop_abr
female_mx_abr <- dxfabr / female_pop_abr
Age_abr       <- names(male_pop_abr)
# Incoming 5-year data test
test_that("movepop handles incoming ages correctly", {
  
  # ---------- projections ----------
  res5_out5 <- movepop(
    initial_date = 1973.58, 
    desired_date = 1973.50,
    male_pop     = male_pop, 
    female_pop   = female_pop,
    male_mx      = male_mx,  
    female_mx    = female_mx,
    asfr         = asfr,
    age_out      = "5-year",
    first_asfr_age = 15,
    annual_net_migrants = -50000
  )
  
  res5_outA <- movepop(
    initial_date = 1973.58, 
    desired_date = 1973.50,
    male_pop     = male_pop, 
    female_pop   = female_pop,
    male_mx      = male_mx,  
    female_mx    = female_mx,
    asfr         = asfr,
    age_out      = "abridged",
    first_asfr_age = 15,
    annual_net_migrants = -50000
  )
  
  res5_outS <- movepop(
    initial_date = 1973.58, 
    desired_date = 1973.50,
    male_pop     = male_pop, 
    female_pop   = female_pop,
    male_mx      = male_mx,  
    female_mx    = female_mx,
    asfr         = asfr,
    age_out      = "single",
    first_asfr_age = 15,
    annual_net_migrants = -50000
  )
  
  # ---------- basic structure tests ----------
  expect_true(is_single(res5_outS$projected_data$age))
  expect_true(is_abridged(res5_outA$projected_data$age))
  expect_true(unique(diff(res5_out5$projected_data$age)) == 5)
  expect_true(min(res5_outS$projected_data$age) == 0)
  expect_true(max(res5_outS$projected_data$age) == 85)
  
})


# Incoming single data
test_that("movepop handles incoming single-year ages correctly", {
  
  # ---- output = single ----
  # first_asfr_age is deduced from vector names
  r_single_charAge <- movepop(
    initial_date = 1973.58,
    desired_date = 1973.50,
    male_pop     = male_pop_single,
    female_pop   = female_pop_single,
    male_mx      = male_mx_single,
    female_mx    = female_mx_single,
    asfr         = asfr_single,
    Age          = Age_single,         # character Ages
    annual_net_migrants = -50000,
    age_out      = "single")
  
  r_single_numAge <- movepop(
    initial_date = 1973.58, 
    desired_date = 1973.50,
    male_pop     = male_pop_single,
    female_pop   = female_pop_single,
    male_mx      = male_mx_single,
    female_mx    = female_mx_single,
    asfr         = asfr_single,
    Age          = as.numeric(Age_single),  # numeric Ages
    annual_net_migrants = -50000,
    age_out      = "single"
  )
  
  # ---- output = 5-year ----
  r_single_out5 <- movepop(
    initial_date = 1973.58, 
    desired_date = 1973.50,
    male_pop     = male_pop_single,
    female_pop   = female_pop_single,
    male_mx      = male_mx_single,
    female_mx    = female_mx_single,
    asfr         = asfr_single,
    Age          = Age_single,
    annual_net_migrants = -50000,
    age_out      = "5-year"
  )
  
  # ---- output = abridged ----
  r_single_outA <- movepop(
    initial_date = 1973.58, 
    desired_date = 1973.50,
    male_pop     = male_pop_single,
    female_pop   = female_pop_single,
    male_mx      = male_mx_single,
    female_mx    = female_mx_single,
    asfr         = asfr_single,
    Age          = Age_single,
    annual_net_migrants = -50000,
    age_out      = "abridged"
  )
  
  
  # ---------- equivalence tests ----------
  expect_equal(r_single_charAge, r_single_numAge)
  
  # ---------- structure checks ----------
  expect_true(is_single(r_single_charAge$projected_data$age))
  expect_true(is_single(r_single_numAge$projected_data$age))
  expect_true(is_abridged(r_single_outA$projected_data$age))
  expect_true(unique(diff(r_single_out5$projected_data$age)) == 5)
  
  expect_true(min(r_single_charAge$projected_data$age) == 0)
  expect_true(max(r_single_charAge$projected_data$age) == 85)
  
})




# Incoming abridged data
test_that("movepop handles incoming asfr correctly", {
  
  names(asfr) <- seq(15, length.out = 7 , by = 5)
  
  # ---- output single ----
  r_abr_to_single <- movepop(
    initial_date = 1973.58, 
    desired_date = 1973.50,
    male_pop     = male_pop_abr,
    female_pop   = female_pop_abr,
    male_mx      = male_mx_abr,
    female_mx    = female_mx_abr,   # NOTE: your female abridged mx had a typo
    asfr         = asfr,
    Age          = NULL,
    annual_net_migrants = -50000,
    age_out      = "single",
    first_asfr_age = NULL # not defined
  )
  
  # ---- output 5-year ----
  r_abr_to_5 <- movepop(
    initial_date = 1973.58, 
    desired_date = 1973.50,
    male_pop     = male_pop_abr,
    female_pop   = female_pop_abr,
    male_mx      = male_mx_abr,
    female_mx    = female_mx_abr,
    asfr         = asfr,
    Age          = Age_abr,
    annual_net_migrants = -50000,
    age_out      = "5-year",
    first_asfr_age = 15 # defined
  )
  
  # ---- output abridged ----
  r_abr_to_abr <- movepop(
    initial_date = 1973.58, 
    desired_date = 1973.50,
    male_pop     = male_pop_abr,
    female_pop   = female_pop_abr,
    male_mx      = male_mx_abr,
    female_mx    = female_mx_abr,
    asfr         = asfr,
    Age          = Age_abr,
    annual_net_migrants = -50000,
    age_out      = "abridged"
  )
  
  # ------------ checks -------------
  expect_true(is_single(r_abr_to_single$projected_data$age))
  expect_true(is_abridged(r_abr_to_abr$projected_data$age))
  expect_true(unique(diff(r_abr_to_5$projected_data$age)) == 5)

  expect_equal(min(r_abr_to_single$projected_data$age), 0)
  expect_equal(max(r_abr_to_single$projected_data$age), 85)
  
})




