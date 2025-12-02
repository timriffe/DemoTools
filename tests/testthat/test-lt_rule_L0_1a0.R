context("test-lt_rule_L0_1a0")

test_that("lt_rule_L0_1a0 returns expected structure", {
  skip_if_not(exists("lt_rule_ak_m0_a0"))
  
  out <- lt_rule_L0_1a0(L0_target = 0.98, l0 = 1, sex = "f")
  
  expect_type(out, "list")
  expect_named(out, c("a0", "m0", "l0", "l1", "L0_target", "sex"))
  
  expect_true(is.numeric(out$a0))
  expect_true(is.numeric(out$m0))
  expect_true(is.numeric(out$l1))
  expect_equal(out$l0, 1)
  expect_equal(out$L0_target, 0.98)
  expect_equal(out$sex, "f")
})

test_that("lt_rule_L0_1a0 errors when L0_target >= l0", {
  expect_error(
    lt_rule_L0_1a0(L0_target = 1, l0 = 1)
  )
  expect_error(
    lt_rule_L0_1a0(L0_target = 5, l0 = 1)
  )
})

test_that("lt_rule_L0_1a0 errors when l0 inconsistent with inferred radix", {
  # Choose L0_target that implies radix 1
  expect_error(
    lt_rule_L0_1a0(L0_target = 0.8, l0 = 100),
    "Input mismatch",
    fixed = FALSE
  )
})

test_that("lt_rule_L0_1a0 works for both sexes", {
  skip_if_not(exists("lt_rule_ak_m0_a0"))
  
  out_f <- lt_rule_L0_1a0(0.97, l0 = 1, sex = "f")
  out_m <- lt_rule_L0_1a0(0.97, l0 = 1, sex = "m")
  
  expect_true(is.numeric(out_f$a0))
  expect_true(is.numeric(out_m$a0))
  expect_true(out_f$a0 != out_m$a0)  # should differ unless model is symmetric
})

test_that("lt_rule_L0_1a0 produces internally consistent L0 reconstruction", {
  skip_if_not(exists("lt_rule_ak_m0_a0"))
  
  L0_in <- 0.965
  out   <- lt_rule_L0_1a0(L0_target = L0_in, l0 = 1, sex = "f")
  
  L0_reconstructed <- out$l1 + out$a0 * (out$l0 - out$l1)
  expect_equal(round(L0_reconstructed, nchar(L0_in) - 1), L0_in, tolerance = 1e-8)
})

test_that("lt_rule_L0_1a0 handles small L0_target values correctly", {
  skip_if_not(exists("lt_rule_ak_m0_a0"))
  
  # adequately not less than ~0.3
  out <- lt_rule_L0_1a0(0.32, l0 = 1, sex = "f")
  expect_true(out$l1 < 1)
  expect_true(out$m0 > 0)
})

test_that("lt_rule_L0_1a0 fails gracefully with invalid types", {
  expect_error(lt_rule_L0_1a0("a", 1, "f"))
  expect_error(lt_rule_L0_1a0(0.95, "a", "f"))
  expect_error(lt_rule_L0_1a0(0.95, 1, 123))
})

test_that("lt_rule_L0_1a0 solution respects uniroot bounds", {
  skip_if_not(exists("lt_rule_ak_m0_a0"))
  
  out <- lt_rule_L0_1a0(0.9, l0 = 1, sex = "f")
  expect_true(out$l1 >= 0)
  expect_true(out$l1 <= out$l0)
})