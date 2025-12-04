context("test-lt_rule_l0_1L0")

test_that("lt_rule_l0_1L0 returns radix correctly for powers of 10", {
  expect_equal(lt_rule_l0_1L0(1e1), 10)
  expect_equal(lt_rule_l0_1L0(1e2), 100)
  expect_equal(lt_rule_l0_1L0(1e4), 10000)
  expect_equal(lt_rule_l0_1L0(1e6), 1000000)
})

test_that("lt_rule_l0_1L0 infers next higher power of 10 for non-radix inputs", {
  expect_equal(lt_rule_l0_1L0(20), 100)
  expect_equal(lt_rule_l0_1L0(999), 1000)
  expect_equal(lt_rule_l0_1L0(10001), 100000)
})

test_that("lt_rule_l0_1L0 handles decimal/non-integer values", {
  expect_equal(lt_rule_l0_1L0(20.3), 100)
  expect_equal(lt_rule_l0_1L0(0.5), 1)
  expect_equal(lt_rule_l0_1L0(2), 10)
})

test_that("lt_rule_l0_1L0 handles boundary values", {
  expect_equal(lt_rule_l0_1L0(1), 1)
  expect_equal(lt_rule_l0_1L0(0), 1)
  expect_equal(lt_rule_l0_1L0(1.000000000), 1)
})

test_that("lt_rule_l0_1L0 handles NA, NaN, Inf", {
  expect_true(is.na(lt_rule_l0_1L0(NA_real_)))
  expect_true(is.na(lt_rule_l0_1L0(NaN)))
  expect_equal(lt_rule_l0_1L0(Inf), NA_real_)
})

test_that("lt_rule_l0_1L0 does not support vector input", {
  expect_error(lt_rule_l0_1L0(c(10, 20)))
})
