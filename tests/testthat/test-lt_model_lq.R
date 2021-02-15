# This just tests that the output is a success
check_normal_output <- function(res) {
  expect_true(class(res) == "lt_model_lq")
}

# Generate combination of arguments to test
combn_args <- t(combn(c("q0_5", "q15_45", "e0", "q0_1", "q15_35"), 2))

# Generate combination of values to test
combn_value_lookup <- t(combn(c(0.05, 0.2, 65, 0.05, 0.15), 2))

# Exclude the combinations that will generate an error
error_tests1 <- combn_args[, 1] %in% c("q0_5") & combn_args[, 2] %in% c("q0_1")
error_tests2 <- combn_args[, 1] %in% c("q15_45") & combn_args[, 2] %in% c("q15_35")

# Only test the non-error rows
passing_tests <- !error_tests1 & !error_tests2
rows <- which(passing_tests)
for (i in rows) {
  # Grab the arguments with their respective values
  first_arg <- setNames(combn_value_lookup[i, 1], combn_args[i, 1])
  second_arg <- setNames(combn_value_lookup[i, 2], combn_args[i, 2])

  # Collapse all arguments into a list. b for both,
  # m for males and f for females.
  all_args_b <- c("Sex" = "b", as.list(c(first_arg, second_arg)))
  all_args_m <- c("Sex" = "m", as.list(c(first_arg, second_arg)))
  all_args_f <- c("Sex" = "f", as.list(c(first_arg, second_arg)))

  # just the test name for test_that (easier for debugging)
  test_both <-
    paste0("lt_model_lq with ", names(first_arg), " and ", names(second_arg),
           "works for Sex = 'b'")

  # Same as above but for males
  test_males <-
    paste0("lt_model_lq with ", names(first_arg), " and ", names(second_arg),
           "works for Sex = 'm'")

  # Same as above but for female
  test_females <-
    paste0("lt_model_lq with ", names(first_arg), " and ", names(second_arg),
           "works for Sex = 'f'")

  test_that(test_both, {
    # Run the lt_model_lq with the arguments in a list
    # THIS is where the test happens
    check_normal_output(do.call(lt_model_lq, all_args_b))
  })

  test_that(test_males, {
    check_normal_output(do.call(lt_model_lq, all_args_m))
  })

  test_that(test_females, {
    check_normal_output(do.call(lt_model_lq, all_args_f))
  })
}

# Some combination of arguments need to return an error
rows_error <- which(!passing_tests)
for (i in rows_error) {
  # Grab the arguments with their respective values
  first_arg <- setNames(combn_value_lookup[i, 1], combn_args[i, 1])
  second_arg <- setNames(combn_value_lookup[i, 2], combn_args[i, 2])
  all_args <- c(Sex = "m", as.list(c(first_arg, second_arg)))

  test_name <- paste0("Checking error in lt_model_lq with ",
                      names(first_arg),
                      " and ",
                      names(second_arg))

  test_that(test_name, expect_error(do.call(lt_model_lq, all_args)))
}

