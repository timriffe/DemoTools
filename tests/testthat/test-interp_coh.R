check_form <- function(x) {
  expect_is(x, "matrix")
  expect_true(nrow(x) == 101)
  expect_true(ncol(x) == 8)
  expect_true(length(colnames(x)) != 0)
}

births <- c(719511L, 760934L, 772973L, 749554L, 760831L,
            828772L, 880543L, 905380L, 919639L)

test_that("interp_coh works without midyear", {

  res <-
    interp_coh(
    country = "Russian Federation",
    sex = "male",
    c1 = pop1m_rus2002,
    c2 = pop1m_rus2010,
    date1 = "2002-10-16",
    date2 = "2010-10-25",
    age1 = 0:100,
    births = births
    )

  check_form(res)
})

test_that("interp_coh works without midyear", {

  res <-
    interp_coh(
      country = "Russian Federation",
      sex = "male",
      c1 = pop1m_rus2002,
      c2 = pop1m_rus2010,
      date1 = "2002-10-16",
      date2 = "2010-10-25",
      age1 = 0:100,
      births = births,
      midyear = TRUE
    )

  check_form(res)
})
