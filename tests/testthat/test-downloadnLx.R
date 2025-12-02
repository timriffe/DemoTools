context("test-downloadnLx")

test_that("User-supplied nLx is returned unchanged (with Age inferred)", {
  
  nLx_data <- matrix(runif(20), nrow = 4)
  colnames(nLx_data) <- c(2000, 2005, 2010, 2015, 2020)
  
  res <- downloadnLx(nLx       = nLx_data,
                     nLxDatesIn = c(2000, 2005, 2010, 2015, 2020),
                     output     = "5-year")
  
  expect_equal(res, nLx_data)
})

test_that("User-supplied nLx with mismatched Age throws error", {

  m <- matrix(1:10, nrow = 5)
  
  expect_error(
    downloadnLx(nLx = m, 
                Age = 1:3),
    "Inconsistent input"
  )
})

test_that("User-supplied nLx with mismatched dates throws error", {
  m <- matrix(1:10, ncol = 5)
  
  expect_error(
    downloadnLx(nLx        = m, 
                nLxDatesIn = 2000:2003),
    "Inconsistent input"
  )
})

# -------------------------------------------------------------------
# MOCKING WPP INSTALLATION
# -------------------------------------------------------------------

test_that("Function errors when no WPP is installed", {
  stub(downloadnLx, "installed.packages", function() matrix("", nrow = 0, ncol = 1))
  expect_error(
    downloadnLx(location = "Argentina", gender = "female", nLxDatesIn = 2000),
    "No wpp package installed"
  )
})

# -------------------------------------------------------------------
# MOCK WPP ≥ 2022 PATH
# -------------------------------------------------------------------

test_that("WPP >= 2022 path returns matrix of correct dimensions", {
  
  res <- downloadnLx(
    location   = "Argentina",
    gender     = "female",
    nLxDatesIn = c(2000:2001),
    output     = "single"
  )
  
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 101)
  expect_equal(ncol(res), 2)
})
# -------------------------------------------------------------------
# OUTPUT TYPE TESTS
# -------------------------------------------------------------------

test_that("Output type = 'abridged' produces collapsed age groups", {
  
  res <- downloadnLx(
    location   = "Argentina",
    gender     = "female",
    nLxDatesIn = 2000,
    output     = "abridged"
  )
  
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 22)
  expect_equal(ncol(res), 1)
})

test_that("Output type = '5-year' produces collapsed age groups", {
  
  res <- downloadnLx(
    location   = "Argentina",
    gender     = "female",
    nLxDatesIn = c(2000, 2001),
    output     = "5-year"
  )
  
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 21)
  expect_equal(ncol(res), 2)
})


# note the older package version mock
test_that("Output type = 'abridged' produces collapsed age groups", {
  
  # Make the function believe only wpp2019 is installed  
  fake_installed <- matrix(
    "wpp2019",
    nrow = 1,
    dimnames = list("wpp2019", NULL)
  )
  
  stub(downloadnLx, "installed.packages", fake_installed)
  
  # Run as normal
  expect_warning(
  res <- downloadnLx(
    location   = "Argentina",
    gender     = "female",
    nLxDatesIn = 2000 ,
    output     = "5-year"
  )
  )
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 21)
  expect_equal(ncol(res), 1)
  
})


test_that("Output type = '5-year' produces collapsed age groups", {
  
  # Make the function believe only wpp2019 is installed  
  fake_installed <- matrix(
    "wpp2019",
    nrow = 1,
    dimnames = list("wpp2019", NULL)
  )
  
  stub(downloadnLx, "installed.packages", fake_installed)
  
  # Run as normal
  expect_warning(
  res <- downloadnLx(
    location   = "Argentina",
    gender     = "female",
    nLxDatesIn = c(2000, 2006),
    output     = "5-year"
  )
  )
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 21)
  expect_equal(ncol(res), 2)
  
})



test_that("Output type = 'single' produces collapsed age groups", {
  
  # Make the function believe only wpp2019 is installed  
  fake_installed <- matrix(
    "wpp2019",
    nrow = 1,
    dimnames = list("wpp2019", NULL)
  )
  
  stub(downloadnLx, "installed.packages", fake_installed)
  
  # Run as normal
  expect_warning(
  res <- downloadnLx(
    location   = "Argentina",
    gender     = "female",
    nLxDatesIn = c(2000, 2006, 2010),
    output     = "single"
  )
  )
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 101)
  expect_equal(ncol(res), 3)
  
})


test_that("downloadnLx errors when WPP < 2022 and output='single'", {
  
  fake_installed <- matrix("wpp2019", nrow = 1, 
                           dimnames = list("wpp2019",NULL))
  stub(downloadnLx, "installed.packages", fake_installed)
  
  expect_warning(
    downloadnLx(
      nLx        = NULL,
      Age        = Age_five,
      location   = location,
      gender     = "female",
      nLxDatesIn = refDate - c(0.5, 7.5),
      output     = "single",   # not allowed in WPP < 2022
      radix      = 1
    ),
    regexp = "No single ages are available in wpp versions"
  )
})




