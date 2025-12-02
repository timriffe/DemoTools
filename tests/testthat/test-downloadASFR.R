context("test-downloadASFR")

# Mocked ASFR data
mock_asfrmat           <- matrix(runif(5 * 3, 0, 0.1), nrow = 5, ncol = 3)
rownames(mock_asfrmat) <- c(15, 20, 25, 30, 35)
colnames(mock_asfrmat) <- c(2000, 2005, 2010)

# Test 1: User-supplied matrix is returned unchanged
test_that("user-supplied Asfrmat is returned correctly", {
  res <- downloadASFR(
    Asfrmat     = mock_asfrmat,
    Age         = c(15, 20, 25, 30, 35),
    AsfrDatesIn = c(2000, 2005, 2010),
    output      = "single"
  )
  expect_true(is.matrix(res))
  expect_equal(dim(res), dim(mock_asfrmat))
  expect_equal(rownames(res), c("15", "20", "25", "30", "35"))
  expect_equal(colnames(res), c("2000", "2005", "2010"))
})

# Test 2: Error when Asfrmat and Age length mismatch
test_that("error when Asfrmat rows and Age length mismatch", {
  expect_error(
    downloadASFR(Asfrmat = mock_asfrmat, 
                 Age     = 1:4),
    "Inconsistent input"
  )
})

# Test 3: Error when Asfrmat cols and AsfrDatesIn length mismatch
test_that("error when Asfrmat cols and AsfrDatesIn length mismatch", {
  expect_error(
    downloadASFR(Asfrmat     = mock_asfrmat, 
                 Age         = 15:19, 
                 AsfrDatesIn = 2000:2004),
    "Inconsistent input"
  )
})

# Test 4: Error when location is missing and Asfrmat is NULL
test_that("error when location is NULL and Asfrmat is NULL", {
  expect_error(
    downloadASFR(Asfrmat     = NULL, 
                 location    = NULL, 
                 AsfrDatesIn = 2000:2005),
    "You need to provide a location"
  )
})

# Test: 5-year ASFR output matrix has correct dimensions
test_that("5-year ASFR output matrix has correct dimensions", {
  
  # Run the function
  res <- downloadASFR(
    Asfrmat     = NULL,
    location    = "Argentina",
    AsfrDatesIn = 2000:2001,
    output      = "5-year",
    method      = "linear"
  )
  
  # Assertions
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 9)  # ages 10, 15,20,25,30,35,40,45, 50 assuming 5-year bins
  expect_equal(ncol(res), 2)  # years 2000, 2001
})


# Test 6: ASFR download with single-year output matrix
test_that("single-year ASFR output matrix has correct dimensions", {
  # reuse previous mocks
  res <- downloadASFR(
    Asfrmat     = NULL,
    location    = "Argentina",
    AsfrDatesIn = 2000:2001,
    output      = "single",
    method      = "linear"
  )
  
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 45)  # ages 10:54
  expect_equal(ncol(res), 2)   # years 2000, 2001
})

# Test 7: ASFR extrapolation beyond WPP range triggers message
test_that("message is printed for years outside WPP range", {
  expect_message(
    downloadASFR(
      Asfrmat     = NULL,
      location    = "Argentina",
      AsfrDatesIn = 1900:1902,
      output      = "single"
    ),
    "Careful, choosing beyond range"
  )
})
