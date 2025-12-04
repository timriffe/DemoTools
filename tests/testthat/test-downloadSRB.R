context("test-downloadSRB")

# Helper to mock installed.packages() matrix
mk_installed_mat <- function(pkgs) {
  
  if(length(pkgs) == 0) {
    
    mat <- matrix(nrow = 0, 
                  ncol = 1)
    
    colnames(mat) <- "Package"
    
    return(mat)
    
  }
  
  mat <- matrix("x", 
                nrow     = length(pkgs), 
                ncol     = 1,
                dimnames = list(pkgs, "Package"))
  return(mat)
}

# Helper to clean up objects assigned in tests
teardown_assign_cleanup <- function(names_vec) {
  
  for(nm in names_vec) {
    
    if(exists(nm, envir = .GlobalEnv)) rm(list = nm, envir = .GlobalEnv)
    
  }
  
}

# Imaginable datasets for tests
sexRatio_5yr <- tibble(
  country_code = c(900L, 32L, 724L),
  name         = c("World", "Argentina", "Spain"),
  `1950-1955`  = c(1.05, 1.05, 1.04),
  `1955-1960`  = c(1.05, 1.05, 1.04)
)

sexRatio1_single <- tibble(
  country_code = c(900L, 32L, 724L),
  name         = c("World", "Argentina", "Spain"),
  `1950`       = c(1.05, 1.05, 1.04),
  `1951`       = c(1.05, 1.05, 1.04),
  `1952`       = c(1.05, 1.06, 1.04)
)

# UNIT TESTS 
test_that("user-supplied SRB returns repeated named vector", {
  res <- downloadSRB(SRB      = 1.05,
                     DatesOut = 1950:1952,
                     verbose  = FALSE)
  expect_named(res, as.character(1950:1952))
  expect_equal(unname(res), rep(1.05, 3))
})

test_that("DatesOut length < 1 triggers error", {
  expect_error(downloadSRB(SRB       = NULL, 
                            location = "Argentina", 
                            DatesOut = NULL),
               "DatesOut must contain at least one date")
})

test_that("no installed WPP triggers a clear error", {
  
  # mock no installed wpp
  installed_mat <- mk_installed_mat(character(0))
  mock_inst     <- mock(installed_mat)
  
  expect_error(
    with_mock(
      `installed.packages`  = mock_inst,
      downloadSRB(SRB       = NULL, 
                   location = "Argentina", 
                   DatesOut = 1950),
      "No WPP package installed"
    )
  )
})

test_that("WPP >= 2022 branch uses sexRatio1 and returns correct SRB", {
  
  
  res <- downloadSRB(SRB       = NULL, 
                     location = "Argentina", 
                     DatesOut = 1950:1952, 
                     verbose  = FALSE)
  
  assign("sexRatio1", sexRatio1_single, envir = .GlobalEnv)
  on.exit(teardown_assign_cleanup(c("sexRatio1")), add = TRUE)
  
  expect_named(res, as.character(1950:1952))
  expect_equal(unname(res), c(1.047, 1.046, 1.047))
})


test_that("missing location returns default SRB", {
  
  # Capture messages
  msgs <- capture_messages(
    res <- downloadSRB(SRB       = NULL,
                       location = "Hyperborea",
                       DatesOut = 1950:1952,
                       verbose  = TRUE)
  )
  
  # Check SRB values
  expect_equal(unname(res), rep(1.047, 3))
  
  # Check that warning/message is printed
  expect_true(any(grepl("not available in wpp. Using default SRB", msgs)))
})

test_that("numeric location code works", {
  
  res <- downloadSRB(SRB       = NULL, 
                      location = "Argentina", 
                      DatesOut = 1950:1951, 
                      verbose  = FALSE)
  
  
  res1 <- downloadSRB(SRB       = NULL,
                       location = 32,
                       DatesOut = 1950:1951, 
                       verbose  = FALSE)
  expect_equal(unname(res), unname(res1))
  
})
