# Author: ilya
###############################################################################
context("AGESEX")

# data --------------------------------------------------------------------

# PAS spreadsheet AGESEX.xlsx
# https://www.census.gov/data/software/pas.html [we may think of a proper cit]

# data from PAS spreadsheet AGESEX.xlsx
Males   <- c(4677000, 4135000, 3825000, 3647000, 3247000, 
             2802000, 2409000, 2212000, 1786000, 1505000, 
             1390000, 984000, 745000, 537000, 346000, 335000)
Females <- c(4544000, 4042000, 3735000, 3647000, 3309000, 
             2793000, 2353000, 2112000, 1691000, 1409000, 
             1241000, 887000, 697000, 525000, 348000, 366000)
Age     <- seq(0, 75, by = 5)


# tests -------------------------------------------------------------------

test_that("ageRatioScore works", {
    
    expect_equal(
        ageRatioScore(
            Value = Males, 
            Age = Age, 
            ageMax = 75
        ), 
        3.907, 
        tolerance = .001
    )
    
    expect_equal(
        ageRatioScore(
            Value = Females, 
            Age = Age, 
            ageMax = 75
        ), 
        3.655, 
        tolerance = .001
    )
})


test_that("sexRatioScore works", {
    
    expect_equal(
        sexRatioScore(
            Males = Males, 
            Females = Females, 
            Age = Age
        ), 
        2.249, 
        tolerance = .001
    )
    
    expect_equal(
        sexRatioScore(
            Males = Males, 
            Females = Females, 
            Age = Age, 
            ageMax = 70
        ), 
        2.249, 
        tolerance = .001
    )
})


test_that("ageSexAccuracy works", {
    
    expect_equal(
        ageSexAccuracy(
            Males = Males, 
            Females = Females, 
            Age = Age
        ), 
        14.308, 
        tolerance = .001
    )
})