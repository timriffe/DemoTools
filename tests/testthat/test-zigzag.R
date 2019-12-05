# Author: ilya
###############################################################################
context("test-zigzag")

# Feeney, G. 2013 "Removing "Zigzag" from Age Data," http://demographer.com/white-papers/2013-removing-zigzag-from-age-data/

test_that("smooth_age_5_zigzag works",{
    
    Age <- c(0, 1, seq(5, 90, by = 5))

    # check result using results frozen in Feeney spreadsheet
    # after fixing probable cell range 'error'
    p.feeney <- c(0.0235802695087692, 0.0724286618207911,
                  0.0242327829742267, 0.0883411499065237)
    
    ans      <- 106.1147411629

    
    expect_equal(
        smooth_age_5_zigzag_min(dth5_zigzag, Age, 40,80,p.feeney),
        ans,
        tolerance = 1e-6
    )
})

